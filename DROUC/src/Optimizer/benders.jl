function return_solution(master_pb, computation_time, instance, S, batch, Time_subproblem, gap, force)
    feasibleSolutionFound = primal_status(master_pb) == MOI.FEASIBLE_POINT

    if feasibleSolutionFound

        solution_is_on = round.(convert(Matrix{Float64}, value.(master_pb[:is_on])))
        solution_start_up = round.(convert(Matrix{Float64}, value.(master_pb[:start_up])))
        solution_start_down = round.(convert(Matrix{Float64}, value.(master_pb[:start_down])))
        solution_power_ref=JuMP.value.(master_pb[:power_first_stage])

        first_stage_solution = first_stage(solution_power_ref, solution_is_on, solution_start_up, solution_start_down)

        return instance.name, computation_time, objective_value(master_pb), objective_bound(master_pb), first_stage_solution, S,batch, Time_subproblem, gap, force, value.(master_pb[:prod_tot_first_stage])
    else
        return Time_subproblem
    end
end

function CCG_RO(instance, options; silent = true, force=1.0, Γ::Int64=0, gap=0.05, timelimit=10, iter=5)

    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]
    Lines=values(instance.Lines)
    thermal_units=values(instance.Thermalunits)
    thermal_units_2=values(instance.Thermalunits)[N1+1:N]
    Demandbus=instance.Demandbus
    uncertainty_0 = instance.uncertainty

    master_pb=options.master_problem(instance, silent=silent)
    oracle_pb=options.oracle_problem(instance, Γ=Γ)
    newgap=gap/100
    set_optimizer_attribute(master_pb, "MIPGap", newgap/3)
    set_optimizer_attribute(oracle_pb, "MIPGap", newgap)
    is_on=master_pb[:is_on]

    infty=1e9
    start = time()

    Time_subproblem=[]
    Time_master=[]

    LB=1
    UB=infty
    k=0

    while time() - start<=timelimit && k<=iter && 100*(UB-LB)/UB>=gap
        k+=1
        set_optimizer_attribute(master_pb, "TimeLimit", timelimit-(time() - start))
        set_optimizer_attribute(oracle_pb, "TimeLimit", min(20,timelimit-(time() - start)))
        start_master = time()
        optimize!(master_pb)
        push!(Time_master, time() - start_master)
        LB=max(JuMP.objective_bound(master_pb), LB)

        solution_is_on = convert(Matrix{Float64}, round.(value.(master_pb[:is_on])))
        solution_start_up = convert(Matrix{Float64}, round.(value.(master_pb[:start_up])))
        solution_start_down = convert(Matrix{Float64}, round.(value.(master_pb[:start_down])))

        solution_x = [solution_is_on, solution_start_up, solution_start_down]

        prod_tot_first_stage_val=value.(master_pb[:prod_tot_first_stage])

        value_subproblem, d_up, d_down, computation_time=options.second_stage(instance, options, oracle_pb, prod_tot_first_stage_val, solution_x; Γ=Γ, force=force)
        push!(Time_subproblem, computation_time)

        # println((UB, JuMP.objective_value(master_pb)+value_subproblem-JuMP.value.(master_pb[:thermal_cost])))

        UB=min(UB, JuMP.objective_value(master_pb)+value_subproblem-JuMP.value.(master_pb[:thermal_cost]))

        # println((k, LB, UB, 100*(UB-LB)/UB))

        if value_subproblem>=JuMP.value.(master_pb[:thermal_cost])+0.1 && k <= iter && 100*(UB-LB)/UB>=gap && time() - start<=timelimit-0.1

            unregister(master_pb, :θ)
            unregister(master_pb, :flow)
            @variable(master_pb, θ[b in Buses, t in 1:T, l in k:k])
            @variable(master_pb, flow[b in Buses, bp in Next[b], t in 1:T, l in k:k])

            @constraint(master_pb, [line in Lines, t in 1:T], flow[line.b1,line.b2,t, k]<=line.Fmax)
            @constraint(master_pb, [line in Lines, t in 1:T], flow[line.b1,line.b2,t,k]>=-line.Fmax)
            @constraint(master_pb, [line in Lines, t in 1:T], flow[line.b1,line.b2,t,k]==line.B12*(θ[line.b1,t, k]-θ[line.b2,t,k]))
            
            unregister(master_pb, :power)
            unregister(master_pb, :power_shedding)
            unregister(master_pb, :power_curtailement)
            @variable(master_pb, power[i in 1:N, t in 0:T, l in k:k]>=0)
            @variable(master_pb, power_shedding[b in Buses, t in 0:T, l in k:k]>=0)
            @variable(master_pb, power_curtailement[b in Buses, t in 0:T, l in k:k]>=0)

            @constraint(master_pb, [i in 1:N1, t in 1:T], power[i, t, k] == master_pb[:power_first_stage][i, t])
            
            @constraint(master_pb,  [unit in thermal_units, t in 0:T], power[unit.name, t, k]>=unit.MinPower*is_on[unit.name, t])
            @constraint(master_pb,  [unit in thermal_units, t in 0:T], power[unit.name, t, k]<=unit.MaxPower*is_on[unit.name, t])
            @constraint(master_pb,  [unit in thermal_units], power[unit.name, 0, k]==unit.InitialPower)

            d0 = [(b,t) for b in Buses for t in 1:T if (b,t) ∉ d_up && (b,t) ∉ d_down]

            @constraint(master_pb,  [(b,t) in d0], sum(power[unit.name, t,k] for unit in thermal_units if unit.Bus==b)+power_shedding[b,t,k]-power_curtailement[b,t,k]==instance.Demandbus[b][t]+sum(flow[b,bp,t,k] for bp in Next[b]))

            @constraint(master_pb,  [(b,t) in d_up], sum(power[unit.name, t,k] for unit in thermal_units if unit.Bus==b)+power_shedding[b,t,k]-power_curtailement[b,t,k]==instance.Demandbus[b][t]+sum(flow[b,bp,t,k] for bp in Next[b])+force*1.96*uncertainty_0.dev[t]*Demandbus[b][t])
            @constraint(master_pb,  [(b,t) in d_down], sum(power[unit.name, t,k] for unit in thermal_units if unit.Bus==b)+power_shedding[b,t,k]-power_curtailement[b,t,k]==instance.Demandbus[b][t]+sum(flow[b,bp,t,k] for bp in Next[b])-force*1.96*uncertainty_0.dev[t]*Demandbus[b][t])

            @constraint(master_pb, master_pb[:thermal_cost]>= sum(unit.LinearTerm*power[unit.name, t,k] for unit in thermal_units_2 for t in 1:T) + sum(SHEDDING_COST*power_shedding[b,t,k]+CURTAILEMENT_COST*power_curtailement[b,t,k] for b in Buses for t in 1:T))

        else
            computation_time = time() - start

            solution_is_on = round.(convert(Matrix{Float64}, value.(master_pb[:is_on])))
            solution_start_up = round.(convert(Matrix{Float64}, value.(master_pb[:start_up])))
            solution_start_down = round.(convert(Matrix{Float64}, value.(master_pb[:start_down])))
            solution_power_ref=JuMP.value.(master_pb[:power_first_stage])

            first_stage_solution = first_stage(solution_power_ref, solution_is_on, solution_start_up, solution_start_down)
            
            return instance.name, computation_time, objective_value(master_pb), objective_bound(master_pb), first_stage_solution, Time_subproblem, Time_master
        end
    end

    computation_time = time() - start

    solution_is_on = round.(convert(Matrix{Float64}, value.(master_pb[:is_on])))
    solution_start_up = round.(convert(Matrix{Float64}, value.(master_pb[:start_up])))
    solution_start_down = round.(convert(Matrix{Float64}, value.(master_pb[:start_down])))
    solution_power_ref=JuMP.value.(master_pb[:power_first_stage])

    first_stage_solution = first_stage(solution_power_ref, solution_is_on, solution_start_up, solution_start_down)
        
    return instance.name, computation_time, objective_value(master_pb), objective_bound(master_pb), first_stage_solution, Time_subproblem, Time_master
end

function benders_RO(instance, options; silent = true, force=1.0, Γ::Int64=0, gap=0.05, timelimit=10, iter=5)

    master_pb=options.master_problem(instance, silent=silent)
    oracle_pb=options.oracle_problem(instance, Γ=Γ)

    infty=1e9
    LB=1
    start = time()

    Time_subproblem=[]

    LB=1
    UB=infty
    k=0

    function my_callback_function(cb_data, cb_where::Cint)

        if cb_where==GRB_CB_MIPSOL
            k+=1
            Gurobi.load_callback_variable_primal(cb_data, cb_where)
            solution_is_on = convert(Matrix{Float64}, round.(callback_value.(cb_data, master_pb[:is_on])))
            solution_start_up = convert(Matrix{Float64}, round.(callback_value.(cb_data, master_pb[:start_up])))
            solution_start_down = convert(Matrix{Float64}, round.(callback_value.(cb_data, master_pb[:start_down])))

            solution_x = [solution_is_on, solution_start_up, solution_start_down]

            second_stage_cost_ub, thermal_cost_scenario = options.add_cut(cb_data, options, master_pb, oracle_pb, instance, Time_subproblem, solution_x; force=force, gap=gap, Γ=Γ)
            # println((callback_value(cb_data, master_pb[:thermal_fixed_cost]), second_stage_cost_ub, UB, 0.999999*UB))
            if second_stage_cost_ub<= 0.999999*UB
                UB=second_stage_cost_ub
                cstr=@build_constraint(master_pb[:obj] <= (1+gap)*second_stage_cost_ub)
                # MOI.submit(master_pb, MOI.LazyConstraint(cb_data), cstr)
                if k>=10
                    vars = options.get_variables(instance, master_pb)
                    vals = options.get_value_variables(instance, master_pb, cb_data, solution_is_on, solution_start_up, solution_start_down, thermal_cost_scenario)

                    MOI.submit(master_pb, MOI.HeuristicSolution(cb_data), vars, vals)
                end
            end
        end
        return
    end

    set_optimizer_attribute(master_pb, "TimeLimit", timelimit-(time() - start))
    set_optimizer_attribute(master_pb, "LazyConstraints", 1)
    MOI.set(master_pb, Gurobi.CallbackFunction(), my_callback_function)
    start = time()
    newgap=gap/100
    set_optimizer_attribute(master_pb, "MIPGap", newgap)
    set_optimizer_attribute(oracle_pb, "MIPGap", newgap)
    optimize!(master_pb)
    computation_time = time() - start


    solution_is_on = round.(convert(Matrix{Float64}, value.(master_pb[:is_on])))
    solution_start_up = round.(convert(Matrix{Float64}, value.(master_pb[:start_up])))
    solution_start_down = round.(convert(Matrix{Float64}, value.(master_pb[:start_down])))
    solution_power_ref=JuMP.value.(master_pb[:power_first_stage])

    first_stage_solution = first_stage(solution_power_ref, solution_is_on, solution_start_up, solution_start_down)

    return instance.name, computation_time, objective_value(master_pb), objective_bound(master_pb), first_stage_solution, Time_subproblem
end

function benders(instance, options; ρ=0, silent=true, force=1, S::Int64, batch=1, gap=0.05, timelimit=10)

    radius = options.compute_radius(S, ρ)
    master_pb=options.master_problem(instance, silent=silent, ρ=radius, S=S)
    oracle_pb=options.oracle_problem(instance)

    infty=1e9
    LB=1
    start = time()

    Time_subproblem=[]

    LB=1
    UB=infty
    k=0

    function my_callback_function(cb_data, cb_where::Cint)


        if options == DCA_l2 || options == Grb_l2
            if cb_where==GRB_CB_MIPNODE
                resultP = Ref{Cint}()
                GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
                if resultP[] != GRB_OPTIMAL
                    return
                end
                Gurobi.load_callback_variable_primal(cb_data, cb_where)
                λ_val=callback_value(cb_data, master_pb[:λ])
                w_val=callback_value(cb_data, master_pb[:w])
                if abs(w_val - 1/λ_val)>=1e-5 && λ_val >=1e-6
                    cstr=@build_constraint(master_pb[:w]>=1/λ_val-1/λ_val^2*(master_pb[:λ]-λ_val))
                    MOI.submit(master_pb, MOI.LazyConstraint(cb_data), cstr)
                end
            end

            if cb_where==GRB_CB_MIPSOL
                Gurobi.load_callback_variable_primal(cb_data, cb_where)
                λ_val=callback_value(cb_data, master_pb[:λ])
                w_val=callback_value(cb_data, master_pb[:w])
                if abs(w_val - 1/λ_val)>=1e-5 && λ_val >=1e-6
                    cstr=@build_constraint(master_pb[:w]>=1/λ_val-1/λ_val^2*(master_pb[:λ]-λ_val))
                    MOI.submit(master_pb, MOI.LazyConstraint(cb_data), cstr)
                end
            end
        end

        if cb_where==GRB_CB_MIPSOL
            k+=1
            Gurobi.load_callback_variable_primal(cb_data, cb_where)
            solution_is_on = convert(Matrix{Float64}, round.(callback_value.(cb_data, master_pb[:is_on])))
            solution_start_up = convert(Matrix{Float64}, round.(callback_value.(cb_data, master_pb[:start_up])))
            solution_start_down = convert(Matrix{Float64}, round.(callback_value.(cb_data, master_pb[:start_down])))

            solution_x = [solution_is_on, solution_start_up, solution_start_down]

            second_stage_cost_ub, thermal_cost_scenario = options.add_cut(cb_data, options, master_pb, oracle_pb, instance, Time_subproblem, solution_x; force=force, S=S, batch=batch, gap=gap, ρ=ρ)
            # println((callback_value(cb_data, master_pb[:thermal_fixed_cost]), second_stage_cost_ub, UB, 0.999999*UB))
            if second_stage_cost_ub<= 0.999999*UB
                UB=second_stage_cost_ub
                # println((UB, callback_value(cb_data, master_pb[:λ])))
                # cstr=@build_constraint(master_pb[:obj] <= (1+gap)*second_stage_cost_ub)
                # MOI.submit(master_pb, MOI.LazyConstraint(cb_data), cstr)
                if k>=10
                    vars = options.get_variables(instance, master_pb, S)
                    vals = options.get_value_variables(instance, master_pb, S, cb_data, solution_is_on, solution_start_up, solution_start_down, thermal_cost_scenario)

                    MOI.submit(master_pb, MOI.HeuristicSolution(cb_data), vars, vals)
                end
            end
        end
        return
    end

    set_optimizer_attribute(master_pb, "TimeLimit", timelimit-(time() - start))
    set_optimizer_attribute(master_pb, "LazyConstraints", 1)
    MOI.set(master_pb, Gurobi.CallbackFunction(), my_callback_function)
    start = time()
    newgap=gap/100
    set_optimizer_attribute(master_pb, "MIPGap", newgap)
    set_optimizer_attribute(oracle_pb, "MIPGap", newgap)
    optimize!(master_pb)
    computation_time = time() - start

    return return_solution(master_pb, computation_time, instance, S, batch, Time_subproblem, gap, force)
end