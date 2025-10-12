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

function benders_RO(instance, options; silent = true, force=1.0, Γ::Int64=0, gap=0.05, timelimit=10)

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
            if second_stage_cost_ub<= 0.999999*UB
                UB=second_stage_cost_ub
                cstr=@build_constraint(master_pb[:obj] <= (1+gap)*second_stage_cost_ub)
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
            if second_stage_cost_ub<= 0.999999*UB
                UB=second_stage_cost_ub
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