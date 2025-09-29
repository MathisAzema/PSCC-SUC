function create_master_problem(instance; silent=true, S=1)
    """
    Create the master problem for the Benders' decomposition.
    """
    model = Model(instance.optimizer)
    if silent
        set_silent(model)
    end

    T= instance.TimeHorizon
    thermal_units_name=keys(instance.Thermalunits)
    thermal_units=values(instance.Thermalunits)

    Next=instance.Next
    Buses=1:size(Next)[1]
    @variable(model, power_ref[unit in thermal_units_name, t in 0:T]>=0)
    @variable(model, reserve_up[unit in thermal_units_name, t in 0:T]>=0)
    @variable(model, reserve_down[unit in thermal_units_name, t in 0:T]>=0)
    @variable(model, is_on[unit in thermal_units_name, t in 0:T], Bin)
    @variable(model, start_up[unit in thermal_units_name, t in 1:T], Bin)
    @variable(model, start_down[unit in thermal_units_name, t in 1:T], Bin)

    # @constraint(model, [t in 0:T], is_on[1, t] ==1)
    # @constraint(model, [t in 2:T], is_on[2, t] ==0)
    # @constraint(model, [t in 1:T], is_on[3, t] ==1)
    # @constraint(model, [t in 0:3], is_on[3, t] ==0)

    # @variable(model, thermal_fuel_cost[s in 1:S]>=0)
    @variable(model, thermal_fixed_cost>=0)
    @variable(model, thermal_cost>=0)

    # @constraint(model,  thermal_cost>=sum(thermal_fuel_cost[s] for s in 1:S)/S)

    thermal_unit_commit_constraints2(model, instance)

    @constraint(model,  [unit in thermal_units, t in 0:T], power_ref[unit.name, t]-reserve_down[unit.name, t]>=unit.MinPower*is_on[unit.name, t])
    @constraint(model,  [unit in thermal_units, t in 0:T], power_ref[unit.name, t]+reserve_up[unit.name, t]<=unit.MaxPower*is_on[unit.name, t])
    @constraint(model,  [unit in thermal_units, s in 1:S; unit.InitialPower!=nothing], power_ref[unit.name, 0]==unit.InitialPower)
    @constraint(model,  [unit in thermal_units], reserve_up[unit.name,0]==0)
    @constraint(model,  [unit in thermal_units], reserve_down[unit.name,0]==0)

    @constraint(model,  [unit in thermal_units, t in 1:T], power_ref[unit.name, t]+reserve_up[unit.name, t]-(power_ref[unit.name, t-1]-reserve_down[unit.name, t-1])<=(-unit.DeltaRampUp)*start_up[unit.name, t]+(unit.MinPower+unit.DeltaRampUp)*is_on[unit.name, t]-(unit.MinPower)*is_on[unit.name, t-1])
    @constraint(model,  [unit in thermal_units, t in 1:T], power_ref[unit.name, t-1]+reserve_up[unit.name, t-1]-(power_ref[unit.name, t]-reserve_down[unit.name, t])<=(-unit.DeltaRampDown)*start_down[unit.name, t]+(unit.MinPower+unit.DeltaRampDown)*is_on[unit.name, t-1]-(unit.MinPower)*is_on[unit.name, t])

    @objective(model, Min, thermal_fixed_cost+thermal_cost)

    return model
end

function create_oracle_problem(instance; silent=true, S=1)
    """
    Create the master problem for the Benders' decomposition.
    """
    model = Model(instance.optimizer)
    if silent
        set_silent(model)
    end

    T= instance.TimeHorizon
    thermal_units_name=keys(instance.Thermalunits)
    thermal_units=values(instance.Thermalunits)

    Next=instance.Next
    Buses=1:size(Next)[1]
    @variable(model, power_ref[unit in thermal_units_name, t in 0:T]>=0)
    @variable(model, reserve_up[unit in thermal_units_name, t in 0:T]==0)
    @variable(model, reserve_down[unit in thermal_units_name, t in 0:T]==0)
    @variable(model, thermal_fuel_cost>=0)


    @variable(model, power_shedding[b in Buses, t in 0:T]>=0)
    @variable(model, power_curtailement[b in Buses, t in 0:T]>=0)

    @variable(model, θ[b in Buses, t in 1:T])
    @variable(model, flow[b in Buses, bp in Next[b], t in 1:T])
    
    Lines=values(instance.Lines)

    @constraint(model, [line in Lines, t in 1:T], flow[line.b1,line.b2,t]<=line.Fmax)
    @constraint(model, [line in Lines, t in 1:T], flow[line.b1,line.b2,t]>=-line.Fmax)
    @constraint(model, [line in Lines, t in 1:T], flow[line.b1,line.b2,t]==line.B12*(θ[line.b1,t]-θ[line.b2,t]))

    @variable(model, reg_up[unit in thermal_units_name, t in 0:T]>=0)
    @variable(model, reg_down[unit in thermal_units_name, t in 0:T]>=0)
    @variable(model, ξ[b in Buses, t in 1:T])


    @constraint(model, [unit in thermal_units_name, t in 0:T], reg_up[unit, t] <= reserve_up[unit, t])
    @constraint(model, [unit in thermal_units_name, t in 0:T], reg_down[unit, t] <= reserve_down[unit, t])

    @constraint(model, thermal_fuel_cost>= 50*sum(reg_up[unit.name, t] + reg_down[unit.name, t] for unit in thermal_units for t in 1:T) + sum(SHEDDING_COST*power_shedding[b,t]+CURTAILEMENT_COST*power_curtailement[b,t] for b in Buses for t in 1:T))

    @constraint(model,  [t in 1:T, b in Buses], sum(power_ref[unit.name, t]+reg_up[unit.name, t]-reg_down[unit.name, t] for unit in thermal_units if unit.Bus==b)+power_shedding[b,t]-power_curtailement[b,t]==instance.Demandbus[b][t]+sum(flow[b,bp,t] for bp in Next[b]) - ξ[b,t])

    @objective(model, Min, thermal_fuel_cost)

    return model
end

function fix_value_callback(cb_data, master_pb::JuMP.Model, oracle_pb::JuMP.Model, instance::Instance)
    """
    Fix the value of the primal variables in the master problem based on the oracle problem solution.
    """
    solution_power_ref=callback_value.(cb_data, master_pb[:power_ref])
    solution_reserve_up=callback_value.(cb_data, master_pb[:reserve_up])
    solution_reserve_down=callback_value.(cb_data, master_pb[:reserve_down])

    JuMP.fix.(oracle_pb[:power_ref], solution_power_ref; force=true)
    JuMP.fix.(oracle_pb[:reserve_up], solution_reserve_up; force=true)
    JuMP.fix.(oracle_pb[:reserve_down], solution_reserve_down; force=true)
end

function new_cut(cb_data, master_pb::JuMP.Model, oracle_pb::JuMP.Model, instance::Instance; force::Float64, S::Int64, batch=1, gap)
    """
    Fix the value of the primal variables in the master problem based on the oracle problem solution.
    """
    T = instance.TimeHorizon
    thermal_units = values(instance.Thermalunits)

    Next=instance.Next
    Buses=1:size(Next)[1]

    List_scenario=instance.Training_set[batch]

    intercept=0.0
    coef_power_ref = zeros(instance.N, T)
    coef_reserve_up = zeros(instance.N, T)
    coef_reserve_down = zeros(instance.N, T)
    for s in 1:S
        ξ_s = [instance.WGscenario[b][:,List_scenario[s]] for b in Buses]

        for b in Buses, t in 1:T
            JuMP.fix(oracle_pb[:ξ][b,t], force * ξ_s[b][t])
        end
        optimize!(oracle_pb)
        intercept += JuMP.objective_value(oracle_pb)/S

        coef_power_ref .+= convert(Matrix{Float64}, JuMP.dual.(JuMP.FixRef.(oracle_pb[:power_ref]))[:,1:T])./S
        coef_reserve_up .+= convert(Matrix{Float64}, JuMP.dual.(JuMP.FixRef.(oracle_pb[:reserve_up]))[:,1:T])./S
        coef_reserve_down .+= convert(Matrix{Float64}, JuMP.dual.(JuMP.FixRef.(oracle_pb[:reserve_down]))[:,1:T])./S
    end

    # println(value.(oracle_pb[:power_shedding]))

    thermal_cost_val=callback_value.(cb_data, master_pb[:thermal_cost])
    thermal_fixed_cost_val=callback_value.(cb_data, master_pb[:thermal_fixed_cost])
    solution_power_ref=callback_value.(cb_data, master_pb[:power_ref])
    solution_reserve_up=callback_value.(cb_data, master_pb[:reserve_up])
    solution_reserve_down=callback_value.(cb_data, master_pb[:reserve_down])

    # println("Intercept: ", intercept)
    # println("Thermal Cost Value: ", thermal_cost_val)
    # println("Thermal Fixed Cost Value: ", thermal_fixed_cost_val + thermal_cost_val, " ", intercept + thermal_fixed_cost_val)
    # println(solution_power_ref)
    # println(coef_power_ref)

    # println(1+"2")

    if intercept >= (100*gap / 100) * thermal_cost_val
        cstr=@build_constraint(intercept + sum(coef_power_ref[i, t] * (master_pb[:power_ref][i,t]-solution_power_ref[i,t]) + coef_reserve_up[i,t] * (master_pb[:reserve_up][i,t]-solution_reserve_up[i,t]) + coef_reserve_down[i,t] * (master_pb[:reserve_down][i,t]-solution_reserve_down[i,t]) for i in 1:instance.N, t in 1:T) <= master_pb[:thermal_cost])
        MOI.submit(master_pb, MOI.LazyConstraint(cb_data), cstr)

        cstr2=@build_constraint(master_pb[:thermal_cost] + master_pb[:thermal_fixed_cost] <= intercept + thermal_fixed_cost_val)
        MOI.submit(master_pb, MOI.LazyConstraint(cb_data), cstr2)
    end
    return intercept + thermal_fixed_cost_val, intercept
end