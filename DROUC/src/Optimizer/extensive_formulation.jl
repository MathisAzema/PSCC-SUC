struct first_stage
    power_ref::Matrix{Float64}
    is_on::Matrix{Int}
    start_up::Matrix{Int}
    start_down::Matrix{Int}
end

function bin_extensive_neutral(instance; silent=true,  force::Float64=1.0, S::Int64=5, batch=1, gap=gap, timelimit=10)
    """
    Solve the two stage SUC with 3-bin extensive formulation
    """

    model = Model(instance.optimizer)
    if silent
        set_silent(model)
    end

    T= instance.TimeHorizon
    thermal_units_name=keys(instance.Thermalunits)
    thermal_units=values(instance.Thermalunits)

    N1=instance.N1
    N2 = instance.N - N1
    thermal_units_1=values(instance.Thermalunits)[1:N1]

    Next=instance.Next
    Buses=1:size(Next)[1]
    @variable(model, power_ref[i in 1:N1, t in 0:T]>=0)
    @variable(model, is_on[unit in thermal_units_name, t in 0:T], Bin)
    @variable(model, start_up[unit in thermal_units_name, t in 1:T], Bin)
    @variable(model, start_down[unit in thermal_units_name, t in 1:T], Bin)

    @variable(model, thermal_fuel_cost[s in 1:S]>=0)
    @variable(model, thermal_fixed_cost>=0)
    @variable(model, thermal_cost>=0)

    @variable(model, power_shedding[b in Buses, t in 0:T, s in 1:S]>=0)
    @variable(model, power_curtailement[b in Buses, t in 0:T, s in 1:S]>=0)

    @constraint(model,  thermal_cost>=sum(thermal_fuel_cost[s] for s in 1:S)/S)

    thermal_unit_commit_constraints(model, instance)

    @variable(model, power[unit in thermal_units_name, t in 0:T, s in 1:S]>=0)

    @constraint(model, [i in 1:N1, t in 0:T, s in 1:S], power[i, t, s] == power_ref[i, t] )

    @constraint(model,  [unit in thermal_units, s in 1:S; unit.InitialPower!=nothing], power[unit.name, 0, s]==unit.InitialPower)
    @constraint(model,  [unit in thermal_units, t in 0:T, s in 1:S], power[unit.name, t, s]>=unit.MinPower*is_on[unit.name, t])
    @constraint(model,  [unit in thermal_units, t in 0:T, s in 1:S], power[unit.name, t, s]<=unit.MaxPower*is_on[unit.name, t])

    @constraint(model,  [unit in thermal_units_1, t in 1:T, s in 1:S], power[unit.name, t, s]-power[unit.name, t-1, s]<=(-unit.DeltaRampUp)*start_up[unit.name, t]+(unit.MinPower+unit.DeltaRampUp)*is_on[unit.name, t]-(unit.MinPower)*is_on[unit.name, t-1])
    @constraint(model,  [unit in thermal_units_1, t in 1:T, s in 1:S], power[unit.name, t-1, s]-power[unit.name, t, s]<=(-unit.DeltaRampDown)*start_down[unit.name, t]+(unit.MinPower+unit.DeltaRampDown)*is_on[unit.name, t-1]-(unit.MinPower)*is_on[unit.name, t])

    @constraint(model, [s in 1:S], thermal_fuel_cost[s]>= sum(unit.LinearTerm*power[unit.name, t, s] for unit in thermal_units for t in 1:T) + sum(SHEDDING_COST*power_shedding[b,t,s]+CURTAILEMENT_COST*power_curtailement[b,t,s] for b in Buses for t in 1:T))

    @variable(model, θ[b in Buses, t in 1:T, s in 1:S])
    @variable(model, flow[b in Buses, bp in Next[b], t in 1:T, s in 1:S])
    flow_constraints_scenarios(model,instance,S)

    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)
    Wpower=instance.WGscenario
    List_scenario=instance.Training_set[batch]

    @constraint(model,  demand[t in 1:T, b in Buses, s in 1:S], sum(power[unit.name, t, s] for unit in thermal_units if unit.Bus==b)+power_shedding[b,t,s]-power_curtailement[b,t,s]==instance.Demandbus[b][t]+sum(flow[b,bp,t,s] for bp in Next[b])-force * sum(Wpower[w][t,List_scenario[s]] for w in 1:NumWindfarms if BusWind[w]==b))

    @objective(model, Min, thermal_fixed_cost+thermal_cost)
  
    set_optimizer_attribute(model, "TimeLimit", timelimit)
    set_optimizer_attribute(model, "Threads", 1)
    set_optimizer_attribute(model, "MIPGap", gap)

    start = time()
    optimize!(model)
    computation_time = time() - start

    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT

    if feasibleSolutionFound

        solution_is_on=JuMP.value.(is_on)
        solution_start_up=JuMP.value.(start_up)
        solution_start_down=JuMP.value.(start_down)
        solution_power_ref=JuMP.value.(power_ref)

        first_stage_solution = first_stage(solution_power_ref, solution_is_on, solution_start_up, solution_start_down)

        return instance.name, computation_time, objective_value(model), objective_bound(model), first_stage_solution, S, batch,  gap, force, value.(power)

    else
        return nothing
    end
end