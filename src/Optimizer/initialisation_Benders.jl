struct Options
    name::String
    master_problem::Function
    oracle_problem::Function
    second_stage::Function
    add_cut::Function
    _add_optimality_cuts::Function
    get_variables::Function
    get_value_variables::Function
    compute_uncertainty::Function
    results_second_stage::Function
    compute_radius::Function
end

struct first_stage
    power_ref::Matrix{Float64}
    is_on::Matrix{Int}
    start_up::Matrix{Int}
    start_down::Matrix{Int}
end

function initialize_master_problem(instance; silent=true, S::Int64)
    model = Model(instance.optimizer)
    if silent
        set_silent(model)
    end

    set_optimizer_attribute(model, "Threads", 1)

    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]

    @variable(model, is_on[i in 1:N, t in 0:T], Bin)
    @variable(model, start_up[i in 1:N, t in 1:T], Bin)
    @variable(model, start_down[i in 1:N, t in 1:T], Bin)

    @variable(model, thermal_fuel_cost[s in 1:S, t in 1:T]>=0)
    @variable(model, thermal_fixed_cost>=0)
    @variable(model, thermal_cost>=0)

    thermal_unit_commit_constraints(model, instance)

    @variable(model, power_first_stage[i in 1:N1, t in 0:T])
    @variable(model, dispatch_first_stage>=0)
    @variable(model, prod_tot_first_stage[b in Buses, t in 1:T]>=0)
    thermal_units_1=values(instance.Thermalunits)[1:N1]

    @constraint(model,  [unit in thermal_units_1, t in 0:T], power_first_stage[unit.name, t]>=unit.MinPower*is_on[unit.name, t])
    @constraint(model,  [unit in thermal_units_1, t in 0:T], power_first_stage[unit.name, t]<=unit.MaxPower*is_on[unit.name, t])
    @constraint(model,  [unit in thermal_units_1], power_first_stage[unit.name, 0]==unit.InitialPower)

    @constraint(model,  [unit in thermal_units_1, t in 1:T], power_first_stage[unit.name, t]-power_first_stage[unit.name, t-1]<=(-unit.DeltaRampUp)*start_up[unit.name, t]+(unit.MinPower+unit.DeltaRampUp)*is_on[unit.name, t]-(unit.MinPower)*is_on[unit.name, t-1])
    @constraint(model,  [unit in thermal_units_1, t in 1:T], power_first_stage[unit.name, t-1]-power_first_stage[unit.name, t]<=(-unit.DeltaRampDown)*start_down[unit.name, t]+(unit.MinPower+unit.DeltaRampDown)*is_on[unit.name, t-1]-(unit.MinPower)*is_on[unit.name, t])

    @constraint(model, dispatch_first_stage >= sum(unit.LinearTerm*power_first_stage[unit.name, t] for unit in thermal_units_1 for t in 1:T))
    @constraint(model, [b in Buses, t in 1:T], prod_tot_first_stage[b, t] == sum(power_first_stage[unit.name, t] for unit in thermal_units_1 if unit.Bus == b))

    @variable(model, obj>=0)
    @objective(model, Min, obj)

    return model
end

function initialize_oracle_problem(instance)
    model = Model(instance.optimizer)
    set_silent(model)
    T= instance.TimeHorizon

    N=instance.N
    N1=instance.N1
    N2 = N - N1
    thermal_units_2=values(instance.Thermalunits)[N1+1:N]

    Next=instance.Next
    Buses=1:size(Next)[1]

    Lines=instance.Lines
    Numlines=length(Lines)

    @variable(model, μₘᵢₙ[i in 1:N2, t in 1:T]>=0)
    @variable(model, μₘₐₓ[i in 1:N2, t in 1:T]>=0)
    @variable(model, ν[b in Buses, t in 1:T])

    @constraint(model,  [i in 1:N2, t in 1:T], μₘᵢₙ[i, t]-μₘₐₓ[i, t]+ν[thermal_units_2[i].Bus, t]==thermal_units_2[i].LinearTerm)

    @constraint(model,  power_shedding[b in Buses, t in 1:T], ν[b, t]<=SHEDDING_COST)
    @constraint(model,  [b in Buses, t in 1:T], ν[b, t]>=-CURTAILEMENT_COST)

    @variable(model, μ1[l in 1:Numlines, t in 1:T])
    @variable(model, μ2[l in 1:Numlines, t in 1:T]>=0)
    @variable(model, μ3[l in 1:Numlines, t in 1:T]>=0)

    Lines=instance.Lines
    @variable(model, network_cost[t in 1:T])
    @constraint(model, [t in 1:T], network_cost[t]==-sum(line.Fmax*(μ2[line.id,t]+μ3[line.id,t]) for line in Lines))
    @constraint(model, [line in Lines, t in 1:T], ν[line.b2, t] - ν[line.b1, t] + μ1[line.id,t] - μ2[line.id,t] + μ3[line.id,t]==0)
    @constraint(model, [b in Buses, t in 1:T], sum(line.B12*μ1[line.id,t] for line in Lines if line.b1==b) - sum(line.B12*μ1[line.id,t] for line in Lines if line.b2==b) == 0)

    set_optimizer_attribute(model, "Threads", 1)
    return model
end

function master_SP_problem(instance; silent=true, ρ=0, S::Int64)
    """
    Initial master problem in the risk neutral extended formulation
    """

    model = initialize_master_problem(instance; silent=silent, S=S)

    T= instance.TimeHorizon
    @constraint(model,  model[:thermal_cost]>=sum(model[:thermal_fuel_cost][s,t] for t in 1:T for s in 1:S)/S)

    @constraint(model, model[:obj] >= model[:thermal_fixed_cost]+model[:thermal_cost] + model[:dispatch_first_stage])
    
    return model
end

mutable struct oracleSP 
    model::Model
    μₘᵢₙ::Matrix{VariableRef}
    μₘₐₓ::Matrix{VariableRef}
    ν::Matrix{VariableRef}
    λ::Vector{VariableRef}
    network_cost::VariableRef
end

function oracle_SP_problem(instance)
    """
    Initial the subproblem in its dual form
    """

    model = initialize_oracle_problem(instance)
    set_optimizer_attribute(model, "TimeLimit", 15)
    set_optimizer_attribute(model, "DualReductions", 0)
    set_optimizer_attribute(model, "Presolve", 0)

    return model
end

function master_DRO_l2_problem(instance; silent=true, ρ=0, S::Int64)

    model = initialize_master_problem(instance; silent=silent, S=S)

    T= instance.TimeHorizon
    @constraint(model,  model[:thermal_cost]>=sum(model[:thermal_fuel_cost][s,t] for t in 1:T for s in 1:S)/S)

    @variable(model, thermal_cost_DRO>=0)

    @variable(model, λ>=0)
    @variable(model, w>=0)
    @constraint(model, 1<=w+λ)

    R=500
    λᵣ=[1e-3*10*i/R for i in 1:R]
    @constraint(model,  [r in 1:R], 1/λᵣ[r]-1/λᵣ[r]^2*(λ-λᵣ[r]) <=w)

    @constraint(model,  thermal_cost_DRO>=sum(ρ^2*w/1e-2))

    @constraint(model, model[:obj] >= model[:thermal_fixed_cost]+model[:thermal_cost_DRO]+model[:thermal_cost] + model[:dispatch_first_stage])
    
    return model
end

function master_DRO_l1_problem(instance; silent=true, ρ=0, S::Int64)

    model = initialize_master_problem(instance; silent=silent, S=S)

    T= instance.TimeHorizon
    @constraint(model,  model[:thermal_cost]>=sum(model[:thermal_fuel_cost][s,t] for t in 1:T for s in 1:S)/S)

    @variable(model, thermal_cost_DRO>=0)

    @variable(model, λ>=0)

    @constraint(model,  thermal_cost_DRO>=sum(ρ*λ*1e-2))

    @constraint(model, model[:obj] >= model[:thermal_fixed_cost]+model[:thermal_cost_DRO]+model[:thermal_cost] + model[:dispatch_first_stage])
    
    return model
end

function oracle_DRO_l2_problem(instance)

    model = initialize_oracle_problem(instance)
    set_optimizer_attribute(model, "TimeLimit", 15)
    
    return model
end

function oracle_DRO_l1_problem(instance)

    model = initialize_oracle_problem(instance)
    set_optimizer_attribute(model, "TimeLimit", 15)

    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)
    T= instance.TimeHorizon

    @variable(model, zp[w in 1:NumWindfarms, t in 1:T], Bin)
    @variable(model, zm[w in 1:NumWindfarms, t in 1:T], Bin)
    @constraint(model, [w in 1:NumWindfarms, t in 1:T], zp[w,t]+zm[w,t]<=1)

    return model
end

function master_AVAR_problem(instance; silent=true, ρ=0,S::Int64)

    """
    Initial master problem in the AVAR extended formulation
    """

    model = initialize_master_problem(instance; silent=silent, S=S)

    T= instance.TimeHorizon

    @variable(model, thermal_fuel_cost_pos[s in 1:S]>=0)
    @variable(model, z_AVAR>=0)

    thermal_unit_commit_constraints(model, instance)

    @constraint(model,  model[:thermal_cost]>=sum(thermal_fuel_cost_pos[s] for s in 1:S)/S)
    @constraint(model,  [s in 1:S], thermal_fuel_cost_pos[s]>=sum(model[:thermal_fuel_cost][s,t] for t in 1:T)-z_AVAR)

    @constraint(model, model[:obj] >= model[:thermal_fixed_cost]+ model[:dispatch_first_stage]+z_AVAR+model[:thermal_cost]/(1-ρ))
    
    return model
end

function master_KL_problem(instance; silent=true, ρ=0,S::Int64)

    model = initialize_master_problem(instance; silent=silent, S=S)

    T= instance.TimeHorizon
    @variable(model, thermal_fuel_cost_KL[s in 1:S]>=0)
    @constraint(model,  model[:thermal_cost]>=sum(thermal_fuel_cost_KL[s] for s in 1:S)/S)

    T= instance.TimeHorizon

    @variable(model, α>=0)
    @variable(model, β)

    thermal_unit_commit_constraints(model, instance)
    @constraint(model,  [s in 1:S], thermal_fuel_cost_KL[s]>=sum(model[:thermal_fuel_cost][s,t] for t in 1:T)- β)

    @constraint(model, model[:obj] >= model[:thermal_fixed_cost]+ model[:dispatch_first_stage] +model[:thermal_cost]+β + α*ρ)

    return model
end

function master_RO_problem_CCG(instance; silent=true)
    """
    Initial master problem in the risk neutral extended formulation
    """

    model = Model(instance.optimizer)
    if silent
        set_silent(model)
    end

    set_optimizer_attribute(model, "Threads", 1)

    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]

    @variable(model, is_on[i in 1:N, t in 0:T], Bin)
    @variable(model, start_up[i in 1:N, t in 1:T], Bin)
    @variable(model, start_down[i in 1:N, t in 1:T], Bin)

    @variable(model, thermal_fixed_cost>=0)
    @variable(model, thermal_cost>=0)

    thermal_unit_commit_constraints(model, instance)

    @variable(model, power_first_stage[i in 1:N1, t in 0:T])
    @variable(model, dispatch_first_stage>=0)
    @variable(model, prod_tot_first_stage[b in Buses, t in 1:T]>=0)
    thermal_units_1=values(instance.Thermalunits)[1:N1]

    @variable(model, power[i in 1:N, t in 0:T]>=0)
    @variable(model, power_shedding[b in Buses, t in 0:T]>=0)
    @variable(model, power_curtailement[b in Buses, t in 0:T]>=0)

    @variable(model, θ[b in Buses, t in 1:T])
    @variable(model, flow[b in Buses, bp in Next[b], t in 1:T])

    @constraint(model,  [unit in thermal_units_1, t in 0:T], power_first_stage[unit.name, t]>=unit.MinPower*is_on[unit.name, t])
    @constraint(model,  [unit in thermal_units_1, t in 0:T], power_first_stage[unit.name, t]<=unit.MaxPower*is_on[unit.name, t])
    @constraint(model,  [unit in thermal_units_1], power_first_stage[unit.name, 0]==unit.InitialPower)

    @constraint(model,  [unit in thermal_units_1, t in 1:T], power_first_stage[unit.name, t]-power_first_stage[unit.name, t-1]<=(-unit.DeltaRampUp)*start_up[unit.name, t]+(unit.MinPower+unit.DeltaRampUp)*is_on[unit.name, t]-(unit.MinPower)*is_on[unit.name, t-1])
    @constraint(model,  [unit in thermal_units_1, t in 1:T], power_first_stage[unit.name, t-1]-power_first_stage[unit.name, t]<=(-unit.DeltaRampDown)*start_down[unit.name, t]+(unit.MinPower+unit.DeltaRampDown)*is_on[unit.name, t-1]-(unit.MinPower)*is_on[unit.name, t])

    @constraint(model, dispatch_first_stage >= sum(unit.LinearTerm*power_first_stage[unit.name, t] for unit in thermal_units_1 for t in 1:T))
    @constraint(model, [b in Buses, t in 1:T], prod_tot_first_stage[b, t] == sum(power_first_stage[unit.name, t] for unit in thermal_units_1 if unit.Bus == b))

    @objective(model, Min, thermal_fixed_cost+thermal_cost + dispatch_first_stage)
    
    return model
end

function oracle_RO_problem(instance; Γ::Int64=0)
    
    T= instance.TimeHorizon

    BusWind=instance.BusWind

    model = initialize_oracle_problem(instance)
    ν=model[:ν]

    @variable(model, γ[b in BusWind, t in 1:T], Bin)
    
    @constraint(model, [t in 1:T],  sum([γ[b, t] for b in BusWind])<=Γ)

    @variable(model, ϵ[b in BusWind, t in 1:T], Bin)
    @variable(model, ν⁺[b in BusWind, t in 1:T]>=0)
    @variable(model, ν⁻[b in BusWind, t in 1:T]>=0)

    ub=SHEDDING_COST
    lb=CURTAILEMENT_COST
    @constraint(model,  [b in BusWind, t in 1:T], ν⁺[b, t]<=ub*ϵ[b, t])
    @constraint(model,  [b in BusWind, t in 1:T], ν⁻[b,t]<=lb*(1-ϵ[b, t]))
    @constraint(model,  [b in BusWind, t in 1:T], ν[b, t]==ν⁺[b,t]-ν⁻[b,t])

    @variable(model, ζ[b in BusWind, t in 1:T]>=0)
    @constraint(model,  [b in BusWind, t in 1:T], ζ[b,t]<=max(lb,ub)*γ[b, t])
    @constraint(model,  [b in BusWind, t in 1:T], ζ[b, t]<=ν⁺[b,t]+ν⁻[b,t])

    set_optimizer_attribute(model, "Threads", 1)
    set_optimizer_attribute(model, "TimeLimit", 15)
    
    return model
end

function master_RO_problem_benders(instance; silent=true)
    """
    Initial master problem in the robust case
    """

    model = Model(instance.optimizer)
    if silent
        set_silent(model)
    end

    set_optimizer_attribute(model, "Threads", 1)

    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]

    @variable(model, is_on[i in 1:N, t in 0:T], Bin)
    @variable(model, start_up[i in 1:N, t in 1:T], Bin)
    @variable(model, start_down[i in 1:N, t in 1:T], Bin)

    @variable(model, thermal_fixed_cost>=0)
    @variable(model, thermal_cost>=0)
    @variable(model, thermal_fuel_cost[t in 1:T]>=0)

    @constraint(model, thermal_cost >= sum(thermal_fuel_cost[t] for t in 1:T))

    thermal_unit_commit_constraints(model, instance)

    @variable(model, power_first_stage[i in 1:N1, t in 0:T])
    @variable(model, dispatch_first_stage>=0)
    @variable(model, prod_tot_first_stage[b in Buses, t in 1:T]>=0)
    thermal_units_1=values(instance.Thermalunits)[1:N1]

    @constraint(model,  [unit in thermal_units_1, t in 0:T], power_first_stage[unit.name, t]>=unit.MinPower*is_on[unit.name, t])
    @constraint(model,  [unit in thermal_units_1, t in 0:T], power_first_stage[unit.name, t]<=unit.MaxPower*is_on[unit.name, t])
    @constraint(model,  [unit in thermal_units_1], power_first_stage[unit.name, 0]==unit.InitialPower)

    @constraint(model,  [unit in thermal_units_1, t in 1:T], power_first_stage[unit.name, t]-power_first_stage[unit.name, t-1]<=(-unit.DeltaRampUp)*start_up[unit.name, t]+(unit.MinPower+unit.DeltaRampUp)*is_on[unit.name, t]-(unit.MinPower)*is_on[unit.name, t-1])
    @constraint(model,  [unit in thermal_units_1, t in 1:T], power_first_stage[unit.name, t-1]-power_first_stage[unit.name, t]<=(-unit.DeltaRampDown)*start_down[unit.name, t]+(unit.MinPower+unit.DeltaRampDown)*is_on[unit.name, t-1]-(unit.MinPower)*is_on[unit.name, t])

    @constraint(model, dispatch_first_stage >= sum(unit.LinearTerm*power_first_stage[unit.name, t] for unit in thermal_units_1 for t in 1:T))
    @constraint(model, [b in Buses, t in 1:T], prod_tot_first_stage[b, t] == sum(power_first_stage[unit.name, t] for unit in thermal_units_1 if unit.Bus == b))
    
    @variable(model, obj>=0)
    @constraint(model, obj >= thermal_fixed_cost+thermal_cost + dispatch_first_stage)
    @objective(model, Min, obj)
    
    return model
end

function oracle_RO_problem_GRB(instance; Γ::Int64=0)

    T= instance.TimeHorizon

    BusWind=instance.BusWind

    model = initialize_oracle_problem(instance)
    ν=model[:ν]

    @variable(model, γ[b in BusWind, t in 1:T], Bin)

    @constraint(model, [t in 1:T],  sum([γ[b, t] for b in BusWind])<=Γ)

    @variable(model, δp[b in BusWind, t in 1:T], Bin)
    @variable(model, δm[b in BusWind, t in 1:T], Bin)

    @constraint(model,  [b in BusWind, t in 1:T], δp[b, t]+δm[b, t]<=γ[b, t])

    set_optimizer_attribute(model, "Threads", 1)
    set_optimizer_attribute(model, "TimeLimit", 15)
    
    return model
end