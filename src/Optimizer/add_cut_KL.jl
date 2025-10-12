function add_cut_KL(cb_data, options::Options, master_pb::JuMP.Model, oracle_pb, instance::Instance, Time_subproblem, solution_x::Vector{Matrix{Float64}}; force::Float64, S::Int64, batch=1, gap, ρ)
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    N2 = N - N1

    thermal_cost=master_pb[:thermal_cost]
    thermal_fixed_cost=master_pb[:thermal_fixed_cost]
    dispatch_first_stage=master_pb[:dispatch_first_stage]
    prod_tot_first_stage=master_pb[:prod_tot_first_stage]

    thermal_cost_val=callback_value(cb_data, thermal_cost)
    thermal_fixed_cost_val=callback_value(cb_data, thermal_fixed_cost)
    dispatch_first_stage_val=callback_value(cb_data, dispatch_first_stage)
    prod_tot_first_stage_val=callback_value.(cb_data, prod_tot_first_stage)

    α=master_pb[:α]
    α_val=callback_value(cb_data, α)::Float64
    β=master_pb[:β]
    β_val=callback_value(cb_data, β)::Float64

    Time_iter=[0.0 for s in 1:S]

    Pmin, Pmax=get_limit_power_solution(instance, solution_x)

    @objective(oracle_pb, Max, sum(oracle_pb[:μₘᵢₙ][i,t]*Pmin[i,t] - oracle_pb[:μₘₐₓ][i,t]*Pmax[i,t] for i in 1:N2 for t in 1:T)+sum(oracle_pb[:network_cost][t] for t in 1:T))

    cut_parameters=[options.second_stage(instance, options, oracle_pb, prod_tot_first_stage_val, Pmin, Pmax; batch=batch, scenario=s, force=force) for s in 1:S]
    for s in 1:S
        Time_iter[s]=cut_parameters[s].computation_time
    end
    push!(Time_subproblem, Time_iter)

    second_stage_cost = zeros(S)
    for s in 1:S
        u = 0.0
        α_val = max(α_val, 1e-6)
        u = (sum(cut_parameters[s].objective_value[t] for t in 1:T)-β_val)/α_val
        u = min(u, 20.0)
        f = exp(u)-1
        second_stage_cost[s] = α_val * f
    end
    if thermal_fixed_cost_val+dispatch_first_stage_val+sum([second_stage_cost[s] for s in 1:S])/S + β_val + α_val*ρ>= (1+0.01*gap/100)*(thermal_fixed_cost_val+dispatch_first_stage_val+thermal_cost_val+ β_val + α_val*ρ)
        options._add_optimality_cuts(cb_data, master_pb, instance, cut_parameters; S=S)
    end

    return thermal_fixed_cost_val+dispatch_first_stage_val+sum([second_stage_cost[s] for s in 1:S])/S + β_val + α_val*ρ, [cut_parameters[s].objective_value for s in 1:S]
end

function _add_optimality_cuts_KL(cb_data, master_pb, instance::Instance, cut_parameters::Vector{oracleResults}; S::Int64)
    T= instance.TimeHorizon
    thermal_fuel_cost=master_pb[:thermal_fuel_cost]
    thermal_fuel_cost_KL=master_pb[:thermal_fuel_cost_KL]
    prod_tot_first_stage=master_pb[:prod_tot_first_stage]
    is_on=master_pb[:is_on]
    α=master_pb[:α]
    α_val=callback_value(cb_data, α)::Float64
    β=master_pb[:β]
    β_val=callback_value(cb_data, β)::Float64
    N=instance.N
    N1=instance.N1
    N2 = N - N1
    Next=instance.Next
    Buses=1:size(Next)[1] 
    thermal_units_2=values(instance.Thermalunits)[N1+1:N]

    for s in 1:S
        for t in 1:T
            cstr=@build_constraint(sum(cut_parameters[s].mumin[i, t]*thermal_units_2[i].MinPower*is_on[thermal_units_2[i].name,t] - cut_parameters[s].mumax[i, t]*thermal_units_2[i].MaxPower*is_on[thermal_units_2[i].name,t] for i in 1:N2) + cut_parameters[s].interceptE[t] - sum(cut_parameters[s].ν[b,t] * prod_tot_first_stage[b,t] for b in Buses) <= thermal_fuel_cost[s,t])
            MOI.submit(master_pb, MOI.LazyConstraint(cb_data), cstr)
        end
        u = 0.0
        if α_val>0
            u = (sum(cut_parameters[s].objective_value[t] for t in 1:T)-β_val)/α_val
        else
            u = (sum(cut_parameters[s].objective_value[t] for t in 1:T)-β_val)*1e6
        end
        u = min(u, 20.0)
        g = exp(u)
        f = exp(u)-1
        cstr2=@build_constraint(α*f + g*(sum(thermal_fuel_cost[s,t] for t in 1:T) - β - α*u) <= thermal_fuel_cost_KL[s])
        MOI.submit(master_pb, MOI.LazyConstraint(cb_data), cstr2)
    end

end