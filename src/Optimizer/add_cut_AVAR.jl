function add_cut_AVAR(cb_data, options::Options, master_pb::JuMP.Model, oracle_pb, instance::Instance, Time_subproblem, solution_x::Vector{Matrix{Float64}}; force::Float64, S::Int64, batch=1, gap, ρ)
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

    z_AVAR=master_pb[:z_AVAR]
    z_AVAR_val=callback_value(cb_data, z_AVAR)

    Time_iter=[0.0 for s in 1:S]

    Pmin, Pmax=get_limit_power_solution(instance, solution_x)

    @objective(oracle_pb, Max, sum(oracle_pb[:μₘᵢₙ][i,t]*Pmin[i,t] - oracle_pb[:μₘₐₓ][i,t]*Pmax[i,t] for i in 1:N2 for t in 1:T)+sum(oracle_pb[:network_cost][t] for t in 1:T))

    cut_parameters=[options.second_stage(instance, options, oracle_pb, prod_tot_first_stage_val, Pmin, Pmax; batch=batch, scenario=s, force=force) for s in 1:S]
    for s in 1:S
        Time_iter[s]=cut_parameters[s].computation_time
    end
    push!(Time_subproblem, Time_iter)

    if thermal_fixed_cost_val+dispatch_first_stage_val+z_AVAR_val+(sum([max(0,sum(cut_parameters[s].objective_value[t] for t in 1:T)-z_AVAR_val) for s in 1:S])/S)/(1-ρ)>= (1+0.01*gap/100)*(thermal_fixed_cost_val+dispatch_first_stage_val+z_AVAR_val+thermal_cost_val/(1-ρ))
        options._add_optimality_cuts(cb_data, master_pb, instance, cut_parameters; S=S)
    end

    return thermal_fixed_cost_val+dispatch_first_stage_val+z_AVAR_val+(sum([max(0,sum(cut_parameters[s].objective_value[t] for t in 1:T)-z_AVAR_val) for s in 1:S])/S)/(1-ρ), [cut_parameters[s].objective_value for s in 1:S]
end

function _add_optimality_cuts_AVAR(cb_data, master_pb, instance::Instance, cut_parameters::Vector{oracleResults}; S::Int64)
    T= instance.TimeHorizon
    thermal_fuel_cost=master_pb[:thermal_fuel_cost]
    prod_tot_first_stage=master_pb[:prod_tot_first_stage]
    is_on=master_pb[:is_on]
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
    end

end