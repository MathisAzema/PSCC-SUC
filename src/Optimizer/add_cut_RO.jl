function add_cut_RO(cb_data, options::Options, master_pb::JuMP.Model, oracle_pb, instance::Instance, Time_subproblem, solution_x::Vector{Matrix{Float64}}; force::Float64, gap, Γ=0)
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

    cut_parameters=options.second_stage(instance, options, oracle_pb, prod_tot_first_stage_val, solution_x; Γ=Γ, force=force)

    push!(Time_subproblem, cut_parameters.computation_time)

    if thermal_fixed_cost_val+dispatch_first_stage_val+sum(cut_parameters.objective_value[t] for t in 1:T) >= (1+0.01*gap/100)*(thermal_fixed_cost_val+dispatch_first_stage_val+thermal_cost_val)
        _add_optimality_cuts_RO(cb_data, master_pb, instance, cut_parameters)
    end

    return thermal_fixed_cost_val+dispatch_first_stage_val+sum(cut_parameters.objective_value[t] for t in 1:T), cut_parameters.objective_value
end
    

function _add_optimality_cuts_RO(cb_data, master_pb, instance::Instance, cut_parameters::oracleResults)
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

    for t in 1:T
        cstr=@build_constraint(sum(cut_parameters.mumin[i, t]*thermal_units_2[i].MinPower*is_on[thermal_units_2[i].name,t] - cut_parameters.mumax[i, t]*thermal_units_2[i].MaxPower*is_on[thermal_units_2[i].name,t] for i in 1:N2) + cut_parameters.interceptE[t] - sum(cut_parameters.ν[b,t] * prod_tot_first_stage[b,t] for b in Buses) <= thermal_fuel_cost[t])
        MOI.submit(master_pb, MOI.LazyConstraint(cb_data), cstr)
    end

end