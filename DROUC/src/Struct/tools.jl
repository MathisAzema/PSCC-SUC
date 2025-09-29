function read_data_forecast(file)
    return CSV.read(file, DataFrame; delim=",", header=0)
end

function simulate_one_arma(; n=24, seed=123)
    phi1 = 0.670634
    theta1 = 0.117051
    c = 0.0
    w = 0.000291563

    rng = MersenneTwister(seed)
    y = zeros(Float64, n)
    ϵ = sqrt(w) * randn(rng, n)
    
    for t in 1:n
        ar_term = t > 1 ? phi1 * y[t-1] : 0.0
        ma_term = t > 1 ? theta1 * ϵ[t-1] : 0.0
        y[t] = c + ar_term + ma_term + ϵ[t]
    end
    return y
end

function simulate_scenarios(Nb_scenario)
    
    trajectories = [simulate_one_arma(; n=24, seed=123+i) for i in 1:Nb_scenario]

    return hcat(trajectories...)
end

function get_limit_power_solution(instance, solution)
    is_on=solution[1]
    start_up=solution[2]
    start_down=solution[3]
    
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    N2 = N - N1
    thermal_units_2=values(instance.Thermalunits)[N1+1:N]
    Pmin=Dict{Tuple{Int, Int}, Float64}()
    Pmax=Dict{Tuple{Int, Int}, Float64}()
    for i in 1:N2
        unit=thermal_units_2[i]
        if unit.InitialPower!=nothing
            Pmin[i,0]=unit.MinPower*is_on[unit.name,0+1]
            Pmax[i,0]=unit.MaxPower*is_on[unit.name,0+1]
        end
        for t in 1:T
            Pmin[i,t]=unit.MinPower*is_on[unit.name,t+1]
            Pmax[i,t]=unit.MaxPower*is_on[unit.name,t+1]
        end
    end
    return Pmin, Pmax
end

function get_variables_SP(instance, master_pb, S)
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]

    if N1>=1
        vars= vcat([master_pb[:is_on][i,t] for i in 1:N for t in 0:T], [master_pb[:start_up][i,t] for i in 1:N for t in 1:T], [master_pb[:start_down][i,t] for i in 1:N for t in 1:T],
            [master_pb[:thermal_fixed_cost], master_pb[:thermal_cost]], [master_pb[:thermal_fuel_cost][s,t] for s in 1:S for t in 1:T], 
            [master_pb[:power_first_stage][i,t] for i in 1:N1 for t in 1:T], [master_pb[:dispatch_first_stage]], [master_pb[:prod_tot_first_stage][b,t] for b in Buses for t in 1:T]
            )

        return vars
    else
        vars= vcat([master_pb[:is_on][i,t] for i in 1:N for t in 0:T], [master_pb[:start_up][i,t] for i in 1:N for t in 1:T], [master_pb[:start_down][i,t] for i in 1:N for t in 1:T], 
            [master_pb[:thermal_fixed_cost], master_pb[:thermal_cost]], [master_pb[:thermal_fuel_cost][s,t] for s in 1:S for t in 1:T], 
            [master_pb[:dispatch_first_stage]], [master_pb[:prod_tot_first_stage][b,t] for b in Buses for t in 1:T]
            )

        return vars
    end
end

function get_value_variables_SP(instance, master_pb, S, cb_data, solution_is_on, solution_start_up, solution_start_down, thermal_cost_scenario)
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]
    if N1>=1
        vals= vcat([solution_is_on[i,t+1] for i in 1:N for t in 0:T], [solution_start_up[i,t] for i in 1:N for t in 1:T], [solution_start_down[i,t] for i in 1:N for t in 1:T],
            [callback_value.(cb_data, master_pb[:thermal_fixed_cost]), sum(thermal_cost_scenario[s][t] for s in 1:S for t in 1:T)/S], [thermal_cost_scenario[s][t] for s in 1:S for t in 1:T],
            [callback_value.(cb_data, master_pb[:power_first_stage][i,t]) for i in 1:N1 for t in 1:T], callback_value.(cb_data, master_pb[:dispatch_first_stage]), [callback_value.(cb_data, master_pb[:prod_tot_first_stage][b,t]) for b in Buses for t in 1:T]
            )
        
        return vals
    else
        vals= vcat([solution_is_on[i,t+1] for i in 1:N for t in 0:T], [solution_start_up[i,t] for i in 1:N for t in 1:T], [solution_start_down[i,t] for i in 1:N for t in 1:T],
            [callback_value.(cb_data, master_pb[:thermal_fixed_cost]), sum(thermal_cost_scenario[s][t] for s in 1:S for t in 1:T)/S], [thermal_cost_scenario[s][t] for s in 1:S for t in 1:T],
            callback_value.(cb_data, master_pb[:dispatch_first_stage]), [callback_value.(cb_data, master_pb[:prod_tot_first_stage][b,t]) for b in Buses for t in 1:T]
            )
        return vals
    end
end

function get_variables_DRO_l2(instance, master_pb, S)
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]
    if N1 >= 1
        vars= vcat([master_pb[:is_on][i,t] for i in 1:N for t in 0:T], [master_pb[:start_up][i,t] for i in 1:N for t in 1:T], [master_pb[:start_down][i,t] for i in 1:N for t in 1:T],
            [master_pb[:thermal_fixed_cost], master_pb[:thermal_cost_DRO], master_pb[:thermal_cost]], [master_pb[:thermal_fuel_cost][s,t] for s in 1:S for t in 1:T], 
            [master_pb[:power_first_stage][i,t] for i in 1:N1 for t in 1:T], [master_pb[:dispatch_first_stage]], [master_pb[:prod_tot_first_stage][b,t] for b in Buses for t in 1:T],
            [master_pb[:λ], master_pb[:w]])
        
        return vars
    else
        vars= vcat([master_pb[:is_on][i,t] for i in 1:N for t in 0:T], [master_pb[:start_up][i,t] for i in 1:N for t in 1:T], [master_pb[:start_down][i,t] for i in 1:N for t in 1:T],
            [master_pb[:thermal_fixed_cost], master_pb[:thermal_cost_DRO], master_pb[:thermal_cost]], [master_pb[:thermal_fuel_cost][s,t] for s in 1:S for t in 1:T], 
            [master_pb[:λ], master_pb[:w]])
        
        return vars
    end
end

function get_value_variables_DRO_l2(instance, master_pb, S, cb_data, solution_is_on, solution_start_up, solution_start_down, thermal_cost_scenario)
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]
    if N1 >= 1
        vals= vcat([solution_is_on[i,t+1] for i in 1:N for t in 0:T], [solution_start_up[i,t] for i in 1:N for t in 1:T], [solution_start_down[i,t] for i in 1:N for t in 1:T],
            [callback_value.(cb_data, master_pb[:thermal_fixed_cost]), callback_value.(cb_data, master_pb[:thermal_cost_DRO]), sum(thermal_cost_scenario[s][t] for s in 1:S for t in 1:T)/S], [thermal_cost_scenario[s][t] for s in 1:S for t in 1:T],
            [callback_value.(cb_data, master_pb[:power_first_stage][i,t]) for i in 1:N1 for t in 1:T], callback_value.(cb_data, master_pb[:dispatch_first_stage]), [callback_value.(cb_data, master_pb[:prod_tot_first_stage][b,t]) for b in Buses for t in 1:T],
            [callback_value(cb_data, master_pb[:λ]), callback_value(cb_data, master_pb[:w])])
            return vals
    else
        vals= vcat([solution_is_on[i,t+1] for i in 1:N for t in 0:T], [solution_start_up[i,t] for i in 1:N for t in 1:T], [solution_start_down[i,t] for i in 1:N for t in 1:T],
            [callback_value.(cb_data, master_pb[:thermal_fixed_cost]), callback_value.(cb_data, master_pb[:thermal_cost_DRO]), sum(thermal_cost_scenario[s][t] for s in 1:S for t in 1:T)/S], [thermal_cost_scenario[s][t] for s in 1:S for t in 1:T],
            [callback_value(cb_data, master_pb[:λ]), callback_value(cb_data, master_pb[:w])])
        return vals
    end
end

function get_variables_AVAR(instance, master_pb, S)
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]
    if N1>=1
        vars= vcat([master_pb[:is_on][i,t] for i in 1:N for t in 0:T], [master_pb[:start_up][i,t] for i in 1:N for t in 1:T], [master_pb[:start_down][i,t] for i in 1:N for t in 1:T],
            [master_pb[:thermal_fixed_cost], master_pb[:thermal_cost]], [master_pb[:thermal_fuel_cost][s,t] for s in 1:S for t in 1:T], [master_pb[:thermal_fuel_cost_pos][s] for s in 1:S],
            [master_pb[:power_first_stage][i,t] for i in 1:N1 for t in 1:T], [master_pb[:dispatch_first_stage]], [master_pb[:prod_tot_first_stage][b,t] for b in Buses for t in 1:T], [master_pb[:z_AVAR]]
            )

        return vars
    else
        vars= vcat([master_pb[:is_on][i,t] for i in 1:N for t in 0:T], [master_pb[:start_up][i,t] for i in 1:N for t in 1:T], [master_pb[:start_down][i,t] for i in 1:N for t in 1:T], 
            [master_pb[:thermal_fixed_cost], master_pb[:thermal_cost]], [master_pb[:thermal_fuel_cost][s,t] for s in 1:S for t in 1:T], [master_pb[:thermal_fuel_cost_pos][s] for s in 1:S],
            [master_pb[:dispatch_first_stage]], [master_pb[:prod_tot_first_stage][b,t] for b in Buses for t in 1:T], [master_pb[:z_AVAR]]
            )

        return vars
    end
end

function get_value_variables_AVAR(instance, master_pb, S, cb_data, solution_is_on, solution_start_up, solution_start_down, thermal_cost_scenario)
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]
    z_AVAR_val=callback_value(cb_data, master_pb[:z_AVAR])
    if N1>=1
        vals= vcat([solution_is_on[i,t+1] for i in 1:N for t in 0:T], [solution_start_up[i,t] for i in 1:N for t in 1:T], [solution_start_down[i,t] for i in 1:N for t in 1:T],
            [callback_value.(cb_data, master_pb[:thermal_fixed_cost]), sum(max(0, sum(thermal_cost_scenario[s][t] for t in 1:T)-z_AVAR_val) for s in 1:S)/S], [thermal_cost_scenario[s][t] for s in 1:S for t in 1:T], [max(0, sum(thermal_cost_scenario[s][t] for t in 1:T)-z_AVAR_val) for s in 1:S],
            [callback_value.(cb_data, master_pb[:power_first_stage][i,t]) for i in 1:N1 for t in 1:T], callback_value.(cb_data, master_pb[:dispatch_first_stage]), [callback_value.(cb_data, master_pb[:prod_tot_first_stage][b,t]) for b in Buses for t in 1:T], z_AVAR_val
            )
        return vals

    else

        vals= vcat([solution_is_on[i,t+1] for i in 1:N for t in 0:T], [solution_start_up[i,t] for i in 1:N for t in 1:T], [solution_start_down[i,t] for i in 1:N for t in 1:T],
            [callback_value.(cb_data, master_pb[:thermal_fixed_cost]), sum(max(0, sum(thermal_cost_scenario[s][t] for t in 1:T)-z_AVAR_val) for s in 1:S)/S], [thermal_cost_scenario[s][t] for s in 1:S for t in 1:T], [max(0, sum(thermal_cost_scenario[s][t] for t in 1:T)-z_AVAR_val) for s in 1:S],
            callback_value.(cb_data, master_pb[:dispatch_first_stage]), [callback_value.(cb_data, master_pb[:prod_tot_first_stage][b,t]) for b in Buses for t in 1:T], z_AVAR_val
            )
        
        return vals 
    end
end

function get_variables_DRO_l1(instance, master_pb, S)
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]
    if N1 >= 1
        vars= vcat([master_pb[:is_on][i,t] for i in 1:N for t in 0:T], [master_pb[:start_up][i,t] for i in 1:N for t in 1:T], [master_pb[:start_down][i,t] for i in 1:N for t in 1:T],
            [master_pb[:thermal_fixed_cost], master_pb[:thermal_cost_DRO], master_pb[:thermal_cost]], [master_pb[:thermal_fuel_cost][s,t] for s in 1:S for t in 1:T], 
            [master_pb[:power_first_stage][i,t] for i in 1:N1 for t in 1:T], [master_pb[:dispatch_first_stage]], [master_pb[:prod_tot_first_stage][b,t] for b in Buses for t in 1:T],
            [master_pb[:λ]])
        
        return vars
    else
        vars= vcat([master_pb[:is_on][i,t] for i in 1:N for t in 0:T], [master_pb[:start_up][i,t] for i in 1:N for t in 1:T], [master_pb[:start_down][i,t] for i in 1:N for t in 1:T],
            [master_pb[:thermal_fixed_cost], master_pb[:thermal_cost_DRO], master_pb[:thermal_cost]], [master_pb[:thermal_fuel_cost][s,t] for s in 1:S for t in 1:T], 
            [master_pb[:λ]])
        
        return vars
    end
end

function get_value_variables_DRO_l1(instance, master_pb, S, cb_data, solution_is_on, solution_start_up, solution_start_down, thermal_cost_scenario)
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]
    if N1 >= 1
        vals= vcat([solution_is_on[i,t+1] for i in 1:N for t in 0:T], [solution_start_up[i,t] for i in 1:N for t in 1:T], [solution_start_down[i,t] for i in 1:N for t in 1:T],
            [callback_value.(cb_data, master_pb[:thermal_fixed_cost]), callback_value.(cb_data, master_pb[:thermal_cost_DRO]), sum(thermal_cost_scenario[s][t] for s in 1:S for t in 1:T)/S], [thermal_cost_scenario[s][t] for s in 1:S for t in 1:T],
            [callback_value.(cb_data, master_pb[:power_first_stage][i,t]) for i in 1:N1 for t in 1:T], callback_value.(cb_data, master_pb[:dispatch_first_stage]), [callback_value.(cb_data, master_pb[:prod_tot_first_stage][b,t]) for b in Buses for t in 1:T],
            [callback_value(cb_data, master_pb[:λ])])
            return vals
    else
        vals= vcat([solution_is_on[i,t+1] for i in 1:N for t in 0:T], [solution_start_up[i,t] for i in 1:N for t in 1:T], [solution_start_down[i,t] for i in 1:N for t in 1:T],
            [callback_value.(cb_data, master_pb[:thermal_fixed_cost]), callback_value.(cb_data, master_pb[:thermal_cost_DRO]), sum(thermal_cost_scenario[s][t] for s in 1:S for t in 1:T)/S], [thermal_cost_scenario[s][t] for s in 1:S for t in 1:T],
            [callback_value(cb_data, master_pb[:λ])])
        return vals
    end
end

function get_variables_KL(instance, master_pb, S)
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]
    if N1>=1
        vars= vcat([master_pb[:is_on][i,t] for i in 1:N for t in 0:T], [master_pb[:start_up][i,t] for i in 1:N for t in 1:T], [master_pb[:start_down][i,t] for i in 1:N for t in 1:T],
            [master_pb[:thermal_fixed_cost], master_pb[:thermal_cost]], [master_pb[:thermal_fuel_cost][s,t] for s in 1:S for t in 1:T], [master_pb[:thermal_fuel_cost_KL][s] for s in 1:S],
            [master_pb[:power_first_stage][i,t] for i in 1:N1 for t in 1:T], [master_pb[:dispatch_first_stage]], [master_pb[:prod_tot_first_stage][b,t] for b in Buses for t in 1:T], [master_pb[:α]], [master_pb[:β]]
            )

        return vars
    else
        vars= vcat([master_pb[:is_on][i,t] for i in 1:N for t in 0:T], [master_pb[:start_up][i,t] for i in 1:N for t in 1:T], [master_pb[:start_down][i,t] for i in 1:N for t in 1:T],
            [master_pb[:thermal_fixed_cost], master_pb[:thermal_cost]], [master_pb[:thermal_fuel_cost][s,t] for s in 1:S for t in 1:T], [master_pb[:thermal_fuel_cost_KL][s] for s in 1:S],
            [master_pb[:dispatch_first_stage]], [master_pb[:prod_tot_first_stage][b,t] for b in Buses for t in 1:T], [master_pb[:α]], [master_pb[:β]]
            )

        return vars
    end
end

function get_value_variables_KL(instance, master_pb, S, cb_data, solution_is_on, solution_start_up, solution_start_down, thermal_cost_scenario)
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]

    second_stage_cost = zeros(S)
    for s in 1:S
        u = 0.0
        if α_val>0
            u = (sum(thermal_cost_scenario[s][t] for t in 1:T)-β_val)/α_val
        else
            u = (sum(thermal_cost_scenario[s][t] for t in 1:T)-β_val)*1e6
        end
        u = min(u, 20.0)
        f = exp(u)-1
        second_stage_cost[s] = α_val * f
    end

    if N1>=1
        vals= vcat([solution_is_on[i,t+1] for i in 1:N for t in 0:T], [solution_start_up[i,t] for i in 1:N for t in 1:T], [solution_start_down[i,t] for i in 1:N for t in 1:T],
            [callback_value.(cb_data, master_pb[:thermal_fixed_cost]), sum(second_stage_cost[s] for s in 1:S)/S], [thermal_cost_scenario[s][t] for s in 1:S for t in 1:T], [second_stage_cost[s] for s in 1:S],
            [callback_value.(cb_data, master_pb[:power_first_stage][i,t]) for i in 1:N1 for t in 1:T], callback_value.(cb_data, master_pb[:dispatch_first_stage]), [callback_value.(cb_data, master_pb[:prod_tot_first_stage][b,t]) for b in Buses for t in 1:T], [callback_value(cb_data, master_pb[:α]), callback_value(cb_data, master_pb[:β])]
            )
        
        return vals
    else
        vals= vcat([solution_is_on[i,t+1] for i in 1:N for t in 0:T], [solution_start_up[i,t] for i in 1:N for t in 1:T], [solution_start_down[i,t] for i in 1:N for t in 1:T],
            [callback_value.(cb_data, master_pb[:thermal_fixed_cost]), sum(second_stage_cost[s] for s in 1:S)/S], [thermal_cost_scenario[s][t] for s in 1:S for t in 1:T], [second_stage_cost[s] for s in 1:S],
            callback_value.(cb_data, master_pb[:dispatch_first_stage]), [callback_value.(cb_data, master_pb[:prod_tot_first_stage][b,t]) for b in Buses for t in 1:T], [callback_value(cb_data, master_pb[:α]), callback_value(cb_data, master_pb[:β])]
            )
        return vals
    end
end

function get_variables_RO(instance, master_pb)
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]

    if N1>=1
        vars= vcat([master_pb[:is_on][i,t] for i in 1:N for t in 0:T], [master_pb[:start_up][i,t] for i in 1:N for t in 1:T], [master_pb[:start_down][i,t] for i in 1:N for t in 1:T],
            [master_pb[:thermal_fixed_cost], master_pb[:thermal_cost]], [master_pb[:thermal_fuel_cost][t] for t in 1:T], 
            [master_pb[:power_first_stage][i,t] for i in 1:N1 for t in 1:T], [master_pb[:dispatch_first_stage]], [master_pb[:prod_tot_first_stage][b,t] for b in Buses for t in 1:T]
            )

        return vars
    else
        vars= vcat([master_pb[:is_on][i,t] for i in 1:N for t in 0:T], [master_pb[:start_up][i,t] for i in 1:N for t in 1:T], [master_pb[:start_down][i,t] for i in 1:N for t in 1:T], 
            [master_pb[:thermal_fixed_cost], master_pb[:thermal_cost]], [master_pb[:thermal_fuel_cost][t] for t in 1:T], 
            [master_pb[:dispatch_first_stage]], [master_pb[:prod_tot_first_stage][b,t] for b in Buses for t in 1:T]
            )

        return vars
    end
end

function get_value_variables_RO(instance, master_pb, cb_data, solution_is_on, solution_start_up, solution_start_down, thermal_cost_scenario)
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    Next=instance.Next
    Buses=1:size(Next)[1]
    if N1>=1
        vals= vcat([solution_is_on[i,t+1] for i in 1:N for t in 0:T], [solution_start_up[i,t] for i in 1:N for t in 1:T], [solution_start_down[i,t] for i in 1:N for t in 1:T],
            [callback_value.(cb_data, master_pb[:thermal_fixed_cost]), sum(thermal_cost_scenario[t] for t in 1:T)], [thermal_cost_scenario[t] for t in 1:T],
            [callback_value.(cb_data, master_pb[:power_first_stage][i,t]) for i in 1:N1 for t in 1:T], callback_value.(cb_data, master_pb[:dispatch_first_stage]), [callback_value.(cb_data, master_pb[:prod_tot_first_stage][b,t]) for b in Buses for t in 1:T]
            )
        
        return vals
    else
        vals= vcat([solution_is_on[i,t+1] for i in 1:N for t in 0:T], [solution_start_up[i,t] for i in 1:N for t in 1:T], [solution_start_down[i,t] for i in 1:N for t in 1:T],
            [callback_value.(cb_data, master_pb[:thermal_fixed_cost]), sum(thermal_cost_scenario[t] for t in 1:T)], [thermal_cost_scenario[t] for t in 1:T],
            callback_value.(cb_data, master_pb[:dispatch_first_stage]), [callback_value.(cb_data, master_pb[:prod_tot_first_stage][b,t]) for b in Buses for t in 1:T]
            )
        return vals
    end
end

function compute_radius_l2(S, ρ)
    f = S^0.25
    radius = ρ/f
    return radius
end

function compute_radius_l1(S, ρ)
    f = S^0.5
    radius = ρ/f
    return radius
end

function compute_radius_KL(S, ρ)
    df = S - 1
    dist = Chisq(df)
    return quantile(dist, ρ)/(2*S)
end

function compute_radius_AVAR(S, ρ)
    return ρ
end

function compute_radius_SP(S, ρ)
    return ρ
end