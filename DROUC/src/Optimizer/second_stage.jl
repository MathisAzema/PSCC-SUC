mutable struct oracleResults
    objective_value::Vector{Float64}
    intercept3bin::Vector{Float64}
    interceptE::Vector{Float64}
    quad::Vector{Float64}
    mumax::Matrix{Float64}
    mumin::Matrix{Float64}
    ν::Matrix{Float64}
    status::Bool
    computation_time::Float64
end

function solve_second_stage(instance, uncertainty, oracle_pb, prod_tot_first_stage; force=1.0)
    T= instance.TimeHorizon
    Next=instance.Next
    Buses=1:size(Next)[1] 

    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)

    uncertainty_0 = instance.uncertainty

    ζ = [[uncertainty_0.dev[t]*uncertainty_0.forecast[w][t]*uncertainty[w][t] for t in 1:T] for w in 1:NumWindfarms]

    Netdemand = [instance.Demandbus[b][t] - force * sum(ζ[w][t] for w in 1:NumWindfarms if BusWind[w] == b;init=0) for b in Buses, t in 1:T]
    
    set_objective_coefficient.(oracle_pb,oracle_pb[:ν], Netdemand-prod_tot_first_stage)

    start = time()
    optimize!(oracle_pb)
    computation_time = time() - start

    return computation_time, Netdemand
end

function second_stage_SP(instance, options, oracle_pb, prod_tot_first_stage, Pmin, Pmax; λ=0.0, batch=1, scenario=1, force=1)

    uncertainty = options.compute_uncertainty(instance; batch=batch, scenario=scenario, force=force)

    computation_time, Netdemand = solve_second_stage(instance, uncertainty, oracle_pb, prod_tot_first_stage; force=force)

    return options.results_second_stage(instance, oracle_pb, prod_tot_first_stage, Pmin, Pmax, computation_time, Netdemand; λ=λ, batch=batch, scenario=scenario, force=force)
end

function second_stage_DCA(instance, options, oracle_pb, uncertainty, prod_tot_first_stage, Pmin, Pmax, λ; batch=1, scenario=1, force=1)

    computation_time, Lindemand = solve_second_stage(instance, uncertainty, oracle_pb, prod_tot_first_stage; force=force)

    return options.results_second_stage(instance, oracle_pb, prod_tot_first_stage, Pmin, Pmax, computation_time, Lindemand; λ=λ, batch=batch, scenario=scenario, force=force)
end

function second_stage_grb_l2(instance, options, oracle_pb, prod_tot_first_stage, Pmin, Pmax; λ=0.0, batch=1, scenario=1, force=1)

    T= instance.TimeHorizon
    Next=instance.Next
    Buses=1:size(Next)[1] 

    ν=oracle_pb[:ν]

    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)

    List_scenario=instance.Training_set[batch]

    uncertainty_0 = instance.uncertainty
    ζ = [[uncertainty_0.dev[t]*uncertainty_0.forecast[w][t]*uncertainty_0.error[w][t, List_scenario[scenario]] for t in 1:T] for w in 1:NumWindfarms]

    Netdemand = [instance.Demandbus[b][t] - force * sum(ζ[w][t] for w in 1:NumWindfarms if BusWind[w] == b;init=0) for b in Buses, t in 1:T]
    set_objective_coefficient.(oracle_pb,ν, Netdemand-prod_tot_first_stage)

    D = [0.25* λ * (force* uncertainty_0.dev[t] * uncertainty_0.forecast[w][t])^2 for w in 1:NumWindfarms, t in 1:T]

    set_objective_coefficient.(oracle_pb,ν[instance.BusWind, 1:T], ν[instance.BusWind, 1:T], D)

    start = time()
    optimize!(oracle_pb)
    computation_time = time() - start

    return options.results_second_stage(instance, oracle_pb, prod_tot_first_stage, Pmin, Pmax, computation_time, Netdemand; λ=λ, batch=batch, scenario=scenario, force=force)
end

function second_stage_grb_l1(instance, options, oracle_pb, prod_tot_first_stage, Pmin, Pmax; λ=0.0, batch=1, scenario=1, force=1)

    T= instance.TimeHorizon
    Next=instance.Next
    Buses=1:size(Next)[1] 

    ν=oracle_pb[:ν]
    zp =oracle_pb[:zp]
    zm =oracle_pb[:zm]

    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)

    List_scenario=instance.Training_set[batch]

    uncertainty_0 = instance.uncertainty
    ζ = [[uncertainty_0.dev[t]*uncertainty_0.forecast[w][t]*uncertainty_0.error[w][t, List_scenario[scenario]] for t in 1:T] for w in 1:NumWindfarms]
    Δp = [max(0, 1.96 - uncertainty_0.error[w][t, List_scenario[scenario]]) for w in 1:NumWindfarms, t in 1:T]
    Δm = [max(0, 1.96 + uncertainty_0.error[w][t, List_scenario[scenario]]) for w in 1:NumWindfarms, t in 1:T]

    Netdemand = [instance.Demandbus[b][t] - force * sum(ζ[w][t] for w in 1:NumWindfarms if BusWind[w] == b;init=0) for b in Buses, t in 1:T]
    set_objective_coefficient.(oracle_pb,ν, Netdemand-prod_tot_first_stage)
    set_objective_coefficient.(oracle_pb,zp[1:NumWindfarms, 1:T], - λ*Δp)
    set_objective_coefficient.(oracle_pb,zm[1:NumWindfarms, 1:T], - λ*Δm)

    coeff_p = [-force * uncertainty_0.dev[t]*uncertainty_0.forecast[w][t]*Δp[w,t] for w in 1:NumWindfarms, t in 1:T]
    coeff_m = [force * uncertainty_0.dev[t]*uncertainty_0.forecast[w][t]*Δm[w,t] for w in 1:NumWindfarms, t in 1:T]
    set_objective_coefficient.(oracle_pb,ν[instance.BusWind, 1:T], zp[1:NumWindfarms, 1:T], coeff_p)
    set_objective_coefficient.(oracle_pb,ν[instance.BusWind, 1:T], zm[1:NumWindfarms, 1:T], coeff_m)

    start = time()
    optimize!(oracle_pb)
    computation_time = time() - start

    return options.results_second_stage(instance, oracle_pb, prod_tot_first_stage, Pmin, Pmax, computation_time, Netdemand; λ=λ, batch=batch, scenario=scenario, force=force)
end

"""
Compute uncertainty functions
"""

function compute_uncertainty_SP(instance; x0=0.0, λ=0.0, batch=1, scenario=1, force=1.0)
    T= instance.TimeHorizon
    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)

    uncertainty_0 = instance.uncertainty
    List_scenario=instance.Training_set[batch]
    ζ = [[uncertainty_0.error[w][t, List_scenario[scenario]] for t in 1:T] for w in 1:NumWindfarms]

    return ζ
end

function compute_uncertainty_l2(instance; x0, λ, batch=1, scenario=1, force=1.0)
    T= instance.TimeHorizon
    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)

    uncertainty_0 = instance.uncertainty
    List_scenario=instance.Training_set[batch]
    ζ = [[uncertainty_0.error[w][t, List_scenario[scenario]] for t in 1:T] for w in 1:NumWindfarms]
    ξ = [[ζ[w][t] - (force * uncertainty_0.dev[t]*uncertainty_0.forecast[w][t] * x0[BusWind[w],t]*λ)/(2) for t in 1:T] for w in 1:NumWindfarms]

    return ξ
end

function compute_uncertainty_l1(instance; x0, λ, batch=1, scenario=1, force=1.0)
    T= instance.TimeHorizon
    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)

    uncertainty_0 = instance.uncertainty
    List_scenario=instance.Training_set[batch]

    ζ = [[uncertainty_0.error[w][t, List_scenario[scenario]] for t in 1:T] for w in 1:NumWindfarms]
    ξ = [[0.0 for t in 1:T] for w in 1:NumWindfarms]

    for w in 1:NumWindfarms
        for t in 1:T
            cost_scenario = -force * uncertainty_0.dev[t]*uncertainty_0.forecast[w][t] * ζ[w][t]*x0[BusWind[w],t]
            cost_plus = -force * uncertainty_0.dev[t]*uncertainty_0.forecast[w][t] * 1.96*x0[BusWind[w],t] - λ*abs(1.96- ζ[w][t])
            cost_m = force * uncertainty_0.dev[t]*uncertainty_0.forecast[w][t] * 1.96*x0[BusWind[w],t] - λ*abs(-1.96- ζ[w][t])
            if cost_plus > cost_scenario && cost_plus > cost_m
                ξ[w][t] = 1.96
            elseif cost_m > cost_scenario && cost_m > cost_plus
                ξ[w][t] = -1.96
            else
                ξ[w][t] = ζ[w][t]
            end
        end
    end

    return ξ
end

"""
return results functions
"""

function results_second_stage_SP(instance, oracle_pb, prod_tot_first_stage, Pmin, Pmax, computation_time, Netdemand; λ=0.0, batch=1, scenario=1, force=1)
    
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    N2 = N - N1
    thermal_units_2=values(instance.Thermalunits)[N1+1:N]
    Next=instance.Next
    Buses=1:size(Next)[1] 
    
    status = termination_status(oracle_pb)!=MOI.DUAL_INFEASIBLE

    νₖ=convert(Matrix{Float64}, JuMP.value.(oracle_pb[:ν]))

    network_cost_val = value.(oracle_pb[:network_cost])
    βₖ2=[sum(νₖ[b,t]*Netdemand[b,t] for b in Buses) + network_cost_val[t] for t in 1:T]

    βₖ= βₖ2

    mumax, mumin = compute_mu_max_min(N2, T, νₖ, thermal_units_2)

    obj = [sum(mumin[i, t]*Pmin[i,t] - mumax[i, t]*Pmax[i,t] for i in 1:N2) + βₖ2[t] - sum(νₖ[b,t] * prod_tot_first_stage[b,t] for b in Buses) for t in 1:T]

    if 100*abs((sum(obj) - objective_value(oracle_pb))/objective_value(oracle_pb)) > 1e-4
        println(("PROBLEM :", sum(obj), objective_value(oracle_pb), abs(sum(obj) - objective_value(oracle_pb)), obj))
    end

    return oracleResults(obj, βₖ, βₖ2, zeros(T), mumax, mumin, νₖ, status, computation_time)
end

function results_second_stage_dca_l2(instance, oracle_pb, prod_tot_first_stage, Pmin, Pmax, computation_time, Lindemand; λ=0.0, batch=1, scenario=1, force=1)
    
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    N2 = N - N1
    thermal_units_2=values(instance.Thermalunits)[N1+1:N]
    Next=instance.Next
    Buses=1:size(Next)[1] 
    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)

    νₖ=convert(Matrix{Float64}, JuMP.value.(oracle_pb[:ν]))

    network_cost_val = value.(oracle_pb[:network_cost])

    mumax, mumin = compute_mu_max_min(N2, T, νₖ, thermal_units_2)

    List_scenario=instance.Training_set[batch]
    uncertainty_0 = instance.uncertainty
    ζ = [[uncertainty_0.dev[t]*uncertainty_0.forecast[w][t]*uncertainty_0.error[w][t, List_scenario[scenario]] for t in 1:T] for w in 1:NumWindfarms]
    Netdemand = [instance.Demandbus[b][t] - force * sum(ζ[w][t] for w in 1:NumWindfarms if BusWind[w] == b;init=0) for b in Buses, t in 1:T]

    βₖ2=[sum(νₖ[b,t]*Netdemand[b,t] for b in Buses) + network_cost_val[t] for t in 1:T]
    βₖ3=[sum(νₖ[b,t]*Lindemand[b,t] for b in Buses) + network_cost_val[t] for t in 1:T]

    γₖ=[0.25*force^2*sum(νₖ[b, t]^2 * uncertainty_0.dev[t]^2 * uncertainty_0.forecast[w][t]^2 for b in instance.BusWind for w in 1:NumWindfarms if BusWind[w] == b;init=0) for t in 1:T]

    objlin = [sum(mumin[i, t]*Pmin[i,t] - mumax[i, t]*Pmax[i,t] for i in 1:N2) + βₖ3[t] - sum(νₖ[b,t] * prod_tot_first_stage[b,t] for b in Buses) for t in 1:T]

    obj = [sum(mumin[i, t]*Pmin[i,t] - mumax[i, t]*Pmax[i,t] for i in 1:N2) + βₖ2[t] + γₖ[t]*λ - sum(νₖ[b,t] * prod_tot_first_stage[b,t] for b in Buses) for t in 1:T]

    if 100*abs((sum(objlin) - objective_value(oracle_pb))/objective_value(oracle_pb)) > 1e-4
        println(("PROBLEM :", sum(objlin), objective_value(oracle_pb), abs(sum(objlin) - objective_value(oracle_pb)), objlin))
    end

    status = true

    return oracleResults(obj, βₖ2, βₖ2, γₖ, mumax, mumin, νₖ, status, computation_time)
end

function results_second_stage_dca_l1(instance, oracle_pb, prod_tot_first_stage, Pmin, Pmax, computation_time, Lindemand; λ=0.0, batch=1, scenario=1, force=1)
    
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    N2 = N - N1
    thermal_units_2=values(instance.Thermalunits)[N1+1:N]
    Next=instance.Next
    Buses=1:size(Next)[1] 
    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)

    uncertainty_0 = instance.uncertainty
    List_scenario=instance.Training_set[batch]

    ζ = [[uncertainty_0.error[w][t, List_scenario[scenario]] for t in 1:T] for w in 1:NumWindfarms]

    uncertainty = [(instance.Demandbus[BusWind[w]][t] - Lindemand[BusWind[w],t])/(force*uncertainty_0.dev[t]*uncertainty_0.forecast[w][t]) for w in 1:NumWindfarms, t in 1:T]

    νₖ=convert(Matrix{Float64}, JuMP.value.(oracle_pb[:ν]))

    network_cost_val = value.(oracle_pb[:network_cost])

    βₖ2=[sum(νₖ[b,t]*Lindemand[b,t] for b in Buses) + network_cost_val[t] for t in 1:T]

    βₖ= βₖ2

    γₖ=-[sum(abs(uncertainty[w,t]-ζ[w][t]) for w in 1:NumWindfarms;init=0) for t in 1:T]

    mumax, mumin = compute_mu_max_min(N2, T, νₖ, thermal_units_2)

    objlin = [sum(mumin[i, t]*Pmin[i,t] - mumax[i, t]*Pmax[i,t] for i in 1:N2) + βₖ2[t] - sum(νₖ[b,t] * prod_tot_first_stage[b,t] for b in Buses) for t in 1:T]

    obj = [sum(mumin[i, t]*Pmin[i,t] - mumax[i, t]*Pmax[i,t] for i in 1:N2) + βₖ2[t] - sum(νₖ[b,t] * prod_tot_first_stage[b,t] for b in Buses) + γₖ[t]*λ for t in 1:T]

    if 100*abs((sum(objlin) - objective_value(oracle_pb))/objective_value(oracle_pb)) > 1e-4 && objective_value(oracle_pb)<= 1e8
        println(("PROBLEM :", sum(objlin), objective_value(oracle_pb), abs(sum(objlin) - objective_value(oracle_pb)), objlin))
    end

    status = true

    return oracleResults(obj, βₖ, βₖ2, γₖ, mumax, mumin, νₖ, status, computation_time)
end

function results_second_stage_grb_l2(instance, oracle_pb, prod_tot_first_stage, Pmin, Pmax, computation_time, Netdemand; λ=0.0, batch=1, scenario=1, force=1)

    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    N2 = N - N1
    thermal_units_2=values(instance.Thermalunits)[N1+1:N]
    Next=instance.Next
    Buses=1:size(Next)[1] 
    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)
    uncertainty_0 = instance.uncertainty

    νₖ=convert(Matrix{Float64}, JuMP.value.(oracle_pb[:ν]))

    network_cost_val = value.(oracle_pb[:network_cost])
    βₖ2=[sum(νₖ[b,t]*Netdemand[b,t] for b in Buses) + network_cost_val[t] for t in 1:T]

    βₖ= βₖ2

    mumax, mumin = compute_mu_max_min(N2, T, νₖ, thermal_units_2)

    γₖ=[0.25*force^2*sum(νₖ[b, t]^2 * uncertainty_0.dev[t]^2 * uncertainty_0.forecast[w][t]^2 for b in instance.BusWind for w in 1:NumWindfarms if BusWind[w] == b;init=0) for t in 1:T]

    obj = [sum(mumin[i, t]*Pmin[i,t] - mumax[i, t]*Pmax[i,t] for i in 1:N2) + βₖ2[t] - sum(νₖ[b,t] * prod_tot_first_stage[b,t] for b in Buses) + γₖ[t]*λ for t in 1:T]

    if 100*abs((sum(obj) - objective_value(oracle_pb))/objective_value(oracle_pb)) > 1e-4 && objective_value(oracle_pb)<= 1e8
        println(("PROBLEM :", sum(obj), objective_value(oracle_pb), abs(sum(obj) - objective_value(oracle_pb)), obj))
    end

    status = true

    return oracleResults(obj, βₖ, βₖ2, γₖ, mumax, mumin, νₖ, status, computation_time)
end

function results_second_stage_grb_l1(instance, oracle_pb, prod_tot_first_stage, Pmin, Pmax, computation_time, Netdemand; λ=0.0, batch=1, scenario=1, force=1)

    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    N2 = N - N1
    thermal_units_2=values(instance.Thermalunits)[N1+1:N]
    Next=instance.Next
    Buses=1:size(Next)[1] 
    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)
    uncertainty_0 = instance.uncertainty

    List_scenario=instance.Training_set[batch]
    Δp = [1.96 - uncertainty_0.error[w][t, List_scenario[scenario]] for w in 1:NumWindfarms, t in 1:T]
    Δm = [1.96 + uncertainty_0.error[w][t, List_scenario[scenario]] for w in 1:NumWindfarms, t in 1:T]

    νₖ=convert(Matrix{Float64}, JuMP.value.(oracle_pb[:ν]))
    zpₖ=convert(Matrix{Float64}, JuMP.value.(oracle_pb[:zp]))
    zmₖ=convert(Matrix{Float64}, JuMP.value.(oracle_pb[:zm]))

    network_cost_val = value.(oracle_pb[:network_cost])
    βₖ2=[sum(νₖ[b,t]*Netdemand[b,t] for b in Buses) + network_cost_val[t] - force * sum(uncertainty_0.dev[t]*uncertainty_0.forecast[w][t]*νₖ[b,t] * (Δp[w,t]*zpₖ[w,t] - Δm[w,t]*zmₖ[w,t]) for b in Buses for w in 1:NumWindfarms if BusWind[w] == b; init=0) for t in 1:T]

    βₖ= βₖ2

    mumax, mumin = compute_mu_max_min(N2, T, νₖ, thermal_units_2)

    γₖ=-[sum(Δp[w,t]*zpₖ[w,t] + Δm[w,t]*zmₖ[w,t] for w in 1:NumWindfarms;init=0) for t in 1:T]

    obj = [sum(mumin[i, t]*Pmin[i,t] - mumax[i, t]*Pmax[i,t] for i in 1:N2) + βₖ2[t] - sum(νₖ[b,t] * prod_tot_first_stage[b,t] for b in Buses) + γₖ[t]*λ for t in 1:T]

    if 100*abs((sum(obj) - objective_value(oracle_pb))/objective_value(oracle_pb)) > 1e-4 && objective_value(oracle_pb)<= 1e8
        println(("PROBLEM :", sum(obj), objective_value(oracle_pb), abs(sum(obj) - objective_value(oracle_pb)), obj))
    end

    status = true

    return oracleResults(obj, βₖ, βₖ2, γₖ, mumax, mumin, νₖ, status, computation_time)
end

function compute_mu_max_min(N2, T, νₖ, thermal_units_2)
    mumax=zeros(Float64, N2, T)
    mumin=zeros(Float64, N2, T)

    for i in 1:N2
        for t in 1:T
            mumax[i,t] = max(0,νₖ[thermal_units_2[i].Bus, t] - thermal_units_2[i].LinearTerm)
            mumin[i,t] = max(0,thermal_units_2[i].LinearTerm - νₖ[thermal_units_2[i].Bus, t])
        end
    end
    return mumax, mumin
end


function DCAlgo(instance, options, oracle_pb, prod_tot_first_stage, Pmin, Pmax; λ=0.0, batch=1, scenario=1, force=1.0)
    k=1
    T= instance.TimeHorizon
    Next=instance.Next
    Buses=1:size(Next)[1]
    priceref=zeros(Buses, 1:T)
    uncertainty = options.compute_uncertainty(instance; batch=batch, scenario=scenario, force=force, x0 = priceref, λ=λ)
    resDCAk=second_stage_DCA(instance, options, oracle_pb, uncertainty, prod_tot_first_stage, Pmin, Pmax, λ; batch=batch, scenario=scenario, force=force)
    cost=zeros(T)

    comp_time=resDCAk.computation_time

    while sum(abs.(resDCAk.objective_value-cost))>=1e-3 && k<=10
        cost=resDCAk.objective_value
        priceref=resDCAk.ν
        
        uncertainty = options.compute_uncertainty(instance; batch=batch, scenario=scenario, force=force, x0 = priceref, λ=λ)
        resDCAk=second_stage_DCA(instance, options, oracle_pb, uncertainty, prod_tot_first_stage, Pmin, Pmax, λ; batch=batch, scenario=scenario, force=force)
        comp_time+=resDCAk.computation_time
        k+=1
    end
    resDCAk.computation_time=comp_time
    return resDCAk
end

function second_stage_RO(instance, options, oracle_pb, prod_tot_first_stage, solution_x::Vector{Matrix{Float64}}; Γ=1.0, force=1.0)

    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    N2 = N - N1
    thermal_units_2=values(instance.Thermalunits)[N1+1:N]
    Next=instance.Next
    Buses=1:size(Next)[1] 

    ν=oracle_pb[:ν]
    ζ=oracle_pb[:ζ]
    γ=oracle_pb[:γ]

    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)

    Pmin, Pmax=get_limit_power_solution(instance, solution_x)

    uncertainty_0 = instance.uncertainty

    @objective(oracle_pb, Max, sum(oracle_pb[:μₘᵢₙ][i,t]*Pmin[i,t] - oracle_pb[:μₘₐₓ][i,t]*Pmax[i,t] for i in 1:N2 for t in 1:T)+sum(oracle_pb[:network_cost][t] for t in 1:T)
    + sum((instance.Demandbus[b][t]-prod_tot_first_stage[b,t])*ν[b,t] for b in Buses, t in 1:T) + sum(force*1.96*uncertainty_0.dev[t]*uncertainty_0.forecast[w][t]*ζ[BusWind[w],t] for w in 1:NumWindfarms, t in 1:T))
  

    start = time()
    optimize!(oracle_pb)
    computation_time = time() - start

    status = termination_status(oracle_pb)!=MOI.DUAL_INFEASIBLE

    νₖ=JuMP.value.(ν)
    γₖ=JuMP.value.(γ)
    ζₖ=JuMP.value.(ζ)
    
    d_up=[(b,t) for b in BusWind, t in 1:T if γₖ[b,t]>1e-5 && νₖ[b,t]>1e-5]
    d_down=[(b,t) for b in BusWind, t in 1:T if γₖ[b,t]>1e-5 && νₖ[b,t]<1e-5]

    # return objective_value(oracle_pb), d_up, d_down, computation_time
    mumax, mumin = compute_mu_max_min(N2, T, νₖ, thermal_units_2)
    network_cost_val = value.(oracle_pb[:network_cost])
    βₖ2=[sum(νₖ[b,t]*instance.Demandbus[b][t] for b in Buses)+sum(force*1.96*uncertainty_0.dev[t]*uncertainty_0.forecast[w][t]*ζₖ[BusWind[w],t] for w in 1:NumWindfarms) + network_cost_val[t] for t in 1:T]

    βₖ= βₖ2
    obj = [sum(mumin[i, t]*Pmin[i,t] - mumax[i, t]*Pmax[i,t] for i in 1:N2) + βₖ2[t] - sum(νₖ[b,t] * prod_tot_first_stage[b,t] for b in Buses) for t in 1:T]
    
    if 100*abs((sum(obj) - objective_value(oracle_pb))/objective_value(oracle_pb)) > 1e-4
        println(("PROBLEM :", sum(obj), objective_value(oracle_pb), abs(sum(obj) - objective_value(oracle_pb)), obj))
    end
    return oracleResults(obj, βₖ, βₖ2, zeros(T), mumax, mumin, νₖ, status, computation_time)
end

function compute_uncertainty_RO(instance; x0, Γ=1, force=1.0)
    T= instance.TimeHorizon
    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)   
    uncertainty_0 = instance.uncertainty
    
    ξ = [zeros(T) for w in 1:NumWindfarms]
    cost = [[force * uncertainty_0.dev[t]*uncertainty_0.forecast[w][t]*abs(x0[BusWind[w],t]) for w in 1:NumWindfarms] for t in 1:T]

    y0 = x0[instance.BusWind, :]

    for t in 1:T
        idx_lin = partialsortperm(cost[t], 1:Γ; rev=true)
        for idx in idx_lin
            if y0[idx, t] < 0
                ξ[idx][t] = 1.96
            else
                ξ[idx][t] = -1.96
            end
        end
    end

    return ξ
end

function DCAlgo_RO(instance, options, oracle_pb, prod_tot_first_stage, solution_x::Vector{Matrix{Float64}}; Γ=1.0,force=1.0)
    k=1
    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    N2 = N - N1
    Next=instance.Next
    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)  
    Buses=1:size(Next)[1]
    priceref=zeros(Buses, 1:T)
    
    Pmin, Pmax=get_limit_power_solution(instance, solution_x)

    @objective(oracle_pb, Max, sum(oracle_pb[:μₘᵢₙ][i,t]*Pmin[i,t] - oracle_pb[:μₘₐₓ][i,t]*Pmax[i,t] for i in 1:N2 for t in 1:T)+sum(oracle_pb[:network_cost][t] for t in 1:T))

    uncertainty = options.compute_uncertainty(instance; x0=priceref,Γ=Γ, force=force)
    resDCAk=second_stage_DCA_RO(instance, options, oracle_pb, uncertainty, prod_tot_first_stage, Pmin, Pmax; force=force)
    cost=zeros(T)

    comp_time=resDCAk.computation_time

    while sum(abs.(resDCAk.objective_value-cost))>=1e-3 && k<=10
        # println((k, sum(cost), sum(resDCAk.objective_value), abs.(resDCAk.objective_value-cost)))
        cost=resDCAk.objective_value
        priceref=resDCAk.ν
        
        uncertainty = options.compute_uncertainty(instance; x0=priceref,Γ=Γ, force=force)
        resDCAk=second_stage_DCA_RO(instance, options, oracle_pb, uncertainty, prod_tot_first_stage, Pmin, Pmax; force=force)
        comp_time+=resDCAk.computation_time
        k+=1
    end
    resDCAk.computation_time=comp_time

    d_up=[(BusWind[w],t) for w in 1:NumWindfarms, t in 1:T if uncertainty[w][t]<-1e-5]
    d_down=[(BusWind[w],t) for w in 1:NumWindfarms, t in 1:T if uncertainty[w][t]>1e-5]

    if options == BD_RO_DCA
        return resDCAk
    end

    return objective_value(oracle_pb), d_up, d_down, comp_time
end

function second_stage_DCA_RO(instance, options, oracle_pb, uncertainty, prod_tot_first_stage, Pmin, Pmax;  force=force)

    computation_time, Netdemand = solve_second_stage(instance, uncertainty, oracle_pb, prod_tot_first_stage; force=force)

    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    N2 = N - N1
    thermal_units_2=values(instance.Thermalunits)[N1+1:N]
    Next=instance.Next
    Buses=1:size(Next)[1] 
    
    status = termination_status(oracle_pb)!=MOI.DUAL_INFEASIBLE

    νₖ=convert(Matrix{Float64}, JuMP.value.(oracle_pb[:ν]))

    network_cost_val = value.(oracle_pb[:network_cost])
    βₖ2=[sum(νₖ[b,t]*Netdemand[b,t] for b in Buses) + network_cost_val[t] for t in 1:T]

    βₖ= βₖ2

    mumax, mumin = compute_mu_max_min(N2, T, νₖ, thermal_units_2)

    obj = [sum(mumin[i, t]*Pmin[i,t] - mumax[i, t]*Pmax[i,t] for i in 1:N2) + βₖ2[t] - sum(νₖ[b,t] * prod_tot_first_stage[b,t] for b in Buses) for t in 1:T]

    if 100*abs((sum(obj) - objective_value(oracle_pb))/objective_value(oracle_pb)) > 1e-4
        println(("PROBLEM :", sum(obj), objective_value(oracle_pb), abs(sum(obj) - objective_value(oracle_pb)), obj))
    end

    return oracleResults(obj, βₖ, βₖ2, zeros(T), mumax, mumin, νₖ, status, computation_time)
end

function second_stage_RO_2(instance, options, oracle_pb, prod_tot_first_stage, solution_x::Vector{Matrix{Float64}}; Γ=1.0, force=1.0)

    T= instance.TimeHorizon
    N=instance.N
    N1=instance.N1
    N2 = N - N1
    thermal_units_2=values(instance.Thermalunits)[N1+1:N]
    Next=instance.Next
    Buses=1:size(Next)[1] 

    ν=oracle_pb[:ν]
    δp=oracle_pb[:δp]
    δm=oracle_pb[:δm]
    γ=oracle_pb[:γ]

    BusWind=instance.BusWind
    NumWindfarms=length(BusWind)

    Pmin, Pmax=get_limit_power_solution(instance, solution_x)

    uncertainty_0 = instance.uncertainty

    @objective(oracle_pb, Max, sum(oracle_pb[:μₘᵢₙ][i,t]*Pmin[i,t] - oracle_pb[:μₘₐₓ][i,t]*Pmax[i,t] for i in 1:N2 for t in 1:T)+sum(oracle_pb[:network_cost][t] for t in 1:T)
    + sum((instance.Demandbus[b][t]-prod_tot_first_stage[b,t])*ν[b,t] for b in Buses, t in 1:T) + sum(force*1.96*uncertainty_0.dev[t]*uncertainty_0.forecast[w][t]*ν[BusWind[w],t]*(δp[BusWind[w],t]-δm[BusWind[w],t]) for w in 1:NumWindfarms, t in 1:T))
  

    start = time()
    optimize!(oracle_pb)
    computation_time = time() - start

    println(computation_time)


    δpₖ=JuMP.value.(δp)
    δmₖ=JuMP.value.(δm)

    d_up=[(b,t) for b in BusWind, t in 1:T if δpₖ[b,t]>1e-5]
    d_down=[(b,t) for b in BusWind, t in 1:T if δmₖ[b,t]>1e-5]

    return objective_value(oracle_pb), d_up, d_down, computation_time
end