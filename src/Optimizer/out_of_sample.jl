function test_out_of_sample(instance, solution, method; force=1.0)
    computational_time = solution[2]
    UB = solution[3]
    LB= solution[4]
    first_stage_solution = solution[5]
    S_train = solution[6]
    batch_train = solution[7]
    gap_train = solution[9]
    radius_train = solution[11]
    Γ_train = solution[12]


    size_test_set=length(instance.Test_set)
    res_obj=[0.1 for s in 1:size_test_set]
    res_fake1=[0.1 for s in 1:size_test_set]
    res_fake2=[0.1 for s in 1:size_test_set]
    flow_vals = []
    power_shedding_vals = []
    T= instance.TimeHorizon
    N=instance.N
    thermal_units_name=keys(instance.Thermalunits)
    thermal_units=instance.Thermalunits

    model = initializeJuMPModel()
    set_silent(model)

    @variable(model, is_on[unit in thermal_units_name, t in 0:T], Bin)
    JuMP.fix.(is_on, first_stage_solution.is_on; force=true)
    @variable(model, start_up[unit in thermal_units_name, t in 1:T], Bin)
    JuMP.fix.(start_up, first_stage_solution.start_up; force=true)
    @variable(model, start_down[unit in thermal_units_name, t in 1:T], Bin)
    JuMP.fix.(start_down, first_stage_solution.start_down; force=true)

    @variable(model, thermal_fuel_cost>=0)
    @variable(model, thermal_fixed_cost>=0)
    @variable(model, thermal_cost>=0)

    Next=instance.Next
    Buses=1:size(Next)[1]

    @variable(model, power_shedding[b in Buses, t in 0:T]>=0)
    @variable(model, power_curtailement[b in Buses, t in 0:T]>=0)

    thermal_unit_commit_constraints(model, instance)

    @constraint(model,  thermal_cost>=thermal_fuel_cost)

    Numlines=length(instance.Lines)

    @variable(model, θ[b in Buses, t in 1:T])
    @variable(model, flow[l in 1:Numlines, t in 1:T])
    
    Lines=instance.Lines

    @constraint(model, [line in Lines, t in 1:T], flow[line.id,t]<=line.Fmax)
    @constraint(model, [line in Lines, t in 1:T], flow[line.id,t]>=-line.Fmax)
    @constraint(model, [line in Lines, t in 1:T], flow[line.id,t]==line.B12*(θ[line.b1,t]-θ[line.b2,t]))


    @variable(model, power[unit in thermal_units_name, t in 0:T]>=0)
    N1=instance.N1
    @constraint(model, [i in 1:N1, t in 1:T], power[i, t] == first_stage_solution.power_ref[i, t+1])

    @constraint(model,  [unit in thermal_units; unit.InitialPower!=nothing], power[unit.name, 0]==unit.InitialPower)
    @constraint(model,  [unit in thermal_units, t in 0:T], power[unit.name, t]>=unit.MinPower*is_on[unit.name, t])
    @constraint(model,  [unit in thermal_units, t in 0:T], power[unit.name, t]<=unit.MaxPower*is_on[unit.name, t])

    @variable(model, ξ[b in Buses, t in 1:T])

    @constraint(model, thermal_fuel_cost>= sum(unit.LinearTerm*power[unit.name, t] for unit in thermal_units for t in 1:T) + sum(SHEDDING_COST*power_shedding[b,t]+CURTAILEMENT_COST*power_curtailement[b,t] for b in Buses for t in 1:T))

    @constraint(model,  [t in 1:T, b in Buses], sum(power[unit.name, t] for unit in thermal_units if unit.Bus==b)+power_shedding[b,t]-power_curtailement[b,t] + sum(flow[line.id,t] for line in Lines if line.b2==b) - sum(flow[line.id,t] for line in Lines if line.b1==b)==instance.Demandbus[b][t]-force * ξ[b,t])

    @objective(model, Min, thermal_fixed_cost+thermal_cost)

    k=0
    BusWind=instance.BusWind

    uncertainty_0 = instance.uncertainty
    for s in instance.Test_set
        k+=1

        # if k % 100 == 0
        #     println("OOS scenario ", k, "/", size_test_set)
        # end

        ξ_s = [zeros(T) for b in Buses]
        j=0
        for b in BusWind
            j+=1
            ξ_s[b] = [uncertainty_0.dev[t]*uncertainty_0.forecast[j][t]*uncertainty_0.error[j][t, s] for t in 1:T]
        end

        for b in Buses, t in 1:T
            JuMP.fix(ξ[b,t], ξ_s[b][t])
        end

        optimize!(model)
        
        feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT

        if feasibleSolutionFound
            res_obj[k]=JuMP.objective_value(model)
            res_fake1[k]=sum(value.(model[:power_shedding]))
            res_fake2[k]=sum(value.(model[:power_curtailement]))
            push!(flow_vals, value.(model[:flow]))
            push!(power_shedding_vals, value.(model[:power_shedding]))
        else
            println((s, "ERROR"))
            return instance.name, 0.0, 0.0, 0.0, 0.0, res_obj, res_fake1, res_fake2
            break
        end
    end
    fixed_cost= value.(thermal_fixed_cost)
    mean_obj=mean(res_obj)
    max_obj=maximum(res_obj)
    mean_fake1=mean(res_fake1)
    mean_fake2=mean(res_fake2)

    OOS = [mean_obj, mean_fake1, mean_fake2, max_obj]

    return return_oos_solution(method, instance, computational_time, LB, UB, S_train, batch_train, gap_train, force, radius_train, Γ_train, OOS)
end

function return_oos_solution(method::String, instance::Instance, computational_time::Float64, LB::Float64, UB::Float64, S_train, batch_train, gap_train, force, radius_train, Γ_train, OOS)
    name_csv = "$(instance.name)_"*string(S_train)*"_"*string(batch_train)*"_"*string(force)*"_"*string(gap_train)*"_"*string(radius_train)*"_"*string(Γ_train)
    # write results to CSV (long format: one datum per line)
    try
        df = DataFrame(
            metric = String[],
            value = String[],
        )
        # scalars
        push!(df, ("name", "$(instance.name)"))
        push!(df, ("S", string(S_train)))
        push!(df, ("batch", string(batch_train)))
        push!(df, ("force", string(force)))
        push!(df, ("Time", string(computational_time)))
        push!(df, ("LB", string(LB)))
        push!(df, ("UB", string(UB)))
        push!(df, ("radius", string(radius_train)))
        push!(df, ("Γ", string(Γ_train)))
        push!(df, ("OOS", string(OOS[1])))
        filepath = joinpath(pwd(), "results/"*method, name_csv*".csv")
        CSV.write(filepath, df)
    catch e
        @warn("Failed to write results CSV", error = e)
    end
    return name_csv, OOS
end