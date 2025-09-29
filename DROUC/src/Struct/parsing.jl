function parse_IEEE_JEAS(folder, optimizer; N1=35, NumWind=91)
    """
    Parse the IEEE 118-bus instnace
    """
    TimeHorizon= 24
    syst = "Data/"*folder
    generators = CSV.read(joinpath(pwd(), syst, "generators.csv"), DataFrame; header=false)
    NumberUnits= parse(Int64,generators[end,1]) 
    N=NumberUnits
    name_instance="IEEE"*string(NumberUnits)
    Thermal_units=Vector{ThermalUnit}(undef, N)
    for i in 2:NumberUnits+1
        unit_name=i-1
        Bus = parse(Int64,generators[i,2])
        ConstTerm = parse(Float64,generators[i,3])
        LinearTerm = parse(Float64,generators[i,4])
        MaxPower = parse(Float64,generators[i,6])
        MinPower = parse(Float64,generators[i,7])
        DeltaRampUp = parse(Float64,generators[i,14])
        DeltaRampDown = parse(Float64,generators[i,14])
        StartUpCost = parse(Float64,generators[i,15])
        StartDownCost = 0.0*parse(Float64,generators[i,15])
        MinUpTime=parse(Int64,generators[i,13]) 
        MinDownTime=parse(Int64,generators[i,12])
        QuadTerm =0.0        
        InitialPower=parse(Float64,generators[i,11])
        InitUpDownTime =parse(Int64,generators[i,10])
        intervals=set_intervals(MinUpTime, InitUpDownTime, InitialPower, MaxPower, MinPower, DeltaRampDown)
        Tup=Int64(1+floor((MaxPower-MinPower)/DeltaRampUp))
        Tdown=Int64(1+floor((MaxPower-MinPower)/DeltaRampDown))
        if unit_name <= N1
            unitgroup = 1
        else
            unitgroup = 2
        end
        unit=ThermalUnit(unit_name, Bus, unitgroup, MinPower, MaxPower, DeltaRampUp, DeltaRampDown, QuadTerm, StartUpCost, StartDownCost, LinearTerm, ConstTerm, InitialPower, InitUpDownTime, MinUpTime, MinDownTime, intervals, Tup, Tdown)
        Thermal_units[unit_name]=unit
    end

    maximum_load = CSV.read(joinpath(pwd(), syst, "maximum_load.csv"), DataFrame; header=false)
    Numbus=maximum_load[end,1]
    Buses=1:Numbus
    load_distribution_profile = CSV.read(joinpath(pwd(), syst, "load_distribution_profile.csv"), DataFrame; header=false)[:,2]/100
    Demandbus=[maximum_load[b,2]*load_distribution_profile for b in Buses]
    df_lines = CSV.read(joinpath(pwd(), syst, "lines.csv"), DataFrame; header=false)
    Numlines=parse(Int64,df_lines[end,1])
    Lines=Dict()
    Next=[[] for b in Buses]
    for i in 2:Numlines+1
        b1 = parse(Int64, df_lines[i,2])
        b2 = parse(Int64, df_lines[i,3])
        fmax = parse(Float64, df_lines[i,6])
        X = 1/parse(Float64, df_lines[i,5])
        Lines[b1, b2]=Line(b1,b2, fmax, X)
        Lines[b2, b1]=Line(b2,b1, fmax, X)
        push!(Next[b1], b2)
        push!(Next[b2], b1)
    end

    Windfarms=1:NumWind
    BusWind=[b for b in Buses if sum(Demandbus[b])>=1][Windfarms]
    Nb_total_scenario=2000

    WGscenario=[zeros(24,Nb_total_scenario) for b in BusWind]
    dev = [0.0 for t in 1:TimeHorizon]

    k=0

    scenarios_ref = simulate_scenarios(Nb_total_scenario*NumWind)
    dev = [std(scenarios_ref[t,:]) for t in 1:TimeHorizon]

    for b in BusWind
        k+=1
        WGscenario[k]=scenarios_ref[:,1+(k-1)*Nb_total_scenario:k*Nb_total_scenario]./dev
    end

    forecast = [[Demandbus[b][t] for t in 1:TimeHorizon] for b in BusWind]
    uncertainty = Uncertainty(dev, forecast, WGscenario)

    Num_batch=10
    Smax=100
    Training_set=[[(batch-1)*Smax+day for day in 1:Smax] for batch in 1:Num_batch]
    Test_set=Num_batch*Smax+1:Nb_total_scenario
    push!(Training_set, [day for day in Test_set])

    return Instance(name_instance, TimeHorizon, N, N1, uncertainty, Thermal_units, Lines, Next, Demandbus, BusWind, Training_set, Test_set, optimizer)
end