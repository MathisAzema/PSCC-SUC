struct ThermalUnit
    name::Int64
    Bus::Int64
    Group::Int64
    MinPower::Float64  
    MaxPower::Float64  
    DeltaRampUp::Float64  
    DeltaRampDown::Float64 
    QuadTerm::Float64  
    StartUpCost::Float64
    StartDownCost::Float64  
    LinearTerm::Float64
    ConstTerm::Float64 
    InitialPower::Float64
    InitUpDownTime::Int64
    MinUpTime::Int64  
    MinDownTime::Int64
    intervals::Vector{Vector{Int64}}
    Tup::Int64
    Tdown::Int64
end

struct Line
    b1::Int64
    b2::Int64 
    Fmax::Float64  
    B12::Float64 
end

struct Uncertainty
    dev ::Vector{Float64}
    forecast::Vector{Vector{Float64}}
    error::Vector{Matrix{Float64}}
end

struct Instance
    name::String
    TimeHorizon::Int64 
    N::Int64
    N1::Int64
    uncertainty::Uncertainty
    Thermalunits::Vector{ThermalUnit}
    Lines::Dict{Tuple{Int64, Int64}, Line}
    Next::Vector{Vector{Int64}}
    Demandbus::Vector{Vector{Float64}}
    BusWind::Vector{Int64}
    Training_set::Vector{Vector{Int64}}
    Test_set::UnitRange{Int64}
    optimizer::Any
end