module UC

using JuMP
using Gurobi
using MosekTools
using Statistics
using NCDatasets
using Plots
using XLSX
using CSV
using DataFrames
using Random
using Distributions

const SHEDDING_COST=1000.0
const CURTAILEMENT_COST=1000.0


GUROBI_ENV = nothing
"""
    set_solver_Gurobi()

Set the Gurobi solver.
"""
function set_solver_Gurobi()
    try
        @eval using Gurobi
        global SOLVER = "Gurobi"
        if isnothing(GUROBI_ENV)
            global GUROBI_ENV = Gurobi.Env()
        end
    catch
        @warn("Unable to use Gurobi. Is it properly installed? Have you called `pkg> add Gurobi`?\nDefaulting to SCIP")
    end
end

function set_solver_Mosek()
    try
        @eval using MosekTools
        global SOLVER = "Mosek"
    catch
        @warn("Unable to use Mosek. Is it properly installed? Have you called `pkg> add Mosek`?\nDefaulting to SCIP")
    end
end

include("Struct/Instance.jl")
include("Struct/tools.jl")
include("Struct/Thermal_unit.jl")
include("Struct/parsing.jl")
include("Optimizer/initialisation_Benders.jl")
include("Optimizer/second_stage.jl")
include("Optimizer/add_cut_DRO.jl")
include("Optimizer/add_cut_SP.jl")
include("Optimizer/add_cut_AVAR.jl")
include("Optimizer/add_cut_KL.jl")
include("Optimizer/add_cut_RO.jl")
include("Optimizer/add_cut_moment.jl")
include("Optimizer/add_cut_DRO_budget.jl")
include("Optimizer/benders.jl")
include("Optimizer/options.jl")
include("Optimizer/out_of_sample.jl")

export set_solver_Gurobi, set_solver_Mosek

end