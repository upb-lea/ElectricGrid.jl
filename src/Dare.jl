module Dare

using ArnoldiMethod #*
using Combinatorics
using CSV
using CUDA
using DataStructures
using ControlSystemsBase
using DifferentialEquations
using Distributions
using DSP
using FFTW
using LinearAlgebra
using LinearMaps #*
using JuMP
using NearestNeighbors #*
using NonNegLeastSquares #*
using PlotlyJS
using Flux
using Graphs
using GraphPlot
using IJulia
import Ipopt
using IntervalSets
using PlotlyJS
using Random
using ReinforcementLearning
using StableRNGs
using SpecialFunctions
using StatsBase

#export create_setup, Classical_Policy, create_agent_ddpg, Source_Initialiser, MultiAgentGridController, data_hook, plot_hook_results, plot_best_results, NodeConstructor, dare_setup, SimEnv

include("./Power_System_Theory.jl")
include("./nodeconstructor.jl")
include("./custom_control.jl")
include("./pv_module.jl")
include("./env.jl")
include("./agent_ddpg.jl")
include("./Classical_Control.jl")
include("./MultiAgentGridController.jl")
include("./plotting.jl")
include("./data_hook.jl")
include("./Dare_Wrapper.jl")
include("./Dare_Logger.jl")
include("./Kernel_Machine.jl")
include("./Machine_Dynamics.jl")
include("./Dif_Map.jl")
include("./Dynamical_Systems.jl")
include("./Complexity.jl")

#code to export all, taken from https://discourse.julialang.org/t/exportall/4970/18
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end # module