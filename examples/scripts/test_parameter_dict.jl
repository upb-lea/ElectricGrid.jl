using DrWatson
@quickactivate "JEG"

using ReinforcementLearning
using IntervalSets
using LinearAlgebra
using ControlSystems
using CUDA
using PlotlyJS


include(srcdir("node_constructor.jl"))
include(srcdir("env.jl"));
include(srcdir("sin_policy.jl"))
include(srcdir("data_hook.jl"))
include(srcdir("solar_module.jl"))

CM = [ 0. 0. 1.
     0. 0. 2
     -1. -2. 0.]

ts = 1e-4


load_list = []
load = Dict()

#load["impedance"] = "RLC"
load["R"] = 14.0;
#load["L"] = 57.042;
#load["C"] = 39.18;
push!(load_list, load);


parameters["load"] = load_list;

ts = 1e-4
env = SimEnv(ts=ts, CM = CM, num_sources = 2, num_loads = 1, parameters = parameters, maxsteps = 500)
env.nc.parameters["load"][1]

