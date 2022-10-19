using DrWatson
@quickactivate "dare"

using ReinforcementLearning
using IntervalSets
using LinearAlgebra
using ControlSystems
using CUDA
using PlotlyJS


include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"));
include(srcdir("sin_policy.jl"))
include(srcdir("data_hook.jl"))
include(srcdir("pv_module.jl"))

CM = [ 0. 0. 1.
     0. 0. 2
     -1. -2. 0.]

ts = 1e-4


source1 = Dict()
source_list = []
source1["fltr"] = "LC"
source1["source_type"] = "pv"
source1["R1"] = 0.4
source1["R_C"] = 0.0006
source1["L1"] = 2.3e-3
#source["R2"] = 0.4022094955070556   # needed for LCL
#source["L2"] = 0.001005523738767639
#source1["C"] = 1e-6;

parameters = Dict()
parameters["source"] = source_list
push!(source_list, source1)#, source2);

env = SimEnv(ts=ts, num_sources = 3, num_loads = 2, maxsteps = 500, parameters=parameters);

