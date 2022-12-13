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


parameters = Dict()
source_list = []
source1 = Dict()
source2 = Dict()

source1["pwr"] = 5000.0
source1["vdc"] = 750
source1["fltr"] = "LC"
source1["R1"] = 0.4
source1["R_C"] = 0.0006
source1["L1"] = 2.3e-3
#source["R2"] = 0.4022094955070556   # needed for LCL
#source["L2"] = 0.001005523738767639
source1["C"] = 1e-6;

source2["pwr"] = 5000.0
source2["vdc"] = 750
source2["fltr"] = "LCL"
source2["R1"] = 0.4
source2["R_C"] = 0.0006
source2["L1"] = 2.3e-3
source2["R2"] = 0.4
source2["L2"] = 2.3e-3
source2["C"] = 1e-6;

cable = Dict()
cable_list = []
cable["R"] = 0.722
cable["L"] = 0.264e-3
cable["C"] = 0.4e-6;
push!(cable_list, cable, cable);
parameters["cable"] = cable_list



push!(source_list, source2, source1);

parameters["source"] = source_list


env = SimEnv(ts=ts, CM = CM, num_sources = 2, num_loads = 1, parameters = parameters, maxsteps = 500)
env.nc.parameters

