using DrWatson
@quickactivate "ElectricGrid"

using ReinforcementLearning
using IntervalSets
using LinearAlgebra
using ControlSystems
using CUDA
using PlotlyJS


include(srcdir("node_constructor.jl"))
include(srcdir("electric_grid_env.jl"));
include(srcdir("sin_policy.jl"))
include(srcdir("data_hook.jl"))
include(srcdir("SolarModule.jl"))

function reward(env)
    #implement your reward function here
    return 1
end

CM = [ 0. 0. 1.
     0. 0. 2
     -1. -2. 0.]



parameters = Dict()
source_list = []
source1 = Dict()
source2 = Dict()

#source["pwr"] = 45000.0
#source["v_rip"] = 0.01556109320329396
#source["vdc"] = 750
#source["i_rip"] = 0.10108821490394984
source1["fltr"] = "LCL"
source1["source_type"] = "pv"
source1["R1"] = 0.4
source1["R_C"] = 0.0006
source1["L1"] = 2.3e-3
#source["R2"] = 0.4022094955070556   # needed for LCL
#source["L2"] = 0.001005523738767639
source1["C"] = 1e-6;

#source2["fltr"] = "L"
source2["source_type"] = "pv"
source2["R1"] = 0.4
source2["R_C"] = 0.0006
source2["L1"] = 2.3e-3
#source["R2"] = 0.4022094955070556   # needed for LCL
#source["L2"] = 0.001005523738767639
source2["C"] = 1e-6;

push!(source_list, source1)#, source2);

load_list = []
load = Dict()

load["impedance"] = "RLC"
#load["impedance"] = "R"
load["R"] = 14.0;
load["L"] = 57.042;
load["C"] = 39.18;
push!(load_list, load);

cable_list = []

cable = Dict()
cable["R"] = 0.722
cable["L"] = 0.264e-3
cable["C"] = 0.4e-6;
push!(cable_list, cable, cable);

parameters["source"] = source_list
parameters["cable"] = cable_list
parameters["load"] = load_list;
parameters["grid"] = Dict("fs" => 10000.0, "phase" => 1, "v_rms" => 230);

ts = 1e-4
env = ElectricGridEnv(reward_function = reward, ts=ts, use_gpu=false
, CM = CM, num_sources = 2, num_loads = 1, parameters = parameters, maxsteps = 500)


#######################################################################################
# Helper logging definitions till history works
N = 1000
mess = zeros(N)
u = zeros(N)


#######################################################################################
# load simple example policy which is able to interact with env in the RL-framwork
# policy has to enherit from AbstractPolicy.
# If learner wanted: replace policy by agent (agent has policy, which will be learned)
policy = sin_policy(action_space=action_space(env))


#######################################################################################
# Define data-logging hook
# define which states to store, to check what states are avalible type GetStateIds(env.nc) into command line

state_ids = GetStateIds(env.nc)
action_ids = GetActionIds(env.nc)

state_ids_agent = filter(x -> split(x, "_")[1] == "source1", state_ids)
action_ids_agent = filter(x -> split(x, "_")[1] == "source1", action_ids)

#plt_state_ids = ["u_f1_a", "u_f1_b", "u_f1_c", "u_f2_a", "u_f2_b", "u_f2_c", "i_f1_a", "i_f1_b", "i_f1_c", "i_f2_a", "i_f2_b", "i_f2_c"]  
plt_state_ids = ["source1_i_L1", "source2_i_L1"]#, "i_2_a", "i_2_b", "i_2_c"]  
# define which states to store, to check what states are avalible type GetActionIds(env.nc) into command line 
plt_action_ids = ["source1_u", "source2_u"]#, "u_v2_a", "u_v2_b", "u_v2_c"]
hook = DataHook(collect_state_ids = plt_state_ids, collect_action_ids = plt_action_ids, collect_vdc_idx = [1, 2])

#######################################################################################
# GOAL: Use run function provided by ReinforcementLearning.jl to be able to interact 
#       with our env in the RL-interface kind of manner to be able to use standard 
#       RL-Algorithms
# 
run(policy, env, StopAfterEpisode(1), hook)

#TODO
# this will be shifted to render.jl soon
RenderHookResults(hook=hook, vdc_to_plot=[1, 2])
RenderHookResults(; hook = hook, states_to_plot = ["source1_i_L1"] )