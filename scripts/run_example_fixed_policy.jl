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


function reward(env)
    #implement your reward function here
    return 1
end

CM = [ 0. 0. 1.
        0. 0. 2
        -1. -2. 0.]

#CM = [0. 1.
#    -1. 0.]


parameters = Dict()
source_list = []
source = Dict()

#source["pwr"] = 45000.0
#source["v_rip"] = 0.01556109320329396
#source["vdc"] = 750
#source["i_rip"] = 0.10108821490394984
source["fltr"] = "LC"
source["R1"] = 0.4
source["R_C"] = 0.0006
source["L1"] = 2.3e-3
#source["R2"] = 0.4022094955070556   # needed for LCL
#source["L2"] = 0.001005523738767639
source["C"] = 1e-6;

push!(source_list, source, source);

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
parameters["grid"] = Dict("fs" => 10000.0, "phase" => 3, "v_rms" => 230);

ts = 1e-4
env = SimEnv(reward_function = reward,  v_dc=300, ts=ts, use_gpu=false
, CM = CM, num_sources = 2, num_loads = 1, parameters = parameters, maxsteps = 100)


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
# define which states to store, to check what states are avalible type get_state_ids(env.nc) into command line
plt_state_ids = ["u_f1_a", "u_f1_b", "u_f1_c", "u_f2_a", "u_f2_b", "u_f2_c"]  
# define which states to store, to check what states are avalible type get_action_ids(env.nc) into command line 
plt_action_ids = ["u_v1_a", "u_v1_b", "u_v1_c", "u_v2_a", "u_v2_b", "u_v2_c"]
hook = DataHook(collect_state_ids = plt_state_ids, collect_action_ids = plt_action_ids)

#######################################################################################
# GOAL: Use run function provided by ReinforcementLearning.jl to be able to interact 
#       with our env in the RL-interface kind of manner to be able to use standard 
#       RL-Algorithms
# 
run(policy, env, StopAfterEpisode(1), hook)

#TODO
# this will be shifted to plotting.jl soon
plot_hook_results(hook=hook)

