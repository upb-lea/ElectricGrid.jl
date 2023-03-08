using DrWatson
@quickactivate "JEG"

using ReinforcementLearning
using IntervalSets
using LinearAlgebra
using ControlSystems
using CUDA
using Plots

include(srcdir("node_constructor.jl"))
include(srcdir("electric_grid_env.jl"));
include(srcdir("sin_policy.jl"))
include(srcdir("data_hook.jl"))


function reward(env)
    #implement your reward function here
    return 1
end

CM = [ 0. 0. 1.
        0. 0. 2
        -1. -2. 0.]


parameters = Dict()
source_list = []
source = Dict()

source["pwr"] = 45000.0 #VA not W
source["v_rip"] = 0.01556109320329396
source["vdc"] = 750
source["i_rip"] = 0.10108821490394984
source["fltr"] = "LC"
source["R1"] = 0.4022094955070556
source["R_C"] = 0.0006447094780419011
source["L1"] = 0.001005523738767639
#source["R2"] = 0.4022094955070556   # needed for LCL
#source["L2"] = 0.001005523738767639
source["C"] = 2.302533850149647e-5;

push!(source_list, source, source);

load_list = []
load = Dict()

load["impedance"] = "RLC"
load["R"] = 30236.0;
load["L"] = 57.042;
load["C"] = 39.18;
push!(load_list, load);

cable_list = []

cable = Dict()
cable["R"] = 6.84059
cable["L"] = 0.00250127
cable["C"] = 3.7898e-6;
push!(cable_list, cable, cable);

parameters["source"] = source_list
parameters["cable"] = cable_list
parameters["load"] = load_list;
parameters["grid"] = Dict("fs" => 10000.0, "phase" => 3, "v_rms" => 230);

#######################################################################################
# Define grid using random initialization
power_grid = NodeConstructor(num_sources=2, num_loads=1, S2S_p=1, S2L_p=1, CM = CM, parameters = parameters);
A, B, C, D = GetSystem(power_grid)
ns = length(A[1,:]) # get num of states
ni = length(B[1,:]) # get num of inputs
x0 = [0.0 for i = 1:ns]
env = ElectricGridEnv(A=A, B=B, C=C, D=D, x0=x0, state_ids=GetStateIds(power_grid), rewardfunction = reward)


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
hook = DataHook(state_ids = ["u_f1", "u_1", "u_l1", "i_f2"])

#######################################################################################
# GOAL: Use run function provided by ReinforcementLearning.jl to be able to interact 
#       with our env in the RL-interface kind of manner to be able to use standard 
#       RL-Algorithms
# 
run(policy, env, StopAfterStep(N), hook)

#= # For running without run()-command
reset!(env)
for i in 1:N
    action = policy(env)  
    env(action)
    
    #######################
    # Logg one example actiona and an exemplary state
    # toDo: shfit logging and plotting into "History" using hooks
    mess[i] = env.state[1]
    u[i] = action[1]
end
display(plot(mess, xlabel="time", ylabel="state_no_1"))
display(plot(u, xlabel="time", ylabel="action_no_1"))
 =#