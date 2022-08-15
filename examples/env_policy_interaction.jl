using DrWatson
@quickactivate "dare"

using ReinforcementLearning
using IntervalSets
using LinearAlgebra
using ControlSystems
using CUDA
using Plots

include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"));
include(srcdir("sin_policy.jl"))

#######################################################################################
# Define grid using random initialization
power_grid = NodeConstructor(num_sources=1, num_loads=1, S2S_p=1, S2L_p=1);
A, B, C, D = get_sys(power_grid)
ns = length(A[1,:]) # get num of states
ni = length(B[1,:]) # get num of inputs
x0 = [0.0 for i = 1:ns]
env = SimEnv(A=A, B=B, C=C, D=D, x0=x0)


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
# GOAL: Use run function provided by ReinforcementLearning.jl to be able to interact 
#       with our env in the RL-interface kind of manner to be able to use standard 
#       RL-Algorithms
# 
run(policy, env, StopAfterStep(10))

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



