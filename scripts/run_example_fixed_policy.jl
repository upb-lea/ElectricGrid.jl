using DrWatson
@quickactivate "dare"

using ReinforcementLearning
using IntervalSets
using LinearAlgebra
using ControlSystems
using CUDA
#using Plots
using PlotlyJS

include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"));
include(srcdir("data_hook.jl"))
include(srcdir("sin_policy.jl"))

function reward(env)
    #implement your reward function here
    return 1
end

#_______________________________________________________________________________
# Parameters - Time simulation
Timestep = 75 #time step in μs ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_final = 0.04 #time in seconds, total simulation run time

#_______________________________________________________________________________
# Environment Calcs

ts = Timestep*1e-6
t = 0:ts:t_final # time

fs = 1/ts

N = length(t)

CM = [0. 1.
   -1. 0.]

#= CM = [ 0. 0. 1.
        0. 0. 2
        -1. -2. 0.] =#

parameters = Dict()
source_list = []
source = Dict()

source["fltr"] = "L"
source["R1"] = 0.4
source["L1"] = 0.00034

push!(source_list, source);

load_list = []
load = Dict()

load["impedance"] = "RL"
load["R"] = 5.289
load["L"] = 0.0126
push!(load_list, load);

cable_list = []

cable = Dict()
cable["R"] = 0.722
cable["L"] = 0.000024
cable["C"] = 0.4e-9;
push!(cable_list, cable);

parameters["source"] = source_list
parameters["cable"] = cable_list
parameters["load"] = load_list;
parameters["grid"] = Dict("fs" => fs, "phase" => 3, "v_rms" => 230);

env = SimEnv(reward_function = reward,  v_dc = 1, ts = ts, use_gpu = false
, CM = CM, num_sources = 1, num_loads = 1, parameters = parameters, maxsteps = N-1)

#######################################################################################
# GOAL: Use run function provided by ReinforcementLearning.jl to be able to interact 
#       with our env in the RL-interface kind of manner to be able to use standard 
#       RL-Algorithms
# 
#%% Starting time simulation

reset!(env)

input_action_a = zeros(N-1)
input_action_b = zeros(N-1)
input_action_c = zeros(N-1)

env_action_a = zeros(N-1)
env_action_b = zeros(N-1)
env_action_c = zeros(N-1)

vout_a = zeros(N-1)
vout_b = zeros(N-1)
vout_c = zeros(N-1)

iout_a = zeros(N-1)
iout_b = zeros(N-1)
iout_c = zeros(N-1)

ns = Int(size(env.sys_d.A, 2)/3) # system size
unbalance_ab = maximum((env.sys_d.A)[1:5, 1:5] .- (env.sys_d.A)[6:10, 6:10]) 
unbalance_bc = maximum((env.sys_d.A)[6:10, 6:10] .- (env.sys_d.A)[11:15, 11:15]) 
unbalance_ac = maximum((env.sys_d.A)[1:5, 1:5] .- (env.sys_d.A)[11:15, 11:15]) 

#get_state_ids(env.nc)
V_poc_loc = [2; 7; 12]
#state_index = findfirst(x -> x == "u_1_a", env.state_ids)
I_poc_loc = [1; 6; 11]
#state_index = findfirst(x -> x == "i_1_a", env.state_ids)

#######################################################################################
plt_state_ids = ["i_1_a", "i_1_b", "i_1_c"] 
#plt_state_ids = ["u_1_a", "u_1_b", "u_1_c", "u_2_a", "u_2_b", "u_2_c"] 
#plt_action_ids = ["u_v1_a", "u_v1_b", "u_v1_c"]
#plt_action_ids = ["u_v1_a", "u_v1_b", "u_v1_c", "u_v2_a", "u_v2_b", "u_v2_c"]
hook = DataHook(collect_state_ids = plt_state_ids#= , collect_action_ids = plt_action_ids =#)

policy = sin_policy(action_space = action_space(env), ts = ts)
run(policy, env, StopAfterEpisode(1), hook)

#= for i in 1:N-1

    #u = [1*sin.(50*2*pi*t[i]) for j = 1:3]
    #action = vcat(u,u)

    action = policy(env)
    env(action)

    s = 1 # select source
    num_sources = 1 # total sources
    input_action_a[i] = action[s + num_sources*(1 - 1)]
    input_action_b[i] = action[s + num_sources*(2 - 1)]
    input_action_c[i] = action[s + num_sources*(3 - 1)]

    env_action_a[i] = env.action[s + num_sources*(1 - 1)]
    env_action_b[i] = env.action[s + num_sources*(2 - 1)]
    env_action_c[i] = env.action[s + num_sources*(3 - 1)]

    vout_a[i] = env.x[V_poc_loc[1, s]]
    vout_b[i] = env.x[V_poc_loc[2, s]]
    vout_c[i] = env.x[V_poc_loc[3, s]]

    iout_a[i] = env.x[I_poc_loc[1, s]]
    iout_b[i] = env.x[I_poc_loc[2, s]]
    iout_c[i] = env.x[I_poc_loc[3, s]]

end =#

plot_hook_results(hook = hook)

#= T_plot_start = 0
T_plot_end = 10
fsys = 50

if T_plot_end > t_final*fsys
    T_plot_end = t_final*fsys
end

Nps = 1/ts
N_plot_start = convert(Int64, round((T_plot_start/fsys  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fsys  - 1/Nps)*Nps))
range = N_plot_start:N_plot_end

v_out = Plots.plot(t[range], vout_a[range], label = "a",
            xlabel = "time", ylabel = "V poc", title = "V POC Env.x")
v_out = Plots.plot!(t[range], vout_b[range], label = "b")
v_out = Plots.plot!(t[range], vout_c[range], label = "c")
display(v_out)
 
i_out = Plots.plot(t[range], iout_a[range], label = "a",
            xlabel = "time", ylabel = "I poc", title = "I POC Env.x")
i_out = Plots.plot!(t[range], iout_b[range], label = "b")
i_out = Plots.plot!(t[range], iout_c[range], label = "c")
display(i_out)

u = Plots.plot(t[range], input_action_a[range], label = "a", 
        xlabel = "time", ylabel = "V inv", title = "Policy Action")
u = Plots.plot!(t[range], input_action_b[range], label = "b")
u = Plots.plot!(t[range], input_action_c[range], label = "c")
display(u)

u = Plots.plot(t[range], env_action_a[range], label = "a",
            xlabel = "time", ylabel = "V inv", title = "Env.Action")
u = Plots.plot!(t[range], env_action_b[range], label = "b")
u = Plots.plot!(t[range], env_action_c[range], label = "c")
display(u) =#

