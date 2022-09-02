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

function reward(env)
    #implement your reward function here
    return 1
end

#_______________________________________________________________________________
# Parameters - Time simulation
Timestep = 10 #time step in μs ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_final = 0.4 #time in seconds, total simulation run time

#_______________________________________________________________________________
# Environment Calcs

ts = Timestep*1e-6
t = 0:ts:t_final # time

N = length(t)

CM = [ 0. 0. 1.
        0. 0. 2
        -1. -2. 0.]

parameters = Dict()
source_list = []
source = Dict()

source["fltr"] = "LC"
source["R1"] = 0.4
source["R_C"] = 0.0006
source["L1"] = 2.3e-3
source["C"] = 1e-6;

push!(source_list, source, source);

load_list = []
load = Dict()

load["impedance"] = "R"
load["R"] = 14.0;
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

env = SimEnv(reward_function = reward,  v_dc = 1000, ts = ts, use_gpu = false
, CM = CM, num_sources = 2, num_loads = 1, parameters = parameters, maxsteps = 100)

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

V_poc_loc = [3 6; 12 15; 21 24]
#state_index = findfirst(x -> x == "u_1_a", env.state_ids)

#######################################################################################
plt_state_ids = ["u_f1_a", "u_f1_b", "u_f1_c", "u_f2_a", "u_f2_b", "u_f2_c"] 
plt_action_ids = ["u_v1_a", "u_v1_b", "u_v1_c", "u_v2_a", "u_v2_b", "u_v2_c"]
hook = DataHook(collect_state_ids = plt_state_ids, collect_action_ids = plt_action_ids)

for i in 1:N-1

    u = [1*sin.(50*2*pi*t[i]) for j = 1:3]
    action = vcat(u,u)

    env(action)

    s = 1 # select source
    num_sources = 2 # total sources
    input_action_a[i] = action[s + num_sources*(1 - 1)]
    input_action_b[i] = action[s + num_sources*(2 - 1)]
    input_action_c[i] = action[s + num_sources*(3 - 1)]

    env_action_a[i] = env.action[s + num_sources*(1 - 1)]
    env_action_b[i] = env.action[s + num_sources*(2 - 1)]
    env_action_c[i] = env.action[s + num_sources*(3 - 1)]

    vout_a[i] = env.state[V_poc_loc[1, s]]
    vout_b[i] = env.state[V_poc_loc[2, s]]
    vout_c[i] = env.state[V_poc_loc[3, s]]

end

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

v_out = plot(t[range], vout_a[range], label = "a",
            xlabel = "time", ylabel = "V poc", title = "Env.State")
v_out = plot!(t[range], vout_b[range], label = "b")
v_out = plot!(t[range], vout_c[range], label = "c")
display(v_out)
 
u = plot(t[range], input_action_a[range], label = "a", 
        xlabel = "time", ylabel = "V inv", title = "Policy Action")
u = plot!(t[range], input_action_b[range], label = "b")
u = plot!(t[range], input_action_c[range], label = "c")
display(u)

u = plot(t[range], env_action_a[range], label = "a",
            xlabel = "time", ylabel = "V inv", title = "Env.Action")
u = plot!(t[range], env_action_b[range], label = "b")
u = plot!(t[range], env_action_c[range], label = "c")
display(u) =#

