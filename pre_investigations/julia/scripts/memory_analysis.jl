using DrWatson
@quickactivate "MicroGridSimWithRL"

using PProf, Profile
# using DifferentialEquations
# using Sundials
using Plots
# using LinearAlgebra
# using ControlSystems
using BenchmarkTools
using ReinforcementLearning
using Flux
using StableRNGs
using IntervalSets
using TimerOutputs
using JSON

# include(srcdir("collect_timing_results.jl"))
include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent.jl"))
include(srcdir("run_timed.jl"))

function reward(env)
    #implement your reward function here

    P_required = 466 # W
    V_required = 230 # V

    u_l1_index = findfirst(x -> x == "u_l1", env.state_ids)

    u_l1 = env.state[u_l1_index]

    P_load = (env.norm_array[u_l1_index] * u_l1)^2 / 14
    
    # P_diff = -abs(P_required - P_load) 
    # reward = exp(P_diff/130) - 1

    reward = -(abs(V_required - (env.norm_array[u_l1_index] * u_l1))/300)

    return reward
end

env_cuda = false
agent_cuda = false

CM = [ 0. 0. 1.
        0. 0. 2
        -1. -2. 0.]

parameters = Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}("L1"=>0.0023, "R_C"=>0.4, "C"=>1.0e-5, "R1"=>0.4, "fltr"=>"LC"),
                    Dict{Any, Any}("v_rip"=>0.015, "L1"=>0.0023, "vdc"=>700, "R1"=>0.4, "i_rip"=>0.12, "pwr"=>10000.0, "fltr"=>"L")
                    ],
    "load"   => Any[
                    Dict{Any, Any}("R"=>14, "impedance"=>"R")
                    #Dict{Any, Any}("C"=>0.381, "L"=>0.2, "R"=>14, "impedance"=>"RLC")
                    ],
    "cable"  => Any[
                    Dict{Any, Any}("C"=>4.0e-7, "L"=>0.000264, "R"=>0.722),
                    Dict{Any, Any}("C"=>4.0e-7, "L"=>0.000264, "R"=>0.722)
                    ],
    "grid"   => Dict{Any, Any}("fs"=>10000.0, "phase"=>1, "v_rms"=>230)
)

nc = NodeConstructor(num_sources = 2, num_loads = 1, CM = CM, parameters = parameters)

A, B, C, D = get_sys(nc)

limits = Dict("i_lim" => 20, "v_lim" => 600)

states = get_state_ids(nc)
norm_array = []
for state_name in states
    if startswith(state_name, "i")
        push!(norm_array, limits["i_lim"])
    elseif startswith(state_name, "u")
        push!(norm_array, limits["v_lim"])
    end
end

ns = length(A[1,:])
na = length(B[1,:])

# time step
ts = 1e-5

V_source = 300

x0 = [ 0.0 for i = 1:length(A[1,:]) ]

if env_cuda
    A = CuArray(A)
    B = CuArray(B)
    C = CuArray(C)
    x0 = CuArray(x0)
end

env = SimEnv(A=A, B=B, C=C, norm_array=norm_array, state_ids = states, rewardfunction = reward, x0=x0, v_dc=V_source, ts=ts, convert_state_to_cpu=true, maxsteps=600)
agent = create_agent_ddpg(na = na, ns = ns, use_gpu = agent_cuda)

hook = DataHook(save_best_NNA = true, plot_rewards = true)