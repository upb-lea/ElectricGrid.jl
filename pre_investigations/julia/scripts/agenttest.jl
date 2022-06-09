using DrWatson
@quickactivate "MicroGridSimWithRL"

# using DifferentialEquations
# using Sundials
using Plots
# using LinearAlgebra
# using ControlSystems
# using BenchmarkTools
using ReinforcementLearning
using Flux
using StableRNGs
# using IntervalSets

include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent.jl"))

CM = [ 0.  1.
        -1.  0.]

parameters = Dict()
# LC filter
parameters["source"] = [Dict("fltr" => "LC", "R" => 0.4, "L1" => 2.3e-3, "C" => 10e-6)]
parameters["cable"] = [Dict("R" => 0.722, "L" => 0.955e-3, "C" => 8e-09)]
parameters["load"] = [Dict("impedance" => "R", "R" => 14)]

nc = NodeConstructor(num_source=1, num_loads=1, CM=CM, parameters=parameters)

#draw_graph(Grid_FC)   ---   not yet implemented

A, B, C, D = get_sys(nc)

limits = Dict("i_lim" => 20, "v_lim" => 600)

norm_array = vcat([limits[i] for j = 1:nc.num_source for i in ["i_lim", "v_lim"]], [limits["i_lim"] for i = 1:nc.num_connections] )
norm_array = vcat( norm_array, [limits["v_lim"] for i = 1:nc.num_loads] )

ns = length(A[1,:])
na = length(B[1,:])

# time step
ts = 1e-5

V_source = 300

global env = SimEnv(A=A, B=B, C=C, norm_array=norm_array, v_dc=V_source, ts=rationalize(ts))
global agent = create_agent(na, ns)

# ----------------------------------------------------------------------------------------
function execute_env(env::SimEnv, agent::Agent, t_len::Int, debug::Bool)
    if debug
        output = zeros(length(env.Ad[1,:]), t_len+1)
    else
        output = 0.0
    end

    RLBase.reset!(env)

    for i = 1:t_len
        action = agent(env)
        env(action)
        println(reward(env))
        if debug output[:,i+1] = env.state.*env.norm_array end
    end

    return output
end

P_required = 500 # W
V_required = 230 # V
PLoad = []
Pdiff = []

function reward_func(method::String, env::SimEnv)

    i_1, u_1, i_c1, u_l1 = env.state

    P_load = (env.norm_array[end] * u_l1)^2 / 14
    
    if method == "Power_exp"
        push!(PLoad, P_load)
        P_diff = -abs(P_required - P_load) 
        push!(Pdiff, P_diff)
        reward = exp(P_diff/130) - 1
    
    elseif method == "Power"
        reward = -abs(P_required - P_load) / (600 * 20)

    elseif method == "Voltage"
        reward = -((V_required - u_l1)/ 600) ^2
    end
    return reward
end

hook = TotalRewardPerEpisode()

No_Episodes = 10
run(
    agent,
    env,
    StopAfterEpisode(50),
    hook
)

plot(hook.rewards, title = "Total reward per episode")