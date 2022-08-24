using DrWatson
@quickactivate "dare"

using ReinforcementLearning
using PlotlyJS

include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent_ddpg.jl"))
include(srcdir("data_hook.jl"))
include(srcdir("plotting.jl"))


function reference(t)
    230 * sin(2*pi*50 * t)
end

function reward(env)
    u_l1_index = findfirst(x -> x == "u_load1", env.state_ids)
    u_l1 = env.state[u_l1_index]

    return -(abs(reference(env.t) - (env.norm_array[u_l1_index] * u_l1))/300)
end

function featurize(x0 = nothing, t0 = nothing; env = nothing) 
    if isnothing(env)
        state = copy(x0)
        push!(state, reference(t0))
    else
        state = copy(env.x)
        push!(state, reference(env.t))
    end
    return state
end

env_cuda = false
agent_cuda = false

CM = [ 0. 1.
      -1. 0.]

parameters = Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}("L1"=>0.0023, "R_C"=>0.4, "C"=>1.0e-5, "R1"=>0.4, "fltr"=>"LC"),
                    #Dict{Any, Any}("v_rip"=>0.015, "L1"=>0.0023, "vdc"=>700, "R1"=>0.4, "i_rip"=>0.12, "pwr"=>10000.0, "fltr"=>"L")
                    ],
    "load"   => Any[
                    Dict{Any, Any}("R"=>14, "impedance"=>"R")
                    #Dict{Any, Any}("C"=>0.381, "L"=>0.2, "R"=>14, "impedance"=>"RLC")
                    ],
    "cable"  => Any[
                    Dict{Any, Any}("C"=>4.0e-7, "L"=>0.000264, "R"=>0.722),
                    #Dict{Any, Any}("C"=>4.0e-7, "L"=>0.000264, "R"=>0.722)
                    ],
    "grid"   => Dict{Any, Any}("fs"=>10000.0, "phase"=>1, "v_rms"=>230)
)


# time step
ts = 1e-4

V_source = 300

env = SimEnv(reward_function = reward, featurize = featurize, v_dc=V_source, ts=ts, use_gpu=env_cuda, CM = CM, num_sources = 1, num_loads = 1, parameters = parameters)

ns = length(env.state_space)
na = length(env.action_space)
agent = create_agent_ddpg(na = na, ns = ns, use_gpu = agent_cuda)

hook = DataHook(save_best_NNA = true, plot_rewards = true)

run(agent, env, StopAfterEpisode(500), hook)





# PLOT rewards in 3D plot over every episode

# plot_rewards_3d(hook)



# PLOT a test run with the best behavior_actor NNA so far

plot_best_results(;agent = agent, env = env, hook = hook, state_ids_to_plot = ["u_f1", "u_1"], plot_reward = true, plot_reference = true)#, "u_load1"])

