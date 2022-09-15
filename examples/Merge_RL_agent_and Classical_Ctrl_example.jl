using DrWatson
@quickactivate "dare"

using ReinforcementLearning
using PlotlyJS

include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent_ddpg.jl"))
include(srcdir("data_hook.jl"))
include(srcdir("Classical_Control.jl"))
include(srcdir("Power_System_Theory.jl"))


function reference(t)
    
    u = [230 * sin.(2*pi*50*t .- 2/3*pi*(i-1)) for i = 1:3]
    #return vcat(u,u)  # to control 2 sources
    return u
end

function reward(env)
    
    u_l1_index = findfirst(x -> x == "u_1_a", env.state_ids)
    u_l2_index = findfirst(x -> x == "u_1_b", env.state_ids)
    u_l3_index = findfirst(x -> x == "u_1_c", env.state_ids)
    
    u_l1 = env.state[u_l1_index]
    u_l2 = env.state[u_l2_index]
    u_l3 = env.state[u_l3_index]

    u = [u_l1, u_l2, u_l3]
    refs = reference(env.t)

    r = -(sum(abs.(refs/600 - u)/3))
    return r
end

function featurize(x0 = nothing, t0 = nothing; env = nothing) 
    if isnothing(env)
        state = copy(x0)
        #push!(state, t0/(1e-4*500))  
        state = vcat(state, reference(t0)/600)
    else
        state = copy(env.state)
        #push!(state, env.t/(env.ts*env.maxsteps))
        #push!(state, reference(env.t)/600)
        state = vcat(state, reference(env.t)/600)
     
    end
    return state
end

env_cuda = false
agent_cuda = false

CM = [ 0. 1.
     -1.  0.]

#CM = [ 0. 0. 1.
#     0. 0. 2
#     -1. -2. 0.]

parameters = Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}("L1"=>0.0023, "R_C"=>0.4, "C"=>1.0e-5, "R1"=>0.4, "fltr"=>"LC"),
                    #Dict{Any, Any}("L1"=>0.0023, "R_C"=>0.4, "C"=>1.0e-5, "R1"=>0.4, "fltr"=>"LC") 
                    ],
    "load"   => Any[
                    #Dict{Any, Any}("R"=>14, "impedance"=>"R")
                    Dict{Any, Any}("C"=>0.381, "L"=>0.2, "R"=>14, "impedance"=>"RLC")
                    ],
    "cable"  => Any[
                    Dict{Any, Any}("C"=>4.0e-7, "L"=>0.000264, "R"=>0.722),
                    #Dict{Any, Any}("C"=>4.0e-7, "L"=>0.000264, "R"=>0.722)
                    ],
    "grid"   => Dict{Any, Any}("fs"=>10000.0, "phase"=>3, "v_rms"=>230)
)


# time step
ts = 1e-4

V_source = 300

env = SimEnv(reward_function = reward, featurize = featurize, 
v_dc=V_source, ts=ts, use_gpu=env_cuda, CM = CM, num_sources = 1, num_loads = 1, parameters = parameters,
maxsteps=1000)

ns = length(env.state_space)
na = length(env.action_space)
agent = create_agent_ddpg(na = na, ns = ns, use_gpu = agent_cuda)


# toDo: action_space half
#Animo = Classical_Policy(action_space = action_space(env), t_final = t_final, 
#fs = fs, num_sources = num_sources)

#Modes = [5]

# tune controller
#Source_Initialiser(env, Animo, Modes)

plt_state_ids = ["u_f1_a", "u_f1_b", "u_f1_c"]#, "u_f2_a", "u_f2_b", "u_f2_c"]
plt_action_ids = ["u_v1_a", "u_v1_b", "u_v1_c"]
hook = DataHook(collect_state_ids = plt_state_ids, collect_action_ids = plt_action_ids, save_best_NNA = true, collect_reference = true)

run(agent, env, StopAfterEpisode(50), hook)


###############################
# Plotting
plot_hook_results(; hook=hook, states_to_plot = ["u_f1_a", "u_f1_b", "u_f1_c"], actions_to_plot = ["u_v1_a", "u_v1_b", "u_v1_c"] ,plot_reward = false, plot_reference = true, episode = 50)

plot_hook_results(hook=hook, episode=38, plot_reference = true)

#plot_best_results(;agent = agent, env = env, hook = hook, plot_reward = true, plot_reference = true)#, "u_load1"])
