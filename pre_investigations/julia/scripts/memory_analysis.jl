using DrWatson
quickactivate("../../.", "dare")

using BenchmarkTools
using ReinforcementLearning
using TimerOutputs
using DataFrames
using JSON
using PlotlyJS

include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent_ddpg.jl"))
include(srcdir("data_hook.jl"))
include(srcdir("plotting.jl"))

function reward(env)
    #implement your reward function here

    P_required = 466 # W
    V_required = 230 # V

    u_l1_index = findfirst(x -> x == "u_load1_c", env.state_ids)
    print(u_l1_index) 
    u_l1 = env.state[u_l1_index]

    P_load = (env.norm_array[u_l1_index] * u_l1)^2 / 14
    
    # P_diff = -abs(P_required - P_load) 
    # reward = exp(P_diff/130) - 1

    reward = -(abs(V_required - (env.norm_array[u_l1_index] * u_l1))/300)

    return reward
end


function collectresults(timer::TimerOutput, node::Int, env_cuda::Bool)

    if env_cuda
        global sections = ["hook_pre_experiment", "agent_pre_experiment", "agent_pre_episode",
                "hook_pre_episode", "get_action", "agent_pre_action", "hook_pre_action", "action_cuarray_conversion",
                "env_step", "agent_post_action", "hook_post_action", "agent_post_episode",
                "hook_post_episode", "hook_post_experiment"]
    else
        global sections = ["hook_pre_experiment", "agent_pre_experiment", "agent_pre_episode",
                "hook_pre_episode", "get_action", "agent_pre_action", "hook_pre_action",
                "env_step", "agent_post_action", "hook_post_action", "agent_post_episode",
                "hook_post_episode", "hook_post_experiment"]
    end

    timerow = []
    allocatedrow = []
    for section in sections
        push!(timerow, TimerOutputs.time(to[section]))
        push!(allocatedrow, TimerOutputs.allocated(to[section]))
    end
    
    global results_time[!, string(node)] = timerow

    global results_allocated[!, string(node)] = allocatedrow
end


function memory_analysis(to, env_cuda, agent_cuda, num_nodes)
    CM_list = JSON.parsefile("pre_investigations/julia/src/CM_matrices/CM_nodes" * string(num_nodes) * ".json")

    CM = reduce(hcat, CM_list[1])'
    CM = convert(Matrix{Int}, CM)

    # nc = NodeConstructor(num_sources = num_nodes, num_loads = num_nodes, CM = CM)

    # A, B, C, D = get_sys(nc)

    # limits = Dict("i_lim" => 20, "v_lim" => 600)

    # states = get_state_ids(nc)
    # norm_array = []
    # for state_name in states
    # if startswith(state_name, "i")
    #     push!(norm_array, limits["i_lim"])
    # elseif startswith(state_name, "u")
    #     push!(norm_array, limits["v_lim"])
    # end
    # end

    # ns = length(A[1,:])
    # na = length(B[1,:])

    # # time step
    # ts = 1e-5

    # V_source = 300

    # x0 = [ 0.0 for i = 1:length(A[1,:]) ]

    # if env_cuda
    # A = CuArray(A)
    # B = CuArray(B)
    # C = CuArray(C)
    # x0 = CuArray(x0)
    # end

    # env = SimEnv(A=A, B=B, C=C, norm_array=norm_array, state_ids = states, rewardfunction = reward, x0=x0, v_dc=V_source, ts=ts, convert_state_to_cpu=true, maxsteps=600)
    
    env = SimEnv(num_sources = num_nodes, num_loads = num_nodes, CM = CM, 
                # parameters = parameters, 
                reward_function = reward, 
                maxsteps=600, 
                use_gpu=env_cuda)

    ns = length(env.sys_d.A[1,:])
    na = length(env.sys_d.B[1,:])

    agent = create_agent_ddpg(na = na, ns = ns, use_gpu = agent_cuda)

    hook = DataHook(save_best_NNA = true)

    println("Hook - Pre Experiment Stage")
    hookbefore = deepcopy(hook)
    hook(PRE_EXPERIMENT_STAGE, agent, env);
    hook = deepcopy(hookbefore)
    @time hook(PRE_EXPERIMENT_STAGE, agent, env);
    hook = deepcopy(hookbefore)
    @timeit to "hook_pre_experiment" hook(PRE_EXPERIMENT_STAGE, agent, env);
    println(""); println(""); println("");


    println("Agent - Pre Experiment Stage")
    agentbefore = deepcopy(agent)
    agent(PRE_EXPERIMENT_STAGE, env);
    agent = deepcopy(agentbefore)
    @time agent(PRE_EXPERIMENT_STAGE, env);
    agent = deepcopy(agentbefore)
    @timeit to "agent_pre_experiment" agent(PRE_EXPERIMENT_STAGE, env);
    println(""); println(""); println("");


    println("Agent - Pre Episode Stage")
    agentbefore = deepcopy(agent)
    agent(PRE_EPISODE_STAGE, env);
    agent = deepcopy(agentbefore)
    @time agent(PRE_EPISODE_STAGE, env);
    agent = deepcopy(agentbefore)
    @timeit to "agent_pre_episode" agent(PRE_EPISODE_STAGE, env);
    println(""); println(""); println("");


    println("Hook - Pre Episode Stage")
    hookbefore = deepcopy(hook)
    hook(PRE_EPISODE_STAGE, agent, env);
    hook = deepcopy(hookbefore)
    @time hook(PRE_EPISODE_STAGE, agent, env);
    hook = deepcopy(hookbefore)
    @timeit to "hook_pre_episode" hook(PRE_EPISODE_STAGE, agent, env);
    println(""); println(""); println("");


    println("Agent - Get Action")
    agentbefore = deepcopy(agent)
    action = agent(env);
    agent = deepcopy(agentbefore)
    @time action = agent(env);
    agent = deepcopy(agentbefore)
    @timeit to "get_action" action = agent(env);
    println(""); println(""); println("");


    println("Agent - Pre Action Stage")
    agentbefore = deepcopy(agent)
    agent(PRE_ACT_STAGE, env, action);
    agent = deepcopy(agentbefore)
    @time agent(PRE_ACT_STAGE, env, action);
    agent = deepcopy(agentbefore)
    @timeit to "agent_pre_action" agent(PRE_ACT_STAGE, env, action);
    println(""); println(""); println("");


    println("Hook - Pre Action Stage")
    hookbefore = deepcopy(hook)
    hook(PRE_ACT_STAGE, agent, env, action);
    hook = deepcopy(hookbefore)
    @time hook(PRE_ACT_STAGE, agent, env, action);
    hook = deepcopy(hookbefore)
    @timeit to "hook_pre_action" hook(PRE_ACT_STAGE, agent, env, action);
    println(""); println(""); println("");


    if env_cuda
        println("Action CuArray Conversion")
        if action isa Array
            action = CuArray(action)
            action = Array(action)
            @time action = CuArray(action)
            action = Array(action)
            @timeit to "action_cuarray_conversion" action = CuArray(action)
        else
            action = CuArray([action])
            action = Float64(Array(action))
            @time action = CuArray([action])
            action = Float64(Array(action))
            @timeit to "action_cuarray_conversion" action = CuArray([action])
        end
        println(""); println(""); println("");
    end


    println("Env Calculation")
    envbefore = deepcopy(env)
    env(action);
    env = deepcopy(envbefore)
    @time env(action);
    env = deepcopy(envbefore)
    @timeit to "env_step" env(action);
    println(""); println(""); println("");


    println("Agent - Post Action Stage")
    agentbefore = deepcopy(agent)
    agent(POST_ACT_STAGE, env);
    agent = deepcopy(agentbefore)
    @time agent(POST_ACT_STAGE, env);
    agent = deepcopy(agentbefore)
    @timeit to "agent_post_action" agent(POST_ACT_STAGE, env);
    println(""); println(""); println("");


    println("Hook - Post Action Stage")
    hookbefore = deepcopy(hook)
    hook(POST_ACT_STAGE, agent, env);
    hook = deepcopy(hookbefore)
    @time hook(POST_ACT_STAGE, agent, env);
    hook = deepcopy(hookbefore)
    @timeit to "hook_post_action" hook(POST_ACT_STAGE, agent, env);
    println(""); println(""); println("");


    println("Agent - Post Episode Stage")
    agentbefore = deepcopy(agent)
    agent(POST_EPISODE_STAGE, env);
    agent = deepcopy(agentbefore)
    @time agent(POST_EPISODE_STAGE, env);
    agent = deepcopy(agentbefore)
    @timeit to "agent_post_episode" agent(POST_EPISODE_STAGE, env);
    println(""); println(""); println("");


    println("Hook - Post Episode Stage")
    hookbefore = deepcopy(hook)
    hook(POST_EPISODE_STAGE, agent, env);
    hook = deepcopy(hookbefore)
    @time hook(POST_EPISODE_STAGE, agent, env);
    hook = deepcopy(hookbefore)
    @timeit to "hook_post_episode" hook(POST_EPISODE_STAGE, agent, env);
    println(""); println(""); println("");


    println("Hook - Post Experiment Stage")
    hookbefore = deepcopy(hook)
    hook(POST_EXPERIMENT_STAGE, agent, env);
    hook = deepcopy(hookbefore)
    @time hook(POST_EXPERIMENT_STAGE, agent, env);
    hook = deepcopy(hookbefore)
    @timeit to "hook_post_experiment" hook(POST_EXPERIMENT_STAGE, agent, env);
    println(""); println(""); println("");
end





env_cuda = false
agent_cuda = false


to = TimerOutput()
sections = []
results_time = DataFrame()
results_allocated = DataFrame()

nodes = [2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]

for i = 1:length(nodes)
    reset_timer!(to)
    memory_analysis(to, env_cuda, agent_cuda, nodes[i])
    collectresults(to, nodes[i], env_cuda)
end


layout = Layout(
            plot_bgcolor="#f1f3f7",
            title = "Timing Analysis of Different Sections",
            scene = attr(
                xaxis_title = "Size of Network",
                yaxis_title = "Time"),
            width = 1000,
            height = 650,
            margin=attr(l=10, r=10, b=10, t=60, pad=10),
            barmode="stack",
            #xaxis_categoryorder="category ascending"
        )

toplot = []

for node in nodes
    push!(toplot, bar(results_time, x=sections, y=Symbol(string(node)), name=string(node)))
end

toplot = Array{GenericTrace}(toplot)

p = plot(toplot, layout)
display(p)

layout["title"] = "Memory Analysis of Different Sections"
layout["scene"] = attr(
                    xaxis_title = "Size of Network",
                    yaxis_title = "Size")
toplot = []

for node in nodes
    push!(toplot, bar(results_allocated, x=sections, y=Symbol(string(node)), name=string(node)))
end

toplot = Array{GenericTrace}(toplot)

p = plot(toplot, layout)
display(p)