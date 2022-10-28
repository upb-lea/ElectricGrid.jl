using ReinforcementLearning

mutable struct dare_setup
    env
    control_agent
    hook
end

function create_setup(;num_sources = 0, num_loads = 0, CM = nothing, parameters = nothing, 
                        t_end = 0.5, env_cuda = false, agent_cuda = false, hook = nothing)

    #-------------------------------------------------------------------------------
    # Auxiliary Functions

    function reward(env, name = nothing)
        u_l1_index = findfirst(x -> x == "source1_v_C", env.state_ids)

        u_l1 = env.state[u_l1_index]

        ref = sqrt(2)*230 * cos.(2*pi*50*env.t)

        r = -abs(ref/600 - u_l1)/3 * (1 - 0.99)
    
        return r
    end

    function featurize(x0 = nothing, t0 = nothing; env = nothing, name = nothing)

        if !isnothing(name)
            state = env.state
            if name == agentname
                global state_ids_agent
                global state_ids
                state = state[findall(x -> x in state_ids_agent, state_ids)]
                state = vcat(state, reference(env.t)/600)
            else
                global state_ids_classic
                global state_ids
                state = env.x[findall(x -> x in state_ids_classic, state_ids)]
            end
        elseif isnothing(env)
            return x0
        else
            return env.state
        end
        return state
    end

    function RLBase.action_space(env::SimEnv, name::String)
        if name == "agent"
            return Space(fill(-1.0..1.0, size(action_ids_agent)))
        else
            return Space(fill(-1.0..1.0, size(action_ids_classic)))
        end
    end

    #-------------------------------------------------------------------------------
    # Finalising the Node Constructor

    if !isempty(parameters)
        num_sources = length(source_list)
        num_loads = length(load_list)
    end

    if parameters["grid"]["fs"] === nothing
        ts = 100e-6
        fs = 1/ts
    else
        ts = 1/fs
    end
        
    maxsteps = Int(t_end / ts) + 1

    #-------------------------------------------------------------------------------
    # Defining the Environment

    env = SimEnv(reward_function = reward, num_sources = num_sources, num_loads = num_loads, 
                CM = CM, parameters = parameters, ts = ts, 
                maxsteps = maxsteps, use_gpu = env_cuda, action_delay = 1)
    
    state_ids = get_state_ids(env.nc)
    action_ids = get_action_ids(env.nc)

    #-------------------------------------------------------------------------------
    # Outputs
    
    if hook === nothing
        hook = DataHook(save_best_NNA = true, plot_rewards = true, 
        collect_state_ids = env.state_ids, collect_action_ids = [])
    end

    #-------------------------------------------------------------------------------
    # Setting up the controls for the Classical Sources

    Animo = NamedPolicy("classic", Classical_Policy(env))

    num_clas_sources = length(Animo.policy.Source_Indices) # number of classically controlled sources
    state_ids_classic = Animo.policy.state_ids
    action_ids_classic = Animo.policy.action_ids

    #-------------------------------------------------------------------------------
    # Setting up the controls for the Reinforcement Learning Sources

    num_RL_sources = env.nc.num_sources - num_clas_sources # number of reinforcement learning controlled sources

    if num_RL_sources == 1

        ns = length(env.sys_d.A[1,:])
        na = length(env.sys_d.B[1,:])

        agent = create_agent_ddpg(na = na, ns = ns, use_gpu = env_cuda)
        agent = Agent(policy = NamedPolicy("agent", agent.policy), trajectory = agent.trajectory)
    end

    #-------------------------------------------------------------------------------
    # Finalising the Control

    num_policies = num_RL_sources + convert(Int64, num_clas_sources > 0)

    Multi_Agents = Dict()

    if num_clas_sources > 0

        polc = Dict()

        polc["policy"] = Animo
        polc["state_ids"] = state_ids_classic
        polc["action_ids"] = action_ids_classic

        Multi_Agents[nameof(Animo)] = polc
    end
    
    Multi_Agent = MultiAgentGridController(Multi_Agents, action_ids)

    return dare_setup(env, Multi_Agent, hook)
end

function run(dare_setup, no_episodes)

    RLBase.run(dare_setup.control_agent, dare_setup.env, StopAfterEpisode(no_episodes), dare_setup.hook)

    return nothing
end