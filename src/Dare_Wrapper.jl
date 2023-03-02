"""
    setup_agents(env, hook; num_episodes = 1, return_Agents = false)

# Description
Initialises up the agents that will be controlling the electrical network.

# Arguments
- `env::SimEnv`: mutable struct containing the environment.

# Return Values
- `Multi_Agent::MultiAgentGridController`: the struct containing the initialised agents

"""
function setup_agents(env)


    #_______________________________________________________________________________
    # Setting up the Agents

    Agents = Dict()

    #-------------------------------------------------------------------------------
    # Setting up the controls for the Classical Sources

    Animo = NamedPolicy("classic", Classical_Policy(env))

    if !isnothing(Animo.policy)
        num_clas_sources = Animo.policy.Source.num_sources # number of classically controlled sources
    else
        num_clas_sources = 0
    end

    #-------------------------------------------------------------------------------
    # Setting up the controls for the Reinforcement Learning Sources

    num_RL_sources = env.nc.num_sources - num_clas_sources # number of reinforcement learning controlled sources

    RL_Source_Indices = Array{Int64, 1}(undef, 0)
    for ns in 1:env.nc.num_sources

        if env.nc.parameters["source"][ns]["control_type"] == "RL"
            RL_Source_Indices = [RL_Source_Indices; Int(ns)]
        end
    end

    ssa = "source".*string.(RL_Source_Indices)
    env.state_ids_RL = filter(x -> !isempty(findall(y -> y == split(x, "_")[1], ssa)), env.state_ids)
    env.action_ids_RL = filter(x -> !isempty(findall(y -> y == split(x, "_")[1], ssa)), env.action_ids)

    if num_RL_sources > 0

        if num_clas_sources > 0
            na = length(env.action_ids) - length(Animo.policy.action_ids)
        else
            na = length(env.action_ids)
        end

        agent = create_agent_ddpg(na = na, ns = length(state(env,"agent")), use_gpu = false)

        for i in 1:3
            agent.policy.behavior_actor.model.layers[i].weight ./= 100
            agent.policy.behavior_actor.model.layers[i].bias ./= 100
            agent.policy.target_actor.model.layers[i].weight ./= 100
            agent.policy.target_actor.model.layers[i].bias ./= 100
        end

        agent = Agent(policy = NamedPolicy("agent", agent.policy), trajectory = agent.trajectory)

        RL_policy = Dict()
        RL_policy["policy"] = agent
        RL_policy["state_ids"] = env.state_ids_RL
        RL_policy["action_ids"] = env.action_ids_RL

        Agents[nameof(agent)] = RL_policy
    end

    #-------------------------------------------------------------------------------
    # Finalising the Control

    num_policies = num_RL_sources + convert(Int64, num_clas_sources > 0)

    if num_clas_sources > 0

        polc = Dict()

        polc["policy"] = Animo
        polc["state_ids"] = Animo.policy.state_ids
        polc["action_ids"] = Animo.policy.action_ids

        Agents[nameof(Animo)] = polc
    end

    Multi_Agent = MultiAgentGridController(Agents, env.action_ids)

    #_______________________________________________________________________________
    # Logging

    if env.verbosity > 1

        @info "Time simulation run time: $((env.maxsteps - 1)*env.ts) [s] ~ $(Int(env.maxsteps)) steps"
    end

    #_______________________________________________________________________________
    # returns

    return Multi_Agent

end

function simulate(Multi_Agent, env; num_episodes = 1, hook = nothing)

    if isnothing(hook) # default hook

        hook = default_data_hook(Multi_Agent, env)

    end

    dare_run(Multi_Agent, env, StopAfterEpisode(num_episodes), hook)

    return hook
end

function learn(Multi_Agent, env; num_episodes = 1, hook = nothing)

    if isnothing(hook) # default hook

        hook = default_data_hook(Multi_Agent, env)

    end

    dare_run(Multi_Agent, env, StopAfterEpisode(num_episodes), hook, true)

    return hook
end

"""
which signals are the default ones if the user does not define a data hook for plotting
"""
function default_data_hook(Multi_Agent, env)

    Source = Multi_Agent.agents["classic"]["policy"].policy.Source
    all_class = collect(1:Source.num_sources)
    all_sources = collect(1:env.nc.num_sources)
    all_loads = collect(1:env.nc.num_loads)
    all_cables = collect(1:env.nc.num_connections)

    hook = data_hook(collect_sources  = all_sources,
                    collect_cables   = all_cables,
                    #collect_loads    = all_loads,
                    v_mag            = all_class,
                    i_mag            = all_class,
                    vdq              = all_class,
                    idq              = all_class,
                    power_pq         = all_class,
                    freq             = all_class,
                    angles           = all_class,
                    i_sat            = all_class,
                    v_sat            = Source.grid_forming,
                    i_err_t          = all_class,
                    v_err_t          = Source.grid_forming,
                    debug            = [])

    return hook
end

function dare_run(policy, env, stop_condition, hook, training = false)

    hook(PRE_EXPERIMENT_STAGE, policy, env, training)
    policy(PRE_EXPERIMENT_STAGE, env, training)
    is_stop = false
    while !is_stop
        RLBase.reset!(env)
        policy(PRE_EPISODE_STAGE, env, training)
        hook(PRE_EPISODE_STAGE, policy, env, training)

        while !is_terminated(env) # one episode
            action = policy(env, training)

            policy(PRE_ACT_STAGE, env, action, training)
            hook(PRE_ACT_STAGE, policy, env, action, training)

            env(action)

            policy(POST_ACT_STAGE, env, training)
            hook(POST_ACT_STAGE, policy, env, training)

            if stop_condition(policy, env)
                is_stop = true
                break
            end
        end # end of an episode

        if is_terminated(env)
            policy(POST_EPISODE_STAGE, env, training)  # let the policy see the last observation
            hook(POST_EPISODE_STAGE, policy, env, training)
        end
    end
    hook(POST_EXPERIMENT_STAGE, policy, env, training)
    hook
end