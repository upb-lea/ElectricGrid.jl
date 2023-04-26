#=
Example for agents dict:

Dict(nameof(agent) => {"policy" => agent,
                       "state_ids" => state_ids_agent,
                       "action_ids" => action_ids_agent},
    nameof(Animo) => {"policy" => Animo,
                      "state_ids" => state_ids_classic,
                      "action_ids" => action_ids_classic})
=#
mutable struct MultiController <: AbstractPolicy
    agents::Dict{Any,Any}
    action_ids
    hook
end

"""
MultiController(agents, action_ids)

# Description
Constructs the MultiController providing different subfunctions mappen the correct env
indices to the corresponding controllers and agents.

# Arguments
- `agents::Dict(NamedPolicy, ids)`: Contains all the agents as named policys with there
    corresponding ids.
- `action_ids::Vector{String}`: All env action indices.

# Return Values
- `Multi_Agent::MultiController`: the struct containing the initialised agents

"""
function MultiController(agents, action_ids)
    hook = DataHook(is_inner_hook_RL = true, plot_rewards = true)

    return MultiController(
        agents,
        action_ids,
        hook
    )
end

Base.getindex(A::MultiController, x) = getindex(A.agents, x)

function (A::MultiController)(env::AbstractEnv, training::Bool = false)
    action = Array{Union{Nothing, Float64}}(nothing, length(A.action_ids))

    for agent in values(A.agents)
        action[findall(x -> x in agent["action_ids"], A.action_ids)] .= agent["policy"](env, training)
    end

    return action
end

function (A::MultiController)(stage::AbstractStage, env::AbstractEnv, training::Bool = false)
    A.hook(stage, A, env, training)

    if training
        for agent in values(A.agents)
            agent["policy"](stage, env, training)
        end
    end
end

function (A::MultiController)(stage::PreActStage, env::AbstractEnv, action, training::Bool = false)
    A.hook(stage, A, env, action, training)

    if training
        for agent in values(A.agents)
            agent["policy"](stage, env, action[findall(x -> x in agent["action_ids"], A.action_ids)], training)
        end
    end
end

function ResetPolicy(A::MultiController)
    for agent in values(A.agents)
        ResetPolicy(agent["policy"])
    end
end

function ResetPolicy(np::NamedPolicy)
    ResetPolicy(np.policy)
end

function ResetPolicy(::AbstractPolicy) end


"""
    SetupAgents(env, custom_agents)

# Description
Initialises up the agents that will be controlling the electrical network.
It depends on the control_type and mode in the parameter dict.
Handles to connect the different controllers to the correcponding sources represended by
there indices (for states and actions).

# Arguments
- `env::ElectricGridEnv`: mutable struct containing the environment.
- `custom_agents::Dict{Agent}`: Dict of Agents. The key is used as name in the MultiController.

# Return Values
- `Multi_Agent::MultiController`: the struct containing the initialised agents

"""
function SetupAgents(env, custom_agents = nothing)


    #_______________________________________________________________________________
    # Setting up the Agents

    Agents = Dict()

    #-------------------------------------------------------------------------------
    # Setting up the controls for the Classical Sources

    Animo = NamedPolicy("classic", ClassicalPolicy(env))

    if !isnothing(Animo.policy)
        num_clas_sources = Animo.policy.Source.num_sources # number of classically controlled sources
    else
        num_clas_sources = 0
    end

    #-------------------------------------------------------------------------------
    # Setting up the controls for the Reinforcement Learning Sources

    for (name, config) in env.agent_dict

        if config["mode"] == "ElectricGrid_ddpg"
            agent = CreateAgentDdpg(na = length(env.agent_dict[name]["action_ids"]),
                                        ns = length(state(env, name)),
                                        use_gpu = false)

            for i in 1:3
                agent.policy.behavior_actor.model.layers[i].weight ./= 3
                agent.policy.behavior_actor.model.layers[i].bias ./= 3
                agent.policy.target_actor.model.layers[i].weight ./= 3
                agent.policy.target_actor.model.layers[i].bias ./= 3
            end
        else
            agent = custom_agents[name]
        end

        agent = Agent(policy = NamedPolicy(name, agent.policy), trajectory = agent.trajectory)

        RL_policy = Dict()
        RL_policy["policy"] = agent
        RL_policy["state_ids"] = env.agent_dict[name]["state_ids"]
        RL_policy["action_ids"] = env.agent_dict[name]["action_ids"]

        Agents[name] = RL_policy
    end

    @eval function action_space(env::ElectricGridEnv, name::String)
        for (pname, config) in env.agent_dict
            if name == pname
                return Space(fill(-1.0..1.0, length(env.agent_dict[name]["action_ids"])))
            end
        end

        return Space(fill(-1.0..1.0, size(env.action_ids_RL)))
    end
    #-------------------------------------------------------------------------------
    # Finalising the Control

    #num_policies = num_RL_sources + convert(Int64, num_clas_sources > 0)

    if num_clas_sources > 0

        polc = Dict()

        polc["policy"] = Animo
        polc["state_ids"] = Animo.policy.state_ids
        polc["action_ids"] = Animo.policy.action_ids

        Agents[nameof(Animo)] = polc
    end

    Multi_Agent = MultiController(Agents, env.action_ids)

    #_______________________________________________________________________________
    # Logging


    #_______________________________________________________________________________
    # returns

    return Multi_Agent

end


"""
hook = Simulate(controllers, env, num_episodes = 1 hook=hook)

# Description
Executed the given controllers in the environments for num_episodes. Fills thereby the hook.
In case of RL agents deterministic actions are chosen, e.g., no action noise is used.

# Arguments
- `controllers::MultiController`: Contains controllers and agents
- `env::ElectricGridEnv`: mutable struct containing the environment.
- `num_episodes::Int`: Number of episodes simulated.
- `hook::DataHook`: Stores the data to measure in a DataFrame.

# Return Values
- `hook::DataHook`: Measured data.

"""
function Simulate(Multi_Agent, env; num_episodes = 1, hook = nothing)

    if isnothing(hook) # default hook

        hook = DefaultDataHook(Multi_Agent, env)

    end

    CustomRun(Multi_Agent, env, StopAfterEpisode(num_episodes), hook)

    return hook
end

"""
hook = Learn(controllers, env, num_episodes = 1 hook=hook)

# Description
Executed the given controllers in the environments for num_episodes. Fills thereby the hook.
Here the RL agents are training and e.g. action noise is applied.

# Arguments
- `controllers::MultiController`: Contains controllers and agents
- `env::ElectricGridEnv`: mutable struct containing the environment.
- `num_episodes::Int`: Number of episodes simulated.
- `hook::DataHook`: Stores the data to measure in a DataFrame.

# Return Values
- `hook::DataHook`: Measured data.

"""
function Learn(Multi_Agent, env; num_episodes = 1, hook = nothing)

    if isnothing(hook) # default hook

        hook = DataHook()

    end

    CustomRun(Multi_Agent, env, StopAfterEpisode(num_episodes), hook, true)

    return hook
end


"""
which signals are the default ones if the user does not define a data hook for plotting
"""
function DefaultDataHook(Multi_Agent, env)

    Source = Multi_Agent.agents["classic"]["policy"].policy.Source
    all_class = collect(1:Source.num_sources)
    all_sources = collect(1:env.nc.num_sources)
    all_loads = collect(1:env.nc.num_loads)
    all_cables = collect(1:env.nc.num_connections)

    hook = DataHook(collect_sources  = all_sources,
                    collect_cables   = all_cables,
                    #collect_loads    = all_loads,
                    v_mag_inv        = all_class,
                    v_mag_poc        = all_class,
                    i_mag_inv        = all_class,
                    i_mag_poc        = all_class,
                    v_dq             = all_class,
                    i_dq             = all_class,
                    power_pq_inv     = all_class,
                    power_pq_poc     = all_class,
                    freq             = all_class,
                    angles           = all_class,
                    i_sat            = all_class,
                    v_sat            = Source.grid_forming,
                    i_err_t          = all_class,
                    v_err_t          = Source.grid_forming,
                    debug            = [])

    return hook
end

"""
Wrapps the Run function form https://juliareinforcementlearning.org/ to enable turning off
the action noise.
"""
function CustomRun(policy, env, stop_condition, hook, training = false)

    @timeit to "Hook" begin
        hook(PRE_EXPERIMENT_STAGE, policy, env, training)
    end
    @timeit to "Policy" begin
        @timeit to "Updates" begin
            policy(PRE_EXPERIMENT_STAGE, env, training)
        end
    end
    is_stop = false
    while !is_stop
        @timeit to "Hook" begin
            RLBase.reset!(env)
        end
        @timeit to "Policy" begin
            @timeit to "Updates" begin
                ResetPolicy(policy)

                policy(PRE_EPISODE_STAGE, env, training)
            end
        end
        @timeit to "Hook" begin
            hook(PRE_EPISODE_STAGE, policy, env, training)
        end

        while !is_terminated(env) # one episode
            @timeit to "Policy" begin
                @timeit to "Action" begin
                    action = policy(env, training)
                end
                @timeit to "Updates" begin
                    policy(PRE_ACT_STAGE, env, action, training)
                end
            end
            @timeit to "Hook" begin
                hook(PRE_ACT_STAGE, policy, env, action, training)
            end

            @timeit to "Env" begin
                env(action)
            end

            @timeit to "Policy" begin
                @timeit to "Updates" begin
                    policy(POST_ACT_STAGE, env, training)
                end
            end
            @timeit to "Hook" begin
                hook(POST_ACT_STAGE, policy, env, training)
            end

            if stop_condition(policy, env)
                is_stop = true
                break
            end
        end # end of an episode

        if is_terminated(env)
            @timeit to "Policy" begin
                @timeit to "Updates" begin
                    policy(POST_EPISODE_STAGE, env, training)  # let the policy see the last observation
                end
            end
            @timeit to "Hook" begin
                hook(POST_EPISODE_STAGE, policy, env, training)
            end
        end
    end

    @timeit to "Policy" begin
        @timeit to "Updates" begin
            policy(POST_EXPERIMENT_STAGE, env, training)
        end
    end
    @timeit to "Hook" begin
        hook(POST_EXPERIMENT_STAGE, policy, env, training)
    end
    return hook
end
