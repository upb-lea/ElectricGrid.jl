export MultiAgentGridController


#=
Example for agents dict:

Dict(nameof(agent) => {"policy" => agent,
                       "state_ids" => state_ids_agent,
                       "action_ids" => action_ids_agent},
    nameof(Animo) => {"policy" => Animo,
                      "state_ids" => state_ids_classic,
                      "action_ids" => action_ids_classic})
=#
mutable struct MultiAgentGridController <: AbstractPolicy
    agents::Dict{Any,Any}
    action_ids
    hook
end

function MultiAgentGridController(agents, action_ids)
    hook = data_hook(is_inner_hook_RL = true, plot_rewards = true)

    return MultiAgentGridController(
        agents,
        action_ids,
        hook
    )
end

Base.getindex(A::MultiAgentGridController, x) = getindex(A.agents, x)

function (A::MultiAgentGridController)(env::AbstractEnv, training::Bool = false)
    action = Array{Union{Nothing, Float64}}(nothing, length(A.action_ids))

    for agent in values(A.agents)
        action[findall(x -> x in agent["action_ids"], A.action_ids)] = agent["policy"](env, training)
    end

    return action
end

function (A::MultiAgentGridController)(stage::AbstractStage, env::AbstractEnv, training::Bool = false)
    A.hook(stage, A, env, training)

    if training
        for agent in values(A.agents)
            agent["policy"](stage, env, training)
        end
    end
end

function (A::MultiAgentGridController)(stage::PreActStage, env::AbstractEnv, action, training::Bool = false)
    A.hook(stage, A, env, action, training)

    if training
        for agent in values(A.agents)
            agent["policy"](stage, env, action[findall(x -> x in agent["action_ids"], A.action_ids)], training)
        end
    end
end

function reset_policy(A::MultiAgentGridController)
    for agent in values(A.agents)
        reset_policy(agent["policy"])
    end
end

function reset_policy(np::NamedPolicy) 
    reset_policy(np.policy)
end

function reset_policy(::AbstractPolicy) end