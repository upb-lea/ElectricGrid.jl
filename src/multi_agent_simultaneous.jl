export MultiAgentGridController

struct MultiAgentGridController <: AbstractPolicy
    agents::Dict{Any,Any} # -> [classicCtrn, RL,...]
    # TODO: belonging to specific agent order deine which action/state-ids belong to which agent using get_action/state_ids("source_x")
end

Base.getindex(A::MultiAgentGridController, x) = getindex(A.agents, x)

"""
MultiAgentGridController(source_controller => action...)

"""
MultiAgentGridController(policies...) =
    MultiAgentGridController(Dict{Any,Any}(nameof(p) => p for p in policies))

function (A::MultiAgentGridController)(env::AbstractEnv) 

    # loop throug all agents in A
    for agent in A.agents
        asd = 1 
    end

end

function (A::MultiAgentGridController)(stage::AbstractStage, env::AbstractEnv)
    for agent in values(A.agents)
        agent(stage, env)
    end
end

function (A::MultiAgentGridController)(stage::PreActStage, env::AbstractEnv, action)
    #TODO: pass all real agents their separate action set
    #agent(stage, env, action)
end