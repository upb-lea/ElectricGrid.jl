using ReinforcementLearning
using DataFrames
using UnicodePlots

Base.@kwdef mutable struct DataHook <: AbstractHook

    save_data_to_hd = false
    dir = "episode_data/"

    state_ids = []
    next_state_ids = []

    df = DataFrame()
    tmp = DataFrame()
    ep = 1

    rewards::Vector{Float64} = Float64[]
    reward::Float64 = 0.0
    bestNNA = nothing
    bestreward = -1000000.0
    currentNNA = nothing
end

function (hook::DataHook)(::PreExperimentStage, agent, env)

    # rest
    #hook.df = DataFrame()
    #hook.ep = 1
    
    if hook.currentNNA === nothing
        hook.currentNNA = deepcopy(agent.policy.behavior_actor)
        hook.bestNNA = deepcopy(agent.policy.behavior_actor)
    end
end

function (hook::DataHook)(::PreActStage, agent, env, action)

    insertcols!(hook.tmp, :episode => hook.ep)
    insertcols!(hook.tmp, :time => Float32(env.t))

    for state_id in hook.state_ids
        state_index = findfirst(x -> x == state_id, env.state_ids)

        insertcols!(hook.tmp, state_id => (env.state[state_index] * env.norm_array[state_index]))
    end

    insertcols!(hook.tmp, :action => Ref(action))
    
end

function (hook::DataHook)(::PostActStage, agent, env)

    for state_id in hook.next_state_ids
        state_index = findfirst(x -> x == state_id, env.state_ids)

        insertcols!(hook.tmp, ("next_state_" * state_id) => (env.state[state_index] * env.norm_array[state_index]))
    end

    insertcols!(hook.tmp, :reward => env.reward)
    insertcols!(hook.tmp, :done => env.done)

    append!(hook.df, hook.tmp)
    hook.tmp = DataFrame()
    
    hook.reward += env.reward
end

function (hook::DataHook)(::PostEpisodeStage, agent, env)

    hook.ep += 1

    if length(hook.rewards) >= 1 && hook.reward > maximum(hook.rewards)
        copyto!(hook.bestNNA, agent.policy.behavior_actor)
        hook.bestreward = hook.reward
    end

    push!(hook.rewards, hook.reward)
    hook.reward = 0

    copyto!(hook.currentNNA, agent.policy.behavior_actor)

end

function (hook::DataHook)(::PostExperimentStage, agent, env)

    if hook.save_data_to_hd
        isdir(hook.dir) || mkdir(hook.dir)
        Arrow.write(hook.dir * "data.arrow", hook.df)
    end

    println(lineplot(hook.rewards, title="Total reward per episode", xlabel="Episode", ylabel="Score"))

end