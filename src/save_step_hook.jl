using ReinforcementLearning
using DataFrames
using CSV

Base.@kwdef mutable struct SaveStep <: AbstractHook

    df = DataFrame()
    tmp = DataFrame()
    episode = 1

    rewards::Vector{Float64} = Float64[]

end

function (hook::SaveStep)(::PreActStage, agent, env, action)

    insertcols!(hook.tmp, 1, :state => Ref(env.state))
    insertcols!(hook.tmp, 2, :action => Ref(action))

    # push!(hook.state, env.state)
    # push!(hook.action, action)
    
end

function (hook::SaveStep)(::PostActStage, agent, env)

    insertcols!(hook.tmp, 3, :next_state => Ref(env.state))
    insertcols!(hook.tmp, 4, :reward => env.reward)
    insertcols!(hook.tmp, 5, :done => env.done)

    append!(hook.df, hook.tmp)

    hook.tmp = DataFrame();
    
end

function (hook::SaveStep)(::PostEpisodeStage, agent, env)

    
    push!(hook.rewards, mean(hook.df.reward))

    CSV.write("episode_data/$(hook.episode).csv", hook.df)
    hook.df = DataFrame()
    hook.episode += 1

end

function (hook::SaveStep)(::PostExperimentStage, agent, env)

    summary = DataFrame(:MeanReward => hook.rewards)
    CSV.write("episode_data/summary.csv", summary)
    hook.episode = 0

end