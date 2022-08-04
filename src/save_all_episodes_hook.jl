using ReinforcementLearning
using DataFrames
using CSV

Base.@kwdef mutable struct SaveAllEpisodes <: AbstractHook

    dir = "episode_data"

    df = DataFrame()
    tmp = DataFrame()
    episode = 1

    rewards::Vector{Float64} = Float64[]

end

function (hook::SaveAllEpisodes)(::PreActStage, agent, env, action)

    insertcols!(hook.tmp, 1, :state => Ref(env.state))
    insertcols!(hook.tmp, 2, :action => Ref(action))
    
end

function (hook::SaveAllEpisodes)(::PostActStage, agent, env)

    insertcols!(hook.tmp, 3, :next_state => Ref(env.state))
    insertcols!(hook.tmp, 4, :reward => env.reward)
    insertcols!(hook.tmp, 5, :done => env.done)

    append!(hook.df, hook.tmp)

    hook.tmp = DataFrame();
    
end

function (hook::SaveAllEpisodes)(::PostEpisodeStage, agent, env)

    push!(hook.rewards, mean(hook.df.reward))

    CSV.write(dir + "/$(hook.episode).csv", hook.df)
    hook.df = DataFrame()
    hook.episode += 1

end

function (hook::SaveAllEpisodes)(::PostExperimentStage, agent, env)

    hook.episode = 0

end