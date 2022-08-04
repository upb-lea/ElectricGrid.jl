using ReinforcementLearning
using DataFrames
using CSV

Base.@kwdef mutable struct Summary <: AbstractHook

    save = false
    save_path = "episode_data/summary.csv"

    rewards::Vector{Float64} = Float64[]
    ep_rewards::Vector{Float64} = Float64[]

end

function (hook::Summary)(::PostActStage, agent, env)

    append!(hook.ep_rewards, env.reward)
    
end

function (hook::Summary)(::PostEpisodeStage, agent, env)

    append!(hook.rewards, mean(hook.ep_rewards))
    hook.ep_rewards = Float64[]

end

function (hook::Summary)(::PostExperimentStage, agent, env)

    if save
        summary = DataFrame(:MeanReward => hook.rewards)
        CSV.write("episode_data/summary.csv", summary)
    end

end