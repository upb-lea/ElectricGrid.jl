using ReinforcementLearning
using DataFrames
using CSV

Base.@kwdef mutable struct SaveAllEpisodes <: AbstractHook

    dir = "episode_data"

    df = DataFrame()
    tmp = DataFrame()
    ep = 1

end

function (hook::SaveAllEpisodes)(::PreExperimentStage, agent, env)

    # rest
    hook.df = DataFrame()
    hook.ep = 1
    
end

function (hook::SaveAllEpisodes)(::PreActStage, agent, env, action)

    insertcols!(hook.tmp, 1, :episode => hook.ep)
    insertcols!(hook.tmp, 2, :state => Ref(env.state))
    insertcols!(hook.tmp, 3, :action => Ref(action))
    
end

function (hook::SaveAllEpisodes)(::PostActStage, agent, env)

    insertcols!(hook.tmp, 4, :next_state => Ref(env.state))
    insertcols!(hook.tmp, 5, :reward => env.reward)
    insertcols!(hook.tmp, 6, :done => env.done)

    append!(hook.df, hook.tmp)
    hook.tmp = DataFrame()
    
end

function (hook::SaveAllEpisodes)(::PostEpisodeStage, agent, env)

    hook.ep += 1

end

function (hook::SaveAllEpisodes)(::PostExperimentStage, agent, env)

    CSV.write(hook.dir * "/test.csv", hook.df)

end