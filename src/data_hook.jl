using ReinforcementLearning
using DataFrames

Base.@kwdef mutable struct DataHook <: AbstractHook

    save_data_to_hd = false
    dir = "episode_data/"

    state_ids = []

    df = DataFrame()
    tmp = DataFrame()
    ep = 1

end

function (hook::DataHook)(::PreExperimentStage, agent, env)

    # rest
    hook.df = DataFrame()
    hook.ep = 1
    
end

function (hook::DataHook)(::PreActStage, agent, env, action)

    insertcols!(hook.tmp, 1, :episode => hook.ep)
    insertcols!(hook.tmp, 2, :state => Ref(env.state))
    #TODO actions implementieren
    #insertcols!(hook.tmp, 3, :action => Ref(action))
    
end

function (hook::DataHook)(::PostActStage, agent, env)

    insertcols!(hook.tmp, 4, :next_state => Ref(env.state))
    insertcols!(hook.tmp, 5, :reward => env.reward)
    insertcols!(hook.tmp, 6, :done => env.done)

    append!(hook.df, hook.tmp)
    hook.tmp = DataFrame()
    
end

function (hook::DataHook)(::PostEpisodeStage, agent, env)

    hook.ep += 1

end

function (hook::DataHook)(::PostExperimentStage, agent, env)

    if save_data_to_hd
        isdir(hook.dir) || mkdir(hook.dir)
        Arrow.write(hook.dir * "data.arrow", hook.df)
    end

end