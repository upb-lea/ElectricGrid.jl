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

    insertcols!(hook.tmp, :episode => hook.ep)

    for state_id in hook.state_ids
        insertcols!(hook.tmp, state_id => env.state[findfirst(x -> x == state_id, env.state_ids)])
    end

    insertcols!(hook.tmp, :action => Ref(action))
    
end

function (hook::DataHook)(::PostActStage, agent, env)

   # insertcols!(hook.tmp, :next_state => Ref(env.state))
    insertcols!(hook.tmp, :reward => env.reward)
    insertcols!(hook.tmp, :done => env.done)

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