using ReinforcementLearning
using DataFrames
using UnicodePlots

include(srcdir("plotting.jl"))

Base.@kwdef mutable struct DataHook <: AbstractHook

    save_data_to_hd = false
    dir = "episode_data/"
    A = nothing
    B= nothing
    collect_state_paras = nothing

    collect_sources = []
    collect_cables = []
    collect_loads = []


    collect_state_ids = []
    collect_next_state_ids = []
    collect_action_ids = []

    df = DataFrame()
    tmp = DataFrame()
    ep = 1

    plot_rewards = false
    rewards::Vector{Float64} = Float64[]
    reward::Float64 = 0.0

    save_best_NNA = false
    bestNNA = nothing
    bestreward = -1000000.0
    bestepisode = 0
    currentNNA = nothing

    collect_reference = false

end

function (hook::DataHook)(::PreExperimentStage, agent, env)

    # rest
    #hook.df = DataFrame()
    #hook.ep = 1
    for source_id in hook.collect_sources
        indizes = findall(x -> occursin("_f$source_id"*"_", x) ||occursin("_$source_id"*"_", x) || endswith(x , "_f$source_id") || endswith(x , "_$source_id"), env.state_ids)
        for id in indizes
            if !(env.state_ids[id] in hook.collect_state_ids)
                push!(hook.collect_state_ids,env.state_ids[id])
            end
        end
    end

    for cable_id in hook.collect_cables
        indizes = findall(x ->occursin("i_cable$cable_id"*"_", x) || x == "i_cable$cable_id", env.state_ids)
        for id in indizes
            if !(env.state_ids[id] in hook.collect_state_ids)
                push!(hook.collect_state_ids,env.state_ids[id])
            end
        end
    end

    for load_id in hook.collect_loads
        indizes = findall(x -> occursin("_load$load_id"*"_", x) || endswith(x , "_load$load_id"), env.state_ids)
        for id in indizes
            if !(env.state_ids[id] in hook.collect_state_ids)
                push!(hook.collect_state_ids,env.state_ids[id])
            end
        end
    end

    
    hook.A,hook.B ,_ ,_ = get_sys(env.nc)
    hook.collect_state_paras = get_state_paras(env.nc)

    if hook.save_best_NNA && hook.currentNNA === nothing
        hook.currentNNA = deepcopy(agent.policy.behavior_actor)
        hook.bestNNA = deepcopy(agent.policy.behavior_actor)
    end
end

function (hook::DataHook)(::PreActStage, agent, env, action)

    insertcols!(hook.tmp, :episode => hook.ep)
    insertcols!(hook.tmp, :time => Float32(env.t))

    if hook.collect_reference
        #insertcols!(hook.tmp, :reference => reference(env.t))
        refs = reference(env.t)
        for i = 1:length(refs)
            insertcols!(hook.tmp, "reference_$i" => refs[i])
        end
    end

    states_x = Vector( env.x )
    opstates=(hook.A * states_x + hook.B * (Vector(env.action)) ) .* (hook.collect_state_paras)

    for state_id in hook.collect_state_ids
        state_index = findfirst(x -> x == state_id, env.state_ids)

        insertcols!(hook.tmp, state_id => (env.x[state_index]))
        insertcols!(hook.tmp, "op_$state_id" => opstates[state_index,1])
    end 

    for action_id in hook.collect_action_ids
        action_index = findfirst(x -> x == action_id, env.action_ids)

        insertcols!(hook.tmp, action_id => (env.action[action_index]))
    end
    
    insertcols!(hook.tmp, :action => Ref(action))
    
end

function (hook::DataHook)(::PostActStage, agent, env)

    for state_id in hook.collect_next_state_ids
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

    if length(hook.rewards) >= 1 && hook.reward > maximum(hook.rewards)
        if hook.save_best_NNA
            copyto!(hook.bestNNA, agent.policy.behavior_actor)
        end
        hook.bestepisode = hook.ep
        hook.bestreward = hook.reward
    end

    hook.ep += 1

    push!(hook.rewards, hook.reward)
    hook.reward = 0

    if hook.save_best_NNA
        copyto!(hook.currentNNA, agent.policy.behavior_actor)
    end

end

function (hook::DataHook)(::PostExperimentStage, agent, env)

    if hook.save_data_to_hd
        isdir(hook.dir) || mkdir(hook.dir)
        Arrow.write(hook.dir * "data.arrow", hook.df)
    end

    if hook.plot_rewards
        println(lineplot(hook.rewards, title="Total reward per episode", xlabel="Episode", ylabel="Score"))
    end

end