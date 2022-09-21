using ReinforcementLearning
using DataFrames
using UnicodePlots

include(srcdir("plotting.jl"))
include(srcdir("MultiAgentGridController.jl"))

Base.@kwdef mutable struct DataHook <: AbstractHook

    save_data_to_hd = false
    dir = "episode_data/"
    A = nothing
    B= nothing
    collect_state_paras = nothing
    extra_state_paras= []
    extra_state_ids= []
    extra_state_names = []

    
    
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
    rewards::Vector{Vector{Float64}} = []
    reward::Vector{Float64} = [0.0]
    policy_names::Vector{String} = []

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
        para = env.nc.parameters["source"][source_id]
        indizes = findall(x -> occursin("source$source_id"*"_", x), env.state_ids)
        for id in indizes
            if !(env.state_ids[id] in hook.collect_state_ids)
                push!(hook.collect_state_ids,env.state_ids[id])
                if occursin("_L1", env.state_ids[id])
                    push!(hook.extra_state_ids,id)
                    push!(hook.extra_state_paras,para["R1"])
                    push!(hook.extra_state_names,replace(env.state_ids[id], "_L1" => "_R1"))
                elseif occursin("u_C", env.state_ids[id])
                    push!(hook.extra_state_ids,id)
                    push!(hook.extra_state_paras,para["R_C"])
                    push!(hook.extra_state_names, replace(env.state_ids[id], "u_C" => "i_R_C"))
                elseif occursin("_L2", env.state_ids[id])
                    push!(hook.extra_state_ids,id)
                    push!(hook.extra_state_paras,para["R2"])
                    push!(hook.extra_state_names, replace(env.state_ids[id], "_L2" => "_R2"))
                end
            end
        end
    end

    for cable_id in hook.collect_cables
        para = env.nc.parameters["cable"][cable_id]
        indizes = findall(x ->occursin("cable$cable_id"*"_", x), env.state_ids)
        for id in indizes
            if !(env.state_ids[id] in hook.collect_state_ids)
                push!(hook.collect_state_ids,env.state_ids[id])
                push!(hook.extra_state_ids,id)
                push!(hook.extra_state_paras,para["R"])
                push!(hook.extra_state_names, replace(env.state_ids[id], "_L" => "_R"))
            end
        end
    end

    for load_id in hook.collect_loads
        para = env.nc.parameters["load"][load_id]
        indizes = findall(x -> occursin("load$load_id"*"_", x), env.state_ids)
        for id in indizes
            if !(env.state_ids[id] in hook.collect_state_ids)
                push!(hook.collect_state_ids,env.state_ids[id])
                if occursin("_u", env.state_ids[id])
                    if occursin("R", para["impedance"]) && occursin("C", para["impedance"])
                        push!(hook.extra_state_ids,id)
                        push!(hook.extra_state_paras,(para["R"],(para["C"])*(para["C"]+get_C_sum_cable_node(env.nc.num_sources+load_id,env.nc))^(-1),(get_C_sum_cable_node(env.nc.num_sources+load_id,env.nc))*(para["C"]+get_C_sum_cable_node(env.nc.num_sources+load_id,env.nc))^(-1)), )
                        push!(hook.extra_state_names,(replace(env.state_ids[id], "_u_C_total" => "_i_R"),replace(env.state_ids[id], "_u_C_total" => "_i_C"),replace(env.state_ids[id], "_u_C_total" => "_i_C_cables")))
                    elseif occursin("R", para["impedance"])
                        push!(hook.extra_state_ids,id)
                        push!(hook.extra_state_paras,(para["R"],0,(1)))
                        push!(hook.extra_state_names,(replace(env.state_ids[id], "_u_C_total" => "_i_R"),replace(env.state_ids[id], "_u_C_total" => "_i_C_cables")))
                    elseif occursin("C", para["impedance"])
                        push!(hook.extra_state_ids,id)
                        push!(hook.extra_state_paras,(0,(para["C"])*(para["C"]+get_C_sum_cable_node(env.nc.num_sources+load_id,env.nc))^(-1),(get_C_sum_cable_node(env.nc.num_sources+load_id,env.nc))*(para["C"]+get_C_sum_cable_node(env.nc.num_sources+load_id,env.nc))^(-1) ))
                        push!(hook.extra_state_names,(replace(env.state_ids[id], "_u_C_total" => "_i_C"),replace(env.state_ids[id], "_u_C_total" => "_i_C_cables")))
                    else
                        push!(hook.extra_state_ids,id)
                        push!(hook.extra_state_paras,(0,0,1))
                        push!(hook.extra_state_names,replace(env.state_ids[id], "_u_C_total" => "_i_C_cables"))
                    end
                end
            end
        end
    end

    print(hook.extra_state_names)
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
    extra_state_cntr= 1
    for state_id in hook.collect_state_ids
        state_index = findfirst(x -> x == state_id, env.state_ids)
        
        insertcols!(hook.tmp, state_id => (env.x[state_index]))
        insertcols!(hook.tmp, replace(state_id, "_i_" => "_u_", "_u_" => "_i_") => opstates[state_index,1])

        if state_index in hook.extra_state_ids
            if occursin("load",state_id)
                if hook.extra_state_paras[extra_state_cntr][1] != 0
                    insertcols!(hook.tmp, hook.extra_state_names[extra_state_cntr][1] => (env.x[state_index])*hook.extra_state_paras[extra_state_cntr][1]^(-1))
                    insertcols!(hook.tmp, replace(hook.extra_state_names[extra_state_cntr][1], "_i_" => "_u_") => (env.x[state_index]))
                end
                if  hook.extra_state_paras[extra_state_cntr][2] != 0
                    insertcols!(hook.tmp, hook.extra_state_names[extra_state_cntr][2] => (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][2])
                    insertcols!(hook.tmp, replace(hook.extra_state_names[extra_state_cntr][2], "_i_" => "_u_") => (env.x[state_index]))
                end
                if  hook.extra_state_paras[extra_state_cntr][3] != 0
                    insertcols!(hook.tmp, hook.extra_state_names[extra_state_cntr][3] => (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][3])
                    insertcols!(hook.tmp, replace(hook.extra_state_names[extra_state_cntr][3], "_i_" => "_u_") => (env.x[state_index]))
                end
            elseif occursin("_i_", state_id)
                insertcols!(hook.tmp, hook.extra_state_names[extra_state_cntr] => (env.x[state_index]))
                insertcols!(hook.tmp, replace(hook.extra_state_names[extra_state_cntr], "_i_" => "_u_") => (env.x[state_index])*hook.extra_state_paras[extra_state_cntr])
            else
                insertcols!(hook.tmp, hook.extra_state_names[extra_state_cntr] => opstates[state_index,1])
                insertcols!(hook.tmp, replace(hook.extra_state_names[extra_state_cntr], "_i_" => "_u_") => (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr])
            end
        extra_state_cntr+=1
        end
    end 

    for action_id in hook.collect_action_ids
        action_index = findfirst(x -> x == action_id, env.action_ids)

        insertcols!(hook.tmp, action_id => (env.action[action_index]))
    end
    
    insertcols!(hook.tmp, :action => Ref(action))
    
end

function (hook::DataHook)(::PostActStage, agent, env)


    states_x = Vector( env.x )
    opstates=(hook.A * states_x + hook.B * (Vector(env.action)) ) .* (hook.collect_state_paras)
    extra_state_cntr= 1
    for state_id in hook.collect_state_ids
        state_index = findfirst(x -> x == state_id, env.state_ids)
        
        insertcols!(hook.tmp, "next_state_"*state_id => (env.x[state_index]))
        insertcols!(hook.tmp, "next_state_"*replace(state_id, "_i_" => "_u_", "_u_" => "_i_") => opstates[state_index,1])

        if state_index in hook.extra_state_ids
            if occursin("load",state_id)
                if hook.extra_state_paras[extra_state_cntr][1] != 0
                    insertcols!(hook.tmp, "next_state_"*hook.extra_state_names[extra_state_cntr][1] => (env.x[state_index])*hook.extra_state_paras[extra_state_cntr][1]^(-1))
                    insertcols!(hook.tmp, "next_state_"*replace(hook.extra_state_names[extra_state_cntr][1], "_i_" => "_u_") => (env.x[state_index]))
                end
                if  hook.extra_state_paras[extra_state_cntr][2] != 0
                    insertcols!(hook.tmp, "next_state_"*hook.extra_state_names[extra_state_cntr][2] => (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][2])
                    insertcols!(hook.tmp, "next_state_"*replace(hook.extra_state_names[extra_state_cntr][2], "_i_" => "_u_") => (env.x[state_index]))
                end
                if  hook.extra_state_paras[extra_state_cntr][3] != 0
                    insertcols!(hook.tmp, "next_state_"*hook.extra_state_names[extra_state_cntr][3] => (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][3])
                    insertcols!(hook.tmp, "next_state_"*replace(hook.extra_state_names[extra_state_cntr][3], "_i_" => "_u_") => (env.x[state_index]))
                end
            elseif occursin("_i_", state_id)
                insertcols!(hook.tmp, "next_state_"*hook.extra_state_names[extra_state_cntr] => (env.x[state_index]))
                insertcols!(hook.tmp, "next_state_"*replace(hook.extra_state_names[extra_state_cntr], "_i_" => "_u_") => (env.x[state_index])*hook.extra_state_paras[extra_state_cntr])
            else
                insertcols!(hook.tmp, "next_state_"*hook.extra_state_names[extra_state_cntr] => opstates[state_index,1])
                insertcols!(hook.tmp, "next_state_"*replace(hook.extra_state_names[extra_state_cntr], "_i_" => "_u_") => (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr])
            end
        extra_state_cntr+=1
        end
    end

    insertcols!(hook.tmp, :reward => env.reward)
    insertcols!(hook.tmp, :done => env.done)

    append!(hook.df, hook.tmp)
    hook.tmp = DataFrame()
    
    if isa(agent, MultiAgentGridController)

        if length(hook.reward) != length(agent.agents)
            hook.reward = zeros(length(agent.agents))
        end
        if length(hook.policy_names) != length(agent.agents)
            hook.policy_names = [s for s in keys(agent.agents)]
        end
    
        i = 1
        for name in hook.policy_names
            hook.reward[i] += reward(env, name)
            i += 1
        end
    else
        hook.reward[1] += env.reward
    end
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
    if isa(agent, MultiAgentGridController)
        hook.reward = zeros(length(agent.agents))
    else
        hook.reward[1] = 0.0
    end

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
        matrix_to_plot = reduce(hcat, hook.rewards)
        p = lineplot(matrix_to_plot[1,:], ylim=(minimum(matrix_to_plot), maximum(matrix_to_plot)), name=hook.policy_names[1], title="Total reward per episode", xlabel="Episode", ylabel="Score")
        for i in 2:length(hook.rewards[1])
            lineplot!(p, matrix_to_plot[i,:], name=hook.policy_names[i])
        end
        println(p)
    end

end