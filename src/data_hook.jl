using DataFrames
using UnicodePlots


Base.@kwdef mutable struct data_hook <: AbstractHook

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

    is_inner_hook_RL = false
    bestNNA = nothing
    bestreward = -1000000.0
    bestepisode = 1
    currentNNA = nothing

    collect_reference = false
    collect_vdc_ids = []

    vdq = []
    v_mag_inv = []
    v_mag_cap = []
    idq = []
    i_mag_inv = []
    i_mag_poc = []
    power_pq_inv = []
    power_pq_poc = []
    freq = []
    angles = []
    i_sat = []
    i_err = []
    i_err_t = []
    v_sat = []
    v_err = []
    v_err_t = []

    debug = []

end

function (hook::data_hook)(::PreExperimentStage, agent, env, training = false)

    # rest
    #hook.df = DataFrame()
    #hook.ep = 1

    # add states of chosen sources to the state and action plotting list
    # with this method, in addition to the states at L and C, one also obtains the states of R


    for source in hook.collect_sources
        para = env.nc.parameters["source"][source]
        indices = get_source_state_indices(env.nc,source)

        for id in indices["source$source"]["state_indices"]

            if !(env.state_ids[id] in hook.collect_state_ids)
                push!(hook.collect_state_ids,env.state_ids[id])
            end

            if occursin("_L1", env.state_ids[id])
                push!(hook.extra_state_ids,id)
                push!(hook.extra_state_paras,para["R1"])
                push!(hook.extra_state_names,replace(env.state_ids[id], "_L1" => "_R1"))
            elseif occursin("v_C", env.state_ids[id])
                if haskey(para, "R_C")
                    push!(hook.extra_state_ids,id)
                    push!(hook.extra_state_paras,para["R_C"])
                    push!(hook.extra_state_names, replace(env.state_ids[id], "v_C" => "i_R_C"))
                end
            elseif occursin("_L2", env.state_ids[id])
                push!(hook.extra_state_ids,id)
                push!(hook.extra_state_paras,para["R2"])
                push!(hook.extra_state_names, replace(env.state_ids[id], "_L2" => "_R2"))
            end
        end

        for id in indices["source$source"]["action_indices"]
            if !(env.action_ids[id] in hook.collect_action_ids)
                push!(hook.collect_action_ids,env.action_ids[id])
            end
        end
    end

    # add states of chosen cables to the state plotting list
    # with this method, in addition to the state at L, one also obtains the states of R
    for cable in hook.collect_cables
        para = env.nc.parameters["cable"][cable]
        indices = get_cable_state_indices(env.nc,cable)
        for id in indices["cable$cable"]["state_indices"]
            if !(env.state_ids[id] in hook.collect_state_ids)
                push!(hook.collect_state_ids,env.state_ids[id])
                push!(hook.extra_state_ids,id)
                push!(hook.extra_state_paras,para["R"])
                push!(hook.extra_state_names, replace(env.state_ids[id], "_L" => "_R"))
            end
        end
    end

    # add states of chosen loads to the state plotting list
    # with this method, in addition to the states at L and C, one also obtains the states of R
    for load in hook.collect_loads
        para = env.nc.parameters["load"][load]
        indices = get_load_state_indices(env.nc,load)
        for id in indices["load$load"]["state_indices"]
            if !(env.state_ids[id] in hook.collect_state_ids)
                push!(hook.collect_state_ids,env.state_ids[id])
                if occursin("_v", env.state_ids[id])
                    if occursin("R", para["impedance"]) && occursin("C", para["impedance"])
                        push!(hook.extra_state_ids,id)
                        push!(hook.extra_state_paras,(para["R"],(para["C"])*(para["C"]+get_C_sum_cable_node(env.nc.num_sources+load,env.nc))^(-1),(get_C_sum_cable_node(env.nc.num_sources+load,env.nc))*(para["C"]+get_C_sum_cable_node(env.nc.num_sources+load,env.nc))^(-1)), )
                        push!(hook.extra_state_names,(replace(env.state_ids[id], "_v_C_total" => "_i_R"),replace(env.state_ids[id], "_v_C_total" => "_i_C"),replace(env.state_ids[id], "_v_C_total" => "_i_C_cables")))
                    elseif occursin("R", para["impedance"])
                        push!(hook.extra_state_ids,id)
                        push!(hook.extra_state_paras,(para["R"],0,(1)))
                        push!(hook.extra_state_names,(replace(env.state_ids[id], "_v_C_total" => "_i_R"),replace(env.state_ids[id], "_v_C_total" => "_i_C_cables")))
                    elseif occursin("C", para["impedance"])
                        push!(hook.extra_state_ids,id)
                        push!(hook.extra_state_paras,(0,(para["C"])*(para["C"]+get_C_sum_cable_node(env.nc.num_sources+load,env.nc))^(-1),(get_C_sum_cable_node(env.nc.num_sources+load,env.nc))*(para["C"]+get_C_sum_cable_node(env.nc.num_sources+load,env.nc))^(-1) ))
                        push!(hook.extra_state_names,(replace(env.state_ids[id], "_v_C_total" => "_i_C"),replace(env.state_ids[id], "_v_C_total" => "_i_C_cables")))
                    else
                        push!(hook.extra_state_ids,id)
                        push!(hook.extra_state_paras,(0,0,1))
                        push!(hook.extra_state_names,replace(env.state_ids[id], "_v_C_total" => "_i_C_cables"))
                    end
                end
            end
        end
    end

    hook.A,hook.B ,_ ,_ = get_sys(env.nc)
    hook.collect_state_paras = get_state_paras(env.nc)

    if hook.is_inner_hook_RL

        if length(hook.policy_names) != length(agent.agents)
            hook.policy_names = [s for s in keys(agent.agents)]
        end

        if isnothing(hook.currentNNA)
            hook.currentNNA = Dict()

            for name in hook.policy_names
                if isa(agent.agents[name]["policy"], Agent)
                    hook.currentNNA[name] = deepcopy(agent.agents[name]["policy"].policy.policy.behavior_actor)
                end
            end
        end

        if isnothing(hook.bestNNA)
            hook.bestNNA = Dict()

            for name in hook.policy_names
                if isa(agent.agents[name]["policy"], Agent)
                    hook.bestNNA[name] = deepcopy(agent.agents[name]["policy"].policy.policy.behavior_actor)
                end
            end
        end

        if training
            for name in hook.policy_names
                if isa(agent.agents[name]["policy"], Agent)
                    agent.agents[name]["policy"].policy.policy.behavior_actor = deepcopy(hook.currentNNA[name])
                end
            end
        else
            for name in hook.policy_names
                if isa(agent.agents[name]["policy"], Agent)
                    agent.agents[name]["policy"].policy.policy.behavior_actor = deepcopy(hook.bestNNA[name])
                end
            end
        end
    end
end

function (hook::data_hook)(::PreActStage, agent, env, action, training = false)

    insertcols!(hook.tmp, :episode => hook.ep)
    insertcols!(hook.tmp, :time => Float32(env.t))

    if length(hook.policy_names) != length(agent.agents)
        hook.policy_names = [s for s in keys(agent.agents)]
    end

    if findfirst(x -> x == "classic", hook.policy_names) !== nothing

        Classical_Policy = agent.agents["classic"]["policy"].policy

        for idx in hook.debug

            insertcols!(hook.tmp, "debug_$(idx)" => Classical_Policy.Source.debug[idx])
        end

        for idx in hook.vdq
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                insertcols!(hook.tmp, "source$(idx)_vd" => Classical_Policy.Source.V_dq0[s_idx, 1])
                insertcols!(hook.tmp, "source$(idx)_vq" => Classical_Policy.Source.V_dq0[s_idx, 2])
            end
        end

        for idx in hook.idq
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                insertcols!(hook.tmp, "source$(idx)_id" => Classical_Policy.Source.I_dq0[s_idx, 1])
                insertcols!(hook.tmp, "source$(idx)_iq" => Classical_Policy.Source.I_dq0[s_idx, 2])
            end
        end

        for idx in hook.power_pq_inv
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                insertcols!(hook.tmp, "source$(idx)_p_inv" => Classical_Policy.Source.p_q_inv[s_idx, 1])
                insertcols!(hook.tmp, "source$(idx)_q_inv" => Classical_Policy.Source.p_q_inv[s_idx, 2])
            end
        end

        for idx in hook.power_pq_poc
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                insertcols!(hook.tmp, "source$(idx)_p_poc" => Classical_Policy.Source.p_q_poc[s_idx, 1])
                insertcols!(hook.tmp, "source$(idx)_q_poc" => Classical_Policy.Source.p_q_poc[s_idx, 2])
            end
        end

        for idx in hook.v_mag_inv
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                v_mag = Clarke_mag((Classical_Policy.Source.Vdc[s_idx]/2)*Classical_Policy.Source.Vd_abc_new[s_idx, :, end])
                insertcols!(hook.tmp, "source$(idx)_v_mag_inv" => v_mag)
                #insertcols!(hook.tmp, "source$(idx)_vrms_a" => Classical_Policy.Source.V_ph[s_idx, 1, 2])
            end
        end

        for idx in hook.v_mag_cap
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                v_mag = Clarke_mag(Classical_Policy.Source.V_filt_cap[s_idx, :, end])
                insertcols!(hook.tmp, "source$(idx)_v_mag_cap" => v_mag)
                #insertcols!(hook.tmp, "source$(idx)_vrms_a" => Classical_Policy.Source.V_ph[s_idx, 1, 2])
            end
        end

        for idx in hook.i_mag_inv
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                i_mag = Clarke_mag(Classical_Policy.Source.I_filt_inv[s_idx, :, end])
                insertcols!(hook.tmp, "source$(idx)_i_mag_inv" => i_mag)
                #insertcols!(hook.tmp, "source$(idx)_irms_a" => Classical_Policy.Source.I_ph[s_idx, 1, 2])
            end
        end

        for idx in hook.i_mag_poc
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                i_mag = Clarke_mag(Classical_Policy.Source.I_filt_poc[s_idx, :, end])
                insertcols!(hook.tmp, "source$(idx)_i_mag_poc" => i_mag)
                #insertcols!(hook.tmp, "source$(idx)_irms_a" => Classical_Policy.Source.I_ph[s_idx, 1, 2])
            end
        end

        for idx in hook.freq
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                freq = Classical_Policy.Source.f_source[s_idx, 1, end]
                insertcols!(hook.tmp, "source$(idx)_freq" => freq)
            end
        end

        if !isempty(hook.angles)
            θ_ref = Classical_Policy.Source.θ_avg[1, end]
        end
        for idx in hook.angles

            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)

            if s_idx !== nothing

                θ_source = (Classical_Policy.Source.θ_source[s_idx, 1, end] - θ_ref
                        + Classical_Policy.Source.ts*π*Classical_Policy.Source.fsys)*180/π

                if θ_source > 180
                    θ_source = θ_source - 360
                elseif θ_source < -180
                    θ_source = θ_source + 360
                end
                insertcols!(hook.tmp, "source$(idx)_θ" => θ_source)
            end
        end

        for idx in hook.i_sat
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                i_sat = sqrt(2)*(Classical_Policy.Source.Vdc[s_idx]/2)*Clarke_mag(Classical_Policy.Source.s_dq0_avg[s_idx, :] .- Classical_Policy.Source.s_lim[s_idx, :])/Classical_Policy.Source.v_max[s_idx]
                insertcols!(hook.tmp, "source$(idx)_i_sat" => i_sat)
            end
        end

        for idx in hook.i_err
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                i_err = Clarke_mag(Classical_Policy.Source.I_err[s_idx, :, end])
                insertcols!(hook.tmp, "source$(idx)_i_err" => i_err)
            end
        end

        for idx in hook.i_err_t
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                i_err_t = Clarke_mag(Classical_Policy.Source.I_err_t[s_idx, :])
                insertcols!(hook.tmp, "source$(idx)_i_err_t" => i_err_t)
            end
        end

        for idx in hook.v_sat
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                v_sat = sqrt(2)*Clarke_mag(Classical_Policy.Source.I_ref_dq0[s_idx, :] .- Classical_Policy.Source.I_lim[s_idx, :])/(0.98*Classical_Policy.Source.i_max[s_idx])
                insertcols!(hook.tmp, "source$(idx)_v_sat" => v_sat)
            end
        end

        for idx in hook.v_err
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                v_err = Clarke_mag(Classical_Policy.Source.V_err[s_idx, :, end])
                insertcols!(hook.tmp, "source$(idx)_v_err" => v_err)
            end
        end

        for idx in hook.i_err_t
            s_idx = findfirst(x -> x == idx, Classical_Policy.Source_Indices)
            if s_idx !== nothing
                v_err_t = Clarke_mag(Classical_Policy.Source.V_err_t[s_idx, :])
                insertcols!(hook.tmp, "source$(idx)_v_err_t" => v_err_t)
            end
        end

    end

    for idx in hook.collect_vdc_ids
        insertcols!(hook.tmp, "source$(idx)_vdc" => env.v_dc[idx])
    end

    if hook.collect_reference
        #insertcols!(hook.tmp, :reference => reference(env.t))
        refs = reference(env.t)
        for i = 1:length(refs)
            insertcols!(hook.tmp, "reference_$i" => refs[i])
        end
    end

    # computation of voltages at Ls and currents at Cs
    states_x = Vector( env.x )
    opstates=(hook.A * states_x + hook.B * (Vector(env.action)) ) .* (hook.collect_state_paras)
    extra_state_cntr= 1

    # saving all states given in the collect_state_ids list + possible R States
    for state_id in hook.collect_state_ids
        state_index = findfirst(x -> x == state_id, env.state_ids)

        insertcols!(hook.tmp, state_id => (env.x[state_index]))
        insertcols!(hook.tmp, replace(state_id, "_i_" => "_v_", "_v_" => "_i_") => opstates[state_index,1])

        if state_index in hook.extra_state_ids
            if occursin("load",state_id)
                if hook.extra_state_paras[extra_state_cntr][1] != 0
                    insertcols!(hook.tmp, hook.extra_state_names[extra_state_cntr][1] => (env.x[state_index])*hook.extra_state_paras[extra_state_cntr][1]^(-1))
                    insertcols!(hook.tmp, replace(hook.extra_state_names[extra_state_cntr][1], "_i_" => "_v_") => (env.x[state_index]))
                end
                if  hook.extra_state_paras[extra_state_cntr][2] != 0
                    insertcols!(hook.tmp, hook.extra_state_names[extra_state_cntr][2] => (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][2])
                    insertcols!(hook.tmp, replace(hook.extra_state_names[extra_state_cntr][2], "_i_" => "_v_") => (env.x[state_index]))
                end
                if  hook.extra_state_paras[extra_state_cntr][3] != 0
                    insertcols!(hook.tmp, hook.extra_state_names[extra_state_cntr][3] => (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][3])
                    insertcols!(hook.tmp, replace(hook.extra_state_names[extra_state_cntr][3], "_i_" => "_v_") => (env.x[state_index]))
                end
            elseif occursin("_i_", state_id)
                insertcols!(hook.tmp, hook.extra_state_names[extra_state_cntr] => (env.x[state_index]))
                insertcols!(hook.tmp, replace(hook.extra_state_names[extra_state_cntr], "_i_" => "_v_") => (env.x[state_index])*hook.extra_state_paras[extra_state_cntr])
            else
                insertcols!(hook.tmp, hook.extra_state_names[extra_state_cntr] => opstates[state_index,1])
                insertcols!(hook.tmp, replace(hook.extra_state_names[extra_state_cntr], "_i_" => "_v_") => (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr])
            end
        extra_state_cntr+=1
        end
    end



    insertcols!(hook.tmp, :action => Ref(action))

end

function (hook::data_hook)(::PostActStage, agent, env, training = false)

    states_x = Vector( env.x )
    opstates = (hook.A * states_x + hook.B * (Vector(env.action)) ) .* (hook.collect_state_paras)
    extra_state_cntr= 1
    for state_id in hook.collect_state_ids
        state_index = findfirst(x -> x == state_id, env.state_ids)

        insertcols!(hook.tmp, "next_state_"*state_id => (env.x[state_index]))
        insertcols!(hook.tmp, "next_state_"*replace(state_id, "_i_" => "_v_", "_v_" => "_i_") => opstates[state_index,1])

        if state_index in hook.extra_state_ids
            if occursin("load",state_id)
                if hook.extra_state_paras[extra_state_cntr][1] != 0
                    insertcols!(hook.tmp, "next_state_"*hook.extra_state_names[extra_state_cntr][1] => (env.x[state_index])*hook.extra_state_paras[extra_state_cntr][1]^(-1))
                    insertcols!(hook.tmp, "next_state_"*replace(hook.extra_state_names[extra_state_cntr][1], "_i_" => "_v_") => (env.x[state_index]))
                end
                if  hook.extra_state_paras[extra_state_cntr][2] != 0
                    insertcols!(hook.tmp, "next_state_"*hook.extra_state_names[extra_state_cntr][2] => (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][2])
                    insertcols!(hook.tmp, "next_state_"*replace(hook.extra_state_names[extra_state_cntr][2], "_i_" => "_v_") => (env.x[state_index]))
                end
                if  hook.extra_state_paras[extra_state_cntr][3] != 0
                    insertcols!(hook.tmp, "next_state_"*hook.extra_state_names[extra_state_cntr][3] => (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][3])
                    insertcols!(hook.tmp, "next_state_"*replace(hook.extra_state_names[extra_state_cntr][3], "_i_" => "_v_") => (env.x[state_index]))
                end
            elseif occursin("_i_", state_id)
                insertcols!(hook.tmp, "next_state_"*hook.extra_state_names[extra_state_cntr] => (env.x[state_index]))
                insertcols!(hook.tmp, "next_state_"*replace(hook.extra_state_names[extra_state_cntr], "_i_" => "_v_") => (env.x[state_index])*hook.extra_state_paras[extra_state_cntr])
            else
                insertcols!(hook.tmp, "next_state_"*hook.extra_state_names[extra_state_cntr] => opstates[state_index,1])
                insertcols!(hook.tmp, "next_state_"*replace(hook.extra_state_names[extra_state_cntr], "_i_" => "_v_") => (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr])
            end
        extra_state_cntr+=1
        end
    end

    for action_id in hook.collect_action_ids
        action_index = findfirst(x -> x == action_id, env.action_ids)

        insertcols!(hook.tmp, action_id => (env.action[action_index]))
    end

    insertcols!(hook.tmp, :reward => env.reward)
    insertcols!(hook.tmp, :done => env.done)

    append!(hook.df, hook.tmp)
    hook.tmp = DataFrame()

    if training
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
        elseif isa(agent.policy, NamedPolicy)
            hook.reward[1] += reward(env, name)
        else
            hook.reward[1] += env.reward
        end
    end
end

function (hook::data_hook)(::PostEpisodeStage, agent, env, training = false)
    hook.ep += 1

    if training
        if isa(agent, MultiAgentGridController)
            if length(hook.rewards) >= 1 && sum(hook.reward) > maximum(sum.(hook.rewards))
                if hook.is_inner_hook_RL
                    for name in hook.policy_names
                        if isa(agent.agents[name]["policy"], Agent)
                            hook.bestNNA[name] = deepcopy(agent.agents[name]["policy"].policy.policy.behavior_actor)
                        end
                    end
                end
                hook.bestepisode = hook.ep
                hook.bestreward = sum(hook.reward)
            end
        else
            if length(hook.rewards) >= 1 && hook.reward > maximum(hook.rewards)
                hook.bestepisode = hook.ep
                hook.bestreward = hook.reward
            end
        end

        push!(hook.rewards, hook.reward)

        if isa(agent, MultiAgentGridController)
            hook.reward = zeros(length(agent.agents))
        else
            hook.reward = [0.0]
        end

        if hook.is_inner_hook_RL
            for name in hook.policy_names
                if isa(agent.agents[name]["policy"], Agent)
                    hook.currentNNA[name] = deepcopy(agent.agents[name]["policy"].policy.policy.behavior_actor)
                end
            end
        end
    end
end

function (hook::data_hook)(::PostExperimentStage, agent, env, training = false)

    if hook.save_data_to_hd
        isdir(hook.dir) || mkdir(hook.dir)
        Arrow.write(hook.dir * "data.arrow", hook.df)
    end

    if hook.plot_rewards
        matrix_to_plot = reduce(hcat, hook.rewards)
        if length(hook.policy_names) > 0
            p = lineplot(matrix_to_plot[1,:], ylim=(minimum(matrix_to_plot), maximum(matrix_to_plot)), name=hook.policy_names[1], title="Total reward per episode", xlabel="Episode", ylabel="Score")
        else
            p = lineplot(matrix_to_plot[1,:], ylim=(minimum(matrix_to_plot), maximum(matrix_to_plot)), title="Total reward per episode", xlabel="Episode", ylabel="Score")
        end
        for i in 2:length(hook.rewards[1])
            lineplot!(p, matrix_to_plot[i,:], name=hook.policy_names[i])
        end
        # println(p)
    end

    if hook.is_inner_hook_RL
        for name in hook.policy_names
            if isa(agent.agents[name]["policy"], Agent)
                agent.agents[name]["policy"].policy.policy.behavior_actor = deepcopy(hook.currentNNA[name])
            end
        end
    end

end
