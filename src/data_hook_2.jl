Base.@kwdef mutable struct DataHook2 <: AbstractHook

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
    tmp = []
    column_names = []
    list_iterator = 0
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

    v_dq = []
    v_mag_inv = []
    v_mag_poc = []
    i_dq = []
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

function append_tmp(hook, value)
    hook.tmp[hook.list_iterator] = value
    hook.list_iterator += 1
end

function (hook::DataHook2)(::PreExperimentStage, agent, env, training = false)
    @timeit to "Pre-Experiment" begin

    # rest
    #hook.df = DataFrame()
    #hook.ep = 1

    # add states of chosen sources to the state and action plotting list
    # with this method, in addition to the states at L and C, one also obtains the states of R


    for source in hook.collect_sources
        para = env.nc.parameters["source"][source]
        indices = GetSourceStateIndices(env.nc,source)

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
        indices = GetCableStateIndices(env.nc,cable)
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
        indices = GetLoadStateIndices(env.nc,load)
        for id in indices["load$load"]["state_indices"]
            if !(env.state_ids[id] in hook.collect_state_ids)
                push!(hook.collect_state_ids,env.state_ids[id])
                if occursin("_v", env.state_ids[id])
                    if occursin("R", para["impedance"]) && occursin("C", para["impedance"])
                        push!(hook.extra_state_ids,id)
                        push!(hook.extra_state_paras,(para["R"],(para["C"])*(para["C"]+GetCSumCableNode(env.nc.num_sources+load,env.nc))^(-1),(GetCSumCableNode(env.nc.num_sources+load,env.nc))*(para["C"]+GetCSumCableNode(env.nc.num_sources+load,env.nc))^(-1)), )
                        push!(hook.extra_state_names,(replace(env.state_ids[id], "_v_C_total" => "_i_R"),replace(env.state_ids[id], "_v_C_total" => "_i_C"),replace(env.state_ids[id], "_v_C_total" => "_i_C_cables")))
                    elseif occursin("R", para["impedance"])
                        push!(hook.extra_state_ids,id)
                        push!(hook.extra_state_paras,(para["R"],0,(1)))
                        push!(hook.extra_state_names,(replace(env.state_ids[id], "_v_C_total" => "_i_R"),replace(env.state_ids[id], "_v_C_total" => "_i_C_cables")))
                    elseif occursin("C", para["impedance"])
                        push!(hook.extra_state_ids,id)
                        push!(hook.extra_state_paras,(0,(para["C"])*(para["C"]+GetCSumCableNode(env.nc.num_sources+load,env.nc))^(-1),(GetCSumCableNode(env.nc.num_sources+load,env.nc))*(para["C"]+GetCSumCableNode(env.nc.num_sources+load,env.nc))^(-1) ))
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

    hook.A,hook.B ,_ ,_ = GetSystem(env.nc)
    hook.collect_state_paras = GetStateParameters(env.nc)

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
end

function (hook::DataHook2)(::PreActStage, agent, env, action, training = false)
    @timeit to "Pre-Act" begin

    if hook.list_iterator == 0
        push!(hook.column_names, :episode)
        push!(hook.tmp, hook.ep)
        push!(hook.column_names, :time)
        push!(hook.tmp, env.t)
    else
        append_tmp(hook, hook.ep)
        append_tmp(hook, env.t)
    end

    if length(hook.policy_names) != length(agent.agents)
        hook.policy_names = [s for s in keys(agent.agents)]
    end

    #TODO: append_tmp for Septimus Code
    if findfirst(x -> x == "classic", hook.policy_names) !== nothing

        ClassicalPolicy = agent.agents["classic"]["policy"].policy

        for idx in hook.debug
            if hook.list_iterator == 0
                push!(hook.column_names, Symbol("debug_$(idx)"))
                push!(hook.tmp, ClassicalPolicy.Source.debug[idx])
            else
                append_tmp(hook, ClassicalPolicy.Source.debug[idx])
            end
        end

        for idx in hook.v_dq
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_v_d"))
                push!(hook.tmp, ClassicalPolicy.Source.V_dq0[s_idx, 1])

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_v_q"))
                push!(hook.tmp, ClassicalPolicy.Source.V_dq0[s_idx, 2])
            end
        end

        for idx in hook.i_dq
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_i_d"))
                push!(hook.tmp, ClassicalPolicy.Source.I_dq0[s_idx, 1])

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_i_q"))
                push!(hook.tmp, ClassicalPolicy.Source.I_dq0[s_idx, 2])
            end
        end

        for idx in hook.power_pq_inv
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_p_inv"))
                push!(hook.tmp, ClassicalPolicy.Source.p_q_inv[s_idx, 1])

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_q_inv"))
                push!(hook.tmp, ClassicalPolicy.Source.p_q_inv[s_idx, 2])
            end
        end

        for idx in hook.power_pq_poc
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_p_poc"))
                push!(hook.tmp, ClassicalPolicy.Source.p_q_poc[s_idx, 1])
                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_q_poc"))
                push!(hook.tmp, ClassicalPolicy.Source.p_q_poc[s_idx, 2])
            end
        end

        for idx in hook.v_mag_inv
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                v_mag = ClarkeMag((ClassicalPolicy.Source.Vdc[s_idx]/2)*ClassicalPolicy.Source.Vd_abc_new[s_idx, :, end])

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_v_mag_inv"))
                push!(hook.tmp, v_mag)

                #hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_vrms_a"))
                #push!(hook.tmp, ClassicalPolicy.Source.V_ph[s_idx, 1, 2])
            end
        end

        for idx in hook.v_mag_poc
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                v_mag = ClarkeMag(ClassicalPolicy.Source.V_filt_cap[s_idx, :, end])

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_v_mag_poc"))
                push!(hook.tmp, v_mag)

                #hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_vrms_a"))
                #push!(hook.tmp, ClassicalPolicy.Source.V_ph[s_idx, 1, 2])
            end
        end

        for idx in hook.i_mag_inv
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                i_mag = ClarkeMag(ClassicalPolicy.Source.I_filt_inv[s_idx, :, end])

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_i_mag_inv"))
                push!(hook.tmp, i_mag)

                #hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_irms_a"))
                #push!(hook.tmp, ClassicalPolicy.Source.I_ph[s_idx, 1, 2])
            end
        end

        for idx in hook.i_mag_poc
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                i_mag = ClarkeMag(ClassicalPolicy.Source.I_filt_poc[s_idx, :, end])

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_i_mag_poc"))
                push!(hook.tmp, i_mag)

                #hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_irms_a"))
                #push!(hook.tmp, ClassicalPolicy.Source.I_ph[s_idx, 1, 2])
            end
        end

        for idx in hook.freq
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                freq = ClassicalPolicy.Source.f_source[s_idx, 1, end]

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_freq"))
                push!(hook.tmp, freq)
            end
        end

        if !isempty(hook.angles)
            θ_ref = ClassicalPolicy.Source.θ_avg[1, end]
        end
        for idx in hook.angles

            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)

            if s_idx !== nothing

                θ_source = (ClassicalPolicy.Source.θ_source[s_idx, 1, end] - θ_ref
                        + ClassicalPolicy.Source.ts*π*ClassicalPolicy.Source.fsys)*180/π

                if θ_source > 180
                    θ_source = θ_source - 360
                elseif θ_source < -180
                    θ_source = θ_source + 360
                end

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_θ"))
                push!(hook.tmp, θ_source)
            end
        end

        for idx in hook.i_sat
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                i_sat = sqrt(2)*(ClassicalPolicy.Source.Vdc[s_idx]/2)*ClarkeMag(ClassicalPolicy.Source.s_dq0_avg[s_idx, :] .- ClassicalPolicy.Source.s_lim[s_idx, :])/ClassicalPolicy.Source.v_max[s_idx]

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_i_sat"))
                push!(hook.tmp, i_sat)
            end
        end

        for idx in hook.i_err
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                i_err = ClarkeMag(ClassicalPolicy.Source.I_err[s_idx, :, end])

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_i_err"))
                push!(hook.tmp, i_err)
            end
        end

        for idx in hook.i_err_t
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                i_err_t = ClarkeMag(ClassicalPolicy.Source.I_err_t[s_idx, :])

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_i_err_t"))
                push!(hook.tmp, i_err_t)
            end
        end

        for idx in hook.v_sat
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                v_sat = sqrt(2)*ClarkeMag(ClassicalPolicy.Source.I_ref_dq0[s_idx, :] .- ClassicalPolicy.Source.I_lim[s_idx, :])/(0.98*ClassicalPolicy.Source.i_max[s_idx])

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_v_sat"))
                push!(hook.tmp, v_sat)
            end
        end

        for idx in hook.v_err
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                v_err = ClarkeMag(ClassicalPolicy.Source.V_err[s_idx, :, end])

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_v_err"))
                push!(hook.tmp, v_err)
            end
        end

        for idx in hook.i_err_t
            s_idx = findfirst(x -> x == idx, ClassicalPolicy.Source_Indices)
            if s_idx !== nothing
                v_err_t = ClarkeMag(ClassicalPolicy.Source.V_err_t[s_idx, :])

                hook.list_iterator == 0 && push!(hook.column_names, Symbol("source$(idx)_v_err_t"))
                push!(hook.tmp, v_err_t)
            end
        end

    end

    for idx in hook.collect_vdc_ids
        if hook.list_iterator == 0
            push!(hook.column_names, Symbol("source$(idx)_vdc"))
            push!(hook.tmp,  env.v_dc[idx])
        else
            append_tmp(hook, env.v_dc[idx])
        end 
        
    end

    if hook.collect_reference
        #push!(hook.tmp, :reference => reference(env.t))
        refs = reference(env.t)
        for i = 1:length(refs)
            if hook.list_iterator == 0
                push!(hook.column_names, Symbol("reference_$i"))
                push!(hook.tmp, refs[i])
            else 
                append_tmp(hook, refs[i])
            end 
            
        end
    end

    # computation of voltages at Ls and currents at Cs
    @timeit to "Matrix Computation" begin
    states_x = Vector( env.x )
    opstates=(hook.A * states_x + hook.B * (Vector(env.action)) ) .* (hook.collect_state_paras)
    extra_state_cntr= 1
    end

    # saving all states given in the collect_state_ids list + possible R States
    @timeit to "Collect States" begin
    for state_id in hook.collect_state_ids
        state_index = findfirst(x -> x == state_id, env.state_ids)

        @timeit to "Normal" begin
        if hook.list_iterator == 0
            @timeit to "Normal firstrun" begin
            push!(hook.column_names, Symbol(state_id))
            push!(hook.tmp, (env.x[state_index]))
            push!(hook.column_names, Symbol(replace(state_id, "_i_" => "_v_", "_v_" => "_i_")))
            push!(hook.tmp, opstates[state_index,1])
            end
        else
            @timeit to "Normal fast" begin
            append_tmp(hook, (env.x[state_index]))
            append_tmp(hook, opstates[state_index,1])
            end
        end
        end
        

        @timeit to "Extra States" begin
        if state_index in hook.extra_state_ids
            if occursin("load",state_id)
                if hook.extra_state_paras[extra_state_cntr][1] != 0

                    if hook.list_iterator == 0
                        push!(hook.column_names, Symbol(hook.extra_state_names[extra_state_cntr][1]))
                        push!(hook.tmp, (env.x[state_index])*hook.extra_state_paras[extra_state_cntr][1]^(-1))
                        push!(hook.column_names, Symbol(replace(hook.extra_state_names[extra_state_cntr][1], "_i_" => "_v_")))
                        push!(hook.tmp, (env.x[state_index]))
                    else
                        append_tmp(hook, (env.x[state_index])*hook.extra_state_paras[extra_state_cntr][1]^(-1))
                        append_tmp(hook, (env.x[state_index]))
                    end
                    
                end
                if  hook.extra_state_paras[extra_state_cntr][2] != 0

                    if hook.list_iterator == 0 
                        push!(hook.column_names, Symbol(hook.extra_state_names[extra_state_cntr][2]))
                        push!(hook.tmp, (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][2])
                        push!(hook.column_names, Symbol(replace(hook.extra_state_names[extra_state_cntr][2], "_i_" => "_v_")))
                        push!(hook.tmp, (env.x[state_index]))
                    else
                        append_tmp(hook, (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][2])
                        append_tmp(hook, (env.x[state_index]))
                    end 
                    
                end
                if  hook.extra_state_paras[extra_state_cntr][3] != 0

                    if hook.list_iterator == 0
                        push!(hook.column_names, Symbol(hook.extra_state_names[extra_state_cntr][3]))
                        push!(hook.tmp, (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][3])
                        push!(hook.column_names, Symbol(replace(hook.extra_state_names[extra_state_cntr][3], "_i_" => "_v_")))
                        push!(hook.tmp, (env.x[state_index]))
                    else
                        append_tmp(hook, (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][3])
                        append_tmp(hook, (env.x[state_index]))
                    end 

                end
            elseif occursin("_i_", state_id)
                
                if hook.list_iterator == 0 
                    push!(hook.column_names, Symbol(hook.extra_state_names[extra_state_cntr]))
                    push!(hook.tmp, (env.x[state_index]))
                    push!(hook.column_names, Symbol(replace(hook.extra_state_names[extra_state_cntr], "_i_" => "_v_")))
                    push!(hook.tmp, (env.x[state_index])*hook.extra_state_paras[extra_state_cntr])
                else
                    append_tmp(hook, (env.x[state_index]))
                    append_tmp(hook, (env.x[state_index])*hook.extra_state_paras[extra_state_cntr])
                end 

            else
                
                if hook.list_iterator == 0
                    push!(hook.column_names, Symbol(hook.extra_state_names[extra_state_cntr]))
                    push!(hook.tmp, opstates[state_index,1])
                    push!(hook.column_names, Symbol(replace(hook.extra_state_names[extra_state_cntr], "_i_" => "_v_")))
                    push!(hook.tmp, (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr])
                else 
                    append_tmp(hook, opstates[state_index,1])
                    append_tmp(hook, (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr])
                end 
                
            end
            extra_state_cntr+=1
        end
        end
    end
    end



    end
end

function (hook::DataHook2)(::PostActStage, agent, env, training = false)
    @timeit to "Post-Act" begin

    states_x = Vector( env.x )
    opstates = (hook.A * states_x + hook.B * (Vector(env.action)) ) .* (hook.collect_state_paras)
    extra_state_cntr= 1
    for state_id in hook.collect_state_ids
        state_index = findfirst(x -> x == state_id, env.state_ids)

        if hook.list_iterator == 0 
            push!(hook.column_names, Symbol("next_state_"*state_id))
            push!(hook.tmp, (env.x[state_index]))
            push!(hook.column_names, Symbol("next_state_"*replace(state_id, "_i_" => "_v_", "_v_" => "_i_")))
            push!(hook.tmp, opstates[state_index,1])
        else 
            append_tmp(hook, (env.x[state_index]))
            append_tmp(hook, opstates[state_index,1])
        end 

        if state_index in hook.extra_state_ids
            if occursin("load",state_id)
                if hook.extra_state_paras[extra_state_cntr][1] != 0

                    if hook.list_iterator == 0 
                        push!(hook.column_names, Symbol("next_state_"*hook.extra_state_names[extra_state_cntr][1]))
                        push!(hook.tmp, (env.x[state_index])*hook.extra_state_paras[extra_state_cntr][1]^(-1))
                        push!(hook.column_names, Symbol("next_state_"*replace(hook.extra_state_names[extra_state_cntr][1], "_i_" => "_v_")))
                        push!(hook.tmp, (env.x[state_index]))
                    else 
                        append_tmp(hook, (env.x[state_index])*hook.extra_state_paras[extra_state_cntr][1]^(-1))
                        append_tmp(hook, (env.x[state_index]))
                    end 

                end
                if  hook.extra_state_paras[extra_state_cntr][2] != 0

                    if hook.list_iterator == 0 
                        push!(hook.column_names, Symbol("next_state_"*hook.extra_state_names[extra_state_cntr][2]))
                        push!(hook.tmp, (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][2])
                        push!(hook.column_names, Symbol("next_state_"*replace(hook.extra_state_names[extra_state_cntr][2], "_i_" => "_v_")))
                        push!(hook.tmp, (env.x[state_index]))
                    else 
                        append_tmp(hook, (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][2])
                        append_tmp(hook, (env.x[state_index]))
                    end 

                end
                if  hook.extra_state_paras[extra_state_cntr][3] != 0

                    if hook.list_iterator == 0 
                        push!(hook.column_names, Symbol("next_state_"*hook.extra_state_names[extra_state_cntr][3]))
                        push!(hook.tmp, (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][3])
                        push!(hook.column_names, Symbol("next_state_"*replace(hook.extra_state_names[extra_state_cntr][3], "_i_" => "_v_")))
                        push!(hook.tmp, (env.x[state_index]))
                    else 
                        append_tmp(hook, (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr][3])
                        append_tmp(hook, (env.x[state_index]))
                    end 

                end
            elseif occursin("_i_", state_id)
                if hook.list_iterator == 0 
                    push!(hook.column_names, Symbol("next_state_"*hook.extra_state_names[extra_state_cntr]))
                    push!(hook.tmp, (env.x[state_index]))
                    push!(hook.column_names, Symbol("next_state_"*replace(hook.extra_state_names[extra_state_cntr], "_i_" => "_v_")))
                    push!(hook.tmp, (env.x[state_index])*hook.extra_state_paras[extra_state_cntr])
                else 
                    append_tmp(hook, (env.x[state_index]))
                    append_tmp(hook, (env.x[state_index])*hook.extra_state_paras[extra_state_cntr])
                end 

            else
                if hook.list_iterator == 0 
                    push!(hook.column_names, Symbol("next_state_"*hook.extra_state_names[extra_state_cntr]))
                    push!(hook.tmp, opstates[state_index,1])
                    push!(hook.column_names, Symbol("next_state_"*replace(hook.extra_state_names[extra_state_cntr], "_i_" => "_v_")))
                    push!(hook.tmp, (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr])
                else 
                    append_tmp(hook, opstates[state_index,1])
                    append_tmp(hook, (opstates[state_index,1])*hook.extra_state_paras[extra_state_cntr])
                end 

            end
        extra_state_cntr+=1
        end
    end

    for action_id in hook.collect_action_ids
        action_index = findfirst(x -> x == action_id, env.action_ids)

        if hook.list_iterator == 0 
            push!(hook.column_names, Symbol(action_id))
            push!(hook.tmp, (env.action[action_index]))
        else 
            append_tmp(hook, (env.action[action_index]))
        end 
    end

    if hook.list_iterator == 0 
        push!(hook.column_names, :reward)
        push!(hook.tmp, env.reward)
        push!(hook.column_names, :done)
        push!(hook.tmp, env.done)
    else 
        append_tmp(hook, env.reward)
        append_tmp(hook, env.done)
    end 
    
    
    if hook.list_iterator == 0
        hook.df = DataFrame([[i] for i in hook.tmp], hook.column_names; copycols=false)
        hook.list_iterator = 1
    else
        push!(hook.df, hook.tmp)
        hook.list_iterator = 1
    end

    if training
        if isa(agent, MultiController)

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
end

function (hook::DataHook2)(::PostEpisodeStage, agent, env, training = false)
    @timeit to "Post-Episode" begin

    hook.ep += 1

    if training
        if isa(agent, MultiController)
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

        if isa(agent, MultiController)
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
end

function (hook::DataHook2)(::PostExperimentStage, agent, env, training = false)
    @timeit to "Post-Experiment" begin

    if hook.save_data_to_hd
        isdir(hook.dir) || mkdir(hook.dir)
        Arrow.write(hook.dir * "data.arrow", hook.df)
    end

    if hook.plot_rewards
        if hook.is_inner_hook_RL
            if training
                matrix_to_plot = reduce(hcat, hook.rewards)

                p = lineplot(matrix_to_plot[1,:], ylim=(minimum(matrix_to_plot), maximum(matrix_to_plot)), name=hook.policy_names[1], title="Total reward per episode", xlabel="Episode", ylabel="Score")

                for i in 2:length(hook.rewards[1])
                    lineplot!(p, matrix_to_plot[i,:], name=hook.policy_names[i])
                end

                display(p)
            end
        else
            matrix_to_plot = reduce(hcat, hook.rewards)

            p = lineplot(matrix_to_plot[1,:], ylim=(minimum(matrix_to_plot), maximum(matrix_to_plot)), title="Total reward per episode", xlabel="Episode", ylabel="Score")

            for i in 2:length(hook.rewards[1])
                lineplot!(p, matrix_to_plot[i,:], name=hook.policy_names[i])
            end

            display(p)
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
end
