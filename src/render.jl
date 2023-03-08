using DataFrames
using PlotlyJS
using Graphs, GraphPlot

function plot_rewards_3d(hook)
    layout = Layout(
            plot_bgcolor="#f1f3f7",
            title = "Reward over Episodes",
            scene = attr(
                xaxis_title = "Time in Seconds",
                yaxis_title = "Episodes",
                zaxis_title = "Reward"),
            width = 1000,
            height = 650,
            margin=attr(l=10, r=10, b=10, t=60, pad=10)
        )

    p = plot(scatter3d(hook.df, x = :time, y = :episode, z = :reward,
                        marker=attr(size=2, color=:reward, colorscale=[[0, "rgb(255,0,0)"], [1, "rgb(0,255,0)"]]),
                        mode = "markers"),
            config = PlotConfig(scrollZoom=true),
            layout)
    display(p)
end

function plot_best_results(;agent, env, hook, states_to_plot = nothing, actions_to_plot = nothing, plot_reward = true, plot_reference = false, use_best = true, use_noise = 0.0)
    RLBase.reset!(env)

    if isnothing(states_to_plot)
        states_to_plot = hook.collect_state_ids
    end

    if isnothing(actions_to_plot)
        actions_to_plot = hook.collect_action_ids
    end

    if isa(agent, MultiController)
        act_noise_old = Dict()
        for (name, agent_temp) in agent.agents
            if isa(agent_temp["policy"], Agent)
                act_noise_old[name] = agent_temp["policy"].policy.policy.act_noise
                agent_temp["policy"].policy.policy.act_noise = use_noise
                if use_best
                    copyto!(agent_temp["policy"].policy.policy.behavior_actor, hook.bestNNA[name])
                end
            end
        end
    else
        act_noise_old = agent.policy.act_noise
        agent.policy.act_noise = use_noise
        if use_best
            copyto!(agent.policy.behavior_actor, hook.bestNNA)
        end
    end
    
    temphook = DataHook(collect_state_ids = states_to_plot, collect_action_ids = actions_to_plot, collect_reference = true)
    
    if isa(agent, MultiController)
        ma2 = deepcopy(agent)
        for (name, policy) in ma2.agents
            if isa(policy["policy"], Agent)
                policy["policy"] = policy["policy"].policy
            end
        end
        run(ma2, env, StopAfterEpisode(1), temphook)
    else
        run(agent.policy, env, StopAfterEpisode(1), temphook)
    end

    episode_string = use_best ? string(hook.bestepisode) : string(hook.ep - 1)
    
    layout = Layout(
            plot_bgcolor="#f1f3f7",
            title = "Results<br><sub>Run with Behavior-Actor-NNA from Episode " * episode_string * "</sub>",
            xaxis_title = "Time in Seconds",
            yaxis_title = "State values",
            yaxis2 = attr(
                title="Reward",
                overlaying="y",
                side="right",
                titlefont_color="orange",
                range=[-1, 1]
            ),
            width = 1200,
            height = 850,
            margin=attr(l=100, r=200, b=80, t=100, pad=10)
        )


    traces = []

    for state_id in states_to_plot
        push!(traces, scatter(temphook.df, x = :time, y = Symbol(state_id), mode="lines", name = state_id))
    end

    for action_id in actions_to_plot
        push!(traces, scatter(temphook.df, x = :time, y = Symbol(action_id), mode="lines", name = action_id))
    end

    if plot_reward
        push!(traces, scatter(temphook.df, x = :time, y = :reward, yaxis = "y2", mode="lines", name = "Reward"))
    end

    if plot_reference
        push!(traces, scatter(temphook.df, x = :time, y = :reference_1, mode="lines", name = "Reference"))
    end

    traces = Array{GenericTrace}(traces)
    
    p = plot(traces, layout, config = PlotConfig(scrollZoom=true))
    display(p)
    
    if isa(agent, MultiController)
        for (name, agent_temp) in agent.agents
            if isa(agent_temp["policy"], Agent)
                agent_temp["policy"].policy.policy.act_noise = act_noise_old[name]
                if use_best
                    copyto!(agent_temp["policy"].policy.policy.behavior_actor, hook.currentNNA[name])
                end
            end
        end
    else
        if use_best
            copyto!(agent.policy.behavior_actor, hook.currentNNA)
        end
        agent.policy.act_noise = act_noise_old
    end
    
    RLBase.reset!(env)

    return nothing
end

"""
    RenderHookResults(; hook, 
                        states_to_plot = nothing, 
                        actions_to_plot = nothing ,
                        episode = 1, 
                        power_p_inv = [], 
                        power_q_inv = [], 
                        power_p_poc = [], 
                        power_q_poc = [], 
                        v_mag_inv   = [], 
                        v_mag_cap   = [], 
                        i_mag_inv   = [], 
                        i_mag_poc   = [],
                        freq        = [], 
                        angles      = [], 
                        i_sat       = [], 
                        i_err       = [], 
                        i_err_t     = [], 
                        v_sat       = [], 
                        v_err       = [], 
                        v_err_t     = [],
                        plot_reward = false, 
                        plot_reference = false,
                        vdc_to_plot = [],)

# Keyword Arguments
- `episode::Int`: episode to plot
- `power_p_inv::Vector{Int}`: instantaneous real power at the terminals of the source [W]
- `power_q_inv::Vector{Int}`: instantaneousimaginary power at the terminals of the source [VAi]
- `power_p_poc::Vector{Int}`: instantaneousreal power at the point of connection of the source [W]
- `power_q_poc::Vector{Int}`: instantaneous imaginary power at the point of connection of the source [VAi]
- `v_mag_inv::Vector{Int}`: scaled L₂ norm in αβγ coordinates of the 3 phase voltages at the terminals of the source [V] 
- `v_mag_cap::Vector{Int}`: scaled L₂ norm in αβγ coordinates of the 3 phase voltages at the filter capacitor of the source [V] 
- `i_mag_inv::Vector{Int}`: scaled L₂ norm in αβγ coordinates of the 3 phase currents at the terminals of the source [A] 
- `i_mag_poc::Vector{Int}`: scaled L₂ norm in αβγ coordinates of the 3 phase currents flowing into the external network of the source [A]
- `freq::Vector{Int}`: angular velocity of the 3 phase voltages over the capacior nearest to the source [Hz]
- `angles::Vector{Int}`: relative positive phase sequence angles of the 3 phase voltages at the control nodes (measured with respect to the average angle over all the grid-forming sources) [degrees]
- `i_sat::Vector{Int}`: normalised and scaled L₂ norm in αβγ coordinates of the degree of nonlinearity (or anti-windup) of the current controller [V]
- `i_err::Vector{Int}`: L₂ norm in αβγ coordinates of the current error signal (measured current subtracted from reference current in DQ0 coordinates) [A]
- `i_err_t::Vector{Int}`: L₂ norm in αβγ coordinates of the integrated current error signal [As]
- `v_sat::Vector{Int}`: normalised and scaled L₂ norm in αβγ coordinates of the degree of nonlinearity (or anti-windup) of the voltage controller [A]
- `v_err::Vector{Int}`: L₂ norm in αβγ coordinates of the voltage error signal (measured voltage subtracted from reference voltage in DQ0 coordinates) [V]
- `v_err_t::Vector{Int}`: L₂ norm in αβγ coordinates of the integrated voltage error signal [Vs]

# Note
For quantities denoted by "..._mag_..." for balanced symmetrical networks in steady state this quantity equals the root-mean-square (rms) of the signal.
"""
function RenderHookResults(; hook, states_to_plot = nothing, actions_to_plot = nothing ,
    plot_reward = false, plot_reference = false, episode = 1, vdc_to_plot = [],
    power_p_inv = [], power_q_inv = [], 
    power_p_poc = [], power_q_poc = [], 
    v_mag_inv = [], v_mag_cap = [], 
    i_mag_inv = [], i_mag_poc = [],
    freq = [], angles = [], 
    i_sat = [], i_err = [], i_err_t = [], 
    v_sat = [], v_err = [], v_err_t = [],
    v_dq = [], i_dq = [])

    #TODO complete documentation

    if isnothing(states_to_plot)
        states_to_plot = hook.collect_state_ids
    end

    if isnothing(actions_to_plot)
        actions_to_plot = hook.collect_action_ids
    end

    if isnothing(episode)
        df = hook.df
    else
        df = hook.df[hook.df.episode .== episode, :]
    end

    layout = Layout(
        plot_bgcolor="#f1f3f7",
        #title = "Results<br><sub>Run with Behavior-Actor-NNA from Episode " * string(hook.bestepisode) * "</sub>",
        xaxis_title = "Time in Seconds",
        yaxis_title = "State values",
        yaxis2 = attr(
            title="Action values",
            overlaying="y",
            side="right",
            titlefont_color="orange",
            #range=[-1, 1]
        ),
        legend = attr(
            x=1,
            y=1.02,
            yanchor="bottom",
            xanchor="right",
            orientation="h"
        ),
        width = 800,
        height = 550,
        margin=attr(l=100, r=80, b=80, t=100, pad=10)
    )
    
    
    traces = []
    
    for state_id in states_to_plot
        push!(traces, scatter(df, x = :time, y = Symbol(state_id), mode="lines", name = state_id))
    end
    
    for action_id in actions_to_plot
        push!(traces, scatter(df, x = :time, y = Symbol(action_id), mode="lines", name = action_id, yaxis = "y2"))
    end

    for vdc_idx in vdc_to_plot
        #TODO: If the index is not collected this function does print only the steps 1,2,3,4.... on y axis, WHY? How to handle that?
        push!(traces, scatter(df, x = :time, y = Symbol("source$(vdc_idx)_vdc"), mode="lines", name = "source$(vdc_idx)_vdc"))
    end

    if findfirst(x -> x == "classic", hook.policy_names) !== nothing

        for idx in hook.debug
            push!(traces, scatter(df, x = :time, y = Symbol("debug_$(idx)"), mode="lines", name = "debug_$(idx)"))
        end
        
        for idx in v_dq #hook.collect_vdq_ids #
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_v_d"), mode="lines", name = "source$(idx)_v_d"))
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_v_q"), mode="lines", name = "source$(idx)_v_q"))
        end

        for idx in i_dq #hook.collect_idq_ids #
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_i_d"), mode="lines", name = "source$(idx)_i_d"))
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_i_q"), mode="lines", name = "source$(idx)_i_q"))
        end

        for idx in power_p_inv #hook.collect_pq_ids #
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_p_inv"), mode="lines", name = "source$(idx)_p_inv"))
        end

        for idx in power_p_poc #hook.collect_pq_ids #
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_p_poc"), mode="lines", name = "source$(idx)_p_poc"))
        end

        for idx in power_q_inv #hook.collect_pq_ids #
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_q_inv"), mode="lines", name = "source$(idx)_q_inv"))
        end

        for idx in power_q_poc #hook.collect_pq_ids #
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_q_poc"), mode="lines", name = "source$(idx)_q_poc"))
        end

        for idx in v_mag_inv #hook.collect_vrms_ids #
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_v_mag_inv"), mode="lines", name = "source$(idx)_v_mag_inv"))
        end

        for idx in v_mag_cap #hook.collect_vrms_ids #
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_v_mag_cap"), mode="lines", name = "source$(idx)_v_mag_cap"))
        end

        for idx in i_mag_inv #hook.collect_irms_ids #
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_i_mag_inv"), mode="lines", name = "source$(idx)_i_mag_inv"))
        end

        for idx in i_mag_poc #hook.collect_irms_ids #
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_i_mag_poc"), mode="lines", name = "source$(idx)_i_mag_poc"))
        end

        for idx in freq
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_freq"), mode="lines", name = "source$(idx)_freq"))
        end

        for idx in angles
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_θ"), mode="lines", name = "source$(idx)_θ"))
        end

        for idx in i_sat
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_i_sat"), mode="lines", name = "source$(idx)_i_sat"))
        end

        for idx in i_err
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_i_err"), mode="lines", name = "source$(idx)_i_err"))
        end

        for idx in i_err_t
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_i_err_t"), mode="lines", name = "source$(idx)_i_err_t"))
        end

        for idx in v_sat
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_v_sat"), mode="lines", name = "source$(idx)_v_sat"))
        end

        for idx in v_err
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_v_err"), mode="lines", name = "source$(idx)_v_err"))
        end

        for idx in v_err_t
            push!(traces, scatter(df, x = :time, y = Symbol("source$(idx)_v_err_t"), mode="lines", name = "source$(idx)_v_err_t"))
        end
    end


    
    if plot_reference
        #TODO: how to check which refs to plot? 
        push!(traces, scatter(df, x = :time, y = :reference_1, mode="lines", name = "Reference"))
        push!(traces, scatter(df, x = :time, y = :reference_2, mode="lines", name = "Reference"))
        push!(traces, scatter(df, x = :time, y = :reference_3, mode="lines", name = "Reference"))
    end

    if plot_reward
        push!(traces, scatter(df, x = :time, y = :reward, yaxis = "y2", mode="lines", name = "Reward"))
    end
    
    traces = Array{GenericTrace}(traces)
    
    p = plot(traces, layout, config = PlotConfig(scrollZoom=true))
    display(p)

end

function RenderPSource(;env, hook, episode, source_ids)
    layout = Layout(
        plot_bgcolor="#f1f3f7",
        xaxis_title = "Time in Seconds",
        yaxis_title = "Power / W",
        width = 1000,
        height = 650,
        margin=attr(l=100, r=80, b=80, t=100, pad=10)
    )

    stepinteval=((episode-1)*env.maxsteps)+1:((episode)*env.maxsteps)
    time = hook.df[stepinteval, :time]
    powers=[]
    if env.nc.parameters["grid"]["phase"] === 1

        
        for id in source_ids

            if id <= env.nc.num_fltr_LCL
                push!(powers, PlotlyJS.scatter(; x = time, y = hook.df[stepinteval, Symbol("i_$id")] .* hook.df[stepinteval , Symbol("u_$id")], mode="lines", name = "Power_Source_$id"*"_Ep_$episode"))
                
            elseif id <= env.nc.num_fltr_LCL+ env.nc.num_fltr_LC
                push!(powers, PlotlyJS.scatter(; x = time, y = (hook.df[stepinteval, Symbol("i_f$id")] - hook.df[stepinteval, Symbol("op_u_f$id")]).* hook.df[stepinteval , Symbol("u_$id")], mode="lines", name = "Power_Source_$id"))
            
            elseif id <= env.nc.num_fltr_LCL+ env.nc.num_fltr_LC+ env.nc.num_fltr_L
                push!(powers, PlotlyJS.scatter(; x = time, y = hook.df[stepinteval, Symbol("i_$id")] .* hook.df[stepinteval , Symbol("u_$id")], mode="lines", name = "Power_Source_$id"))
            else
                throw("Expect sourc_ids to correspond to the amount of sources, not $id")

            end
        end

    elseif env.nc.parameters["grid"]["phase"] === 3
       
        for id in source_ids

            if id <= env.nc.num_fltr_LCL
                power = (hook.df[stepinteval , Symbol("i_$id"*"_a")] .* hook.df[stepinteval , Symbol("u_$id"*"_a")])+(hook.df[stepinteval, Symbol("i_$id"*"_b")] .* hook.df[stepinteval , Symbol("u_$id"*"_b")])+(hook.df[stepinteval , Symbol("i_$id"*"_c")] .* hook.df[stepinteval , Symbol("u_$id"*"_c")])
                push!(powers, PlotlyJS.scatter(; x = time, y = power, mode="lines", name = "Power_Source_$id"))

            elseif id <= env.nc.num_fltr_LCL+ env.nc.num_fltr_LC
                power = (hook.df[stepinteval , Symbol("i_f$id"*"_a")] - hook.df[stepinteval , Symbol("op_u_f$id"*"_a")]).* hook.df[stepinteval , Symbol("u_$id"*"_a")]+ (hook.df[stepinteval , Symbol("i_f$id"*"_b")] - hook.df[stepinteval , Symbol("op_u_f$id"*"_b")]).* hook.df[stepinteval , Symbol("u_$id"*"_b")]+(hook.df[stepinteval , Symbol("i_f$id"*"_c")] - hook.df[stepinteval , Symbol("op_u_f$id"*"_c")]).* hook.df[stepinteval , Symbol("u_$id"*"_c")]
                push!(powers, PlotlyJS.scatter(; x = time, y = power, mode="lines", name = "Power_Source_$id"))
            
            elseif id <= env.nc.num_fltr_LCL+ env.nc.num_fltr_LC+ env.nc.num_fltr_L
                power = hook.df[stepinteval , Symbol("i_$id"*"_a")] .* hook.df[stepinteval , Symbol("u_$id"*"_a")]+hook.df[stepinteval , Symbol("i_$id"*"_b")] .* hook.df[stepinteval , Symbol("u_$id"*"_b")]+hook.df[stepinteval, Symbol("i_$id"*"_c")] .* hook.df[stepinteval, Symbol("u_$id"*"_c")]
                push!(powers, PlotlyJS.scatter(; x = time, y = power, mode="lines", name = "Power_Source_$id"))
            else
                throw("Expect source_ids to correspond to the amount of sources, not $id")

            end
        end
    end

    powers = Array{GenericTrace}(powers)
    
    p = PlotlyJS.plot(powers, layout, config = PlotConfig(scrollZoom=true))
    display(p)


end

function RenderPLoad(;env, hook, episode, load_ids)
    layout = Layout(
        plot_bgcolor="#f1f3f7",
        xaxis_title = "Time in Seconds",
        yaxis_title = "Power / W",
        width = 1000,
        height = 650,
        margin=attr(l=100, r=80, b=80, t=100, pad=10)
    )

    stepinteval=((episode-1)*env.maxsteps)+1:((episode)*env.maxsteps)
    time = hook.df[stepinteval, :time]
    powers=[]
    if env.nc.parameters["grid"]["phase"] === 1

        for id in load_ids

            if id <= env.nc.num_loads_RLC
                power= hook.df[stepinteval, Symbol("u_load$id")] .* (hook.df[stepinteval, Symbol("i_load$id")] + hook.df[stepinteval, Symbol("u_load$id")] *(env.nc.parameters["load"][id]["R"])^(-1)+ (env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+2*id-1])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id")])
                push!(powers, PlotlyJS.scatter(; x = time, y = power, mode="lines", name = "Power_Load_$id"*"_Ep_$episode"))
                
            elseif id <= env.nc.num_loads_RLC+ env.nc.num_loads_LC
                power= hook.df[stepinteval, Symbol("u_load$id")] .* (hook.df[stepinteval, Symbol("i_load$id")] + (env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+2*id-1])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id")])
                push!(powers, PlotlyJS.scatter(; x = time, y = power , mode="lines", name = "Power_Load_$id"*"_Ep_$episode"))
            
            elseif id <= env.nc.num_loads_RLC+ env.nc.num_loads_LC+ env.nc.num_loads_RL
                power= hook.df[stepinteval, Symbol("u_load$id")] .* (hook.df[stepinteval, Symbol("i_load$id")] + hook.df[stepinteval, Symbol("u_load$id")] *(env.nc.parameters["load"][id]["R"])^(-1))
                push!(powers, PlotlyJS.scatter(; x = time, y = power, mode="lines", name = "Power_Load_$id"*"_Ep_$episode"))
            
            elseif id <= env.nc.num_loads_RLC+ env.nc.num_loads_LC+ env.nc.num_loads_RL+ env.nc.num_loads_L
                power= hook.df[stepinteval, Symbol("u_load$id")] .* hook.df[stepinteval, Symbol("i_load$id")]
                push!(powers, PlotlyJS.scatter(; x = time, y = power, mode="lines", name = "Power_Load_$id"*"_Ep_$episode"))
            
            elseif id <= env.nc.num_loads_RLC+ env.nc.num_loads_LC+ env.nc.num_loads_RL+ env.nc.num_loads_L+ env.nc.num_loads_RC
                power= hook.df[stepinteval, Symbol("u_load$id")] .* (hook.df[stepinteval, Symbol("u_load$id")] *(env.nc.parameters["load"][id]["R"])^(-1)+ (env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+(env.nc.num_loads_RLC + env.nc.num_loads_LC + env.nc.num_loads_RL + env.nc.num_loads_L)+id])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id")])
            
            elseif id <= env.nc.num_loads_RLC+ env.nc.num_loads_LC+ env.nc.num_loads_RL+ env.nc.num_loads_L+ env.nc.num_loads_RC+ env.nc.num_loads_C
                power= hook.df[stepinteval, Symbol("u_load$id")] .* ((env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+(env.nc.num_loads_RLC + env.nc.num_loads_LC + env.nc.num_loads_RL + env.nc.num_loads_L)+id])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id")])
            
            elseif id <= env.nc.num_loads_RLC+ env.nc.num_loads_LC+ env.nc.num_loads_RL+ env.nc.num_loads_L+ env.nc.num_loads_RC+ env.nc.num_loads_C+ env.nc.num_loads_R
                power= hook.df[stepinteval, Symbol("u_load$id")] .* (hook.df[stepinteval, Symbol("u_load$id")] *(env.nc.parameters["load"][id]["R"])^(-1))

            else
                throw("Expect laod_ids to correspond to the amount of loads, not $id")

            end
        end
        

    elseif env.nc.parameters["grid"]["phase"] === 3
        
        for id in load_ids

            if id <= env.nc.num_loads_RLC
                power_a= hook.df[stepinteval, Symbol("u_load$id"*"_a")] .* (hook.df[stepinteval, Symbol("i_load$id"*"_a")] + hook.df[stepinteval, Symbol("u_load$id"*"_a")] *(env.nc.parameters["load"][id]["R"])^(-1)+ (env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+2*id-1])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id"*"_a")])
                power_b= hook.df[stepinteval, Symbol("u_load$id"*"_b")] .* (hook.df[stepinteval, Symbol("i_load$id"*"_b")] + hook.df[stepinteval, Symbol("u_load$id"*"_b")] *(env.nc.parameters["load"][id]["R"])^(-1)+ (env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+2*id-1])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id"*"_b")])
                power_c= hook.df[stepinteval, Symbol("u_load$id"*"_c")] .* (hook.df[stepinteval, Symbol("i_load$id"*"_c")] + hook.df[stepinteval, Symbol("u_load$id"*"_c")] *(env.nc.parameters["load"][id]["R"])^(-1)+ (env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+2*id-1])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id"*"_c")])
                push!(powers, PlotlyJS.scatter(; x = time, y = power_a + power_b + power_c, mode="lines", name = "Power_Load_$id"*"_Ep_$episode"))
                
            elseif id <= env.nc.num_loads_RLC+ env.nc.num_loads_LC
                power_a= hook.df[stepinteval, Symbol("u_load$id"*"_a")] .* (hook.df[stepinteval, Symbol("i_load$id"*"_a")] + (env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+2*id-1])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id"*"_a")])
                power_b= hook.df[stepinteval, Symbol("u_load$id"*"_b")] .* (hook.df[stepinteval, Symbol("i_load$id"*"_b")] + (env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+2*id-1])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id"*"_b")])
                power_c= hook.df[stepinteval, Symbol("u_load$id"*"_c")] .* (hook.df[stepinteval, Symbol("i_load$id"*"_c")] + (env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+2*id-1])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id"*"_c")])
                push!(powers, PlotlyJS.scatter(; x = time, y = power_a + power_b + power_c , mode="lines", name = "Power_Load_$id"*"_Ep_$episode"))
            
            elseif id <= env.nc.num_loads_RLC+ env.nc.num_loads_LC+ env.nc.num_loads_RL
                power_a= hook.df[stepinteval, Symbol("u_load$id"*"_a")] .* (hook.df[stepinteval, Symbol("i_load$id"*"_a")] + hook.df[stepinteval, Symbol("u_load$id"*"_a")] *(env.nc.parameters["load"][id]["R"])^(-1))
                power_b= hook.df[stepinteval, Symbol("u_load$id"*"_b")] .* (hook.df[stepinteval, Symbol("i_load$id"*"_b")] + hook.df[stepinteval, Symbol("u_load$id"*"_b")] *(env.nc.parameters["load"][id]["R"])^(-1))
                power_c= hook.df[stepinteval, Symbol("u_load$id"*"_c")] .* (hook.df[stepinteval, Symbol("i_load$id"*"_c")] + hook.df[stepinteval, Symbol("u_load$id"*"_c")] *(env.nc.parameters["load"][id]["R"])^(-1))
                push!(powers, PlotlyJS.scatter(; x = time, y = power_a + power_b + power_c , mode="lines", name = "Power_Load_$id"*"_Ep_$episode"))
            
            elseif id <= env.nc.num_loads_RLC+ env.nc.num_loads_LC+ env.nc.num_loads_RL+ env.nc.num_loads_L
                power_a= hook.df[stepinteval, Symbol("u_load$id"*"_a")] .* hook.df[stepinteval, Symbol("i_load$id"*"_a")]
                power_b= hook.df[stepinteval, Symbol("u_load$id"*"_b")] .* hook.df[stepinteval, Symbol("i_load$id"*"_b")]
                power_c= hook.df[stepinteval, Symbol("u_load$id"*"_c")] .* hook.df[stepinteval, Symbol("i_load$id"*"_c")]
                push!(powers, PlotlyJS.scatter(; x = time, y = power_a + power_b + power_c, mode="lines", name = "Power_Load_$id"*"_Ep_$episode"))
            
            elseif id <= env.nc.num_loads_RLC+ env.nc.num_loads_LC+ env.nc.num_loads_RL+ env.nc.num_loads_L+ env.nc.num_loads_RC
                power_a= hook.df[stepinteval, Symbol("u_load$id"*"_a")] .* (hook.df[stepinteval, Symbol("u_load$id"*"_a")] *(env.nc.parameters["load"][id]["R"])^(-1)+ (env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+(env.nc.num_loads_RLC + env.nc.num_loads_LC + env.nc.num_loads_RL + env.nc.num_loads_L)+id])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id"*"_a")])
                power_b= hook.df[stepinteval, Symbol("u_load$id"*"_b")] .* (hook.df[stepinteval, Symbol("u_load$id"*"_b")] *(env.nc.parameters["load"][id]["R"])^(-1)+ (env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+(env.nc.num_loads_RLC + env.nc.num_loads_LC + env.nc.num_loads_RL + env.nc.num_loads_L)+id])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id"*"_b")])
                power_c= hook.df[stepinteval, Symbol("u_load$id"*"_c")] .* (hook.df[stepinteval, Symbol("u_load$id"*"_c")] *(env.nc.parameters["load"][id]["R"])^(-1)+ (env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+(env.nc.num_loads_RLC + env.nc.num_loads_LC + env.nc.num_loads_RL + env.nc.num_loads_L)+id])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id"*"_c")])
                push!(powers, PlotlyJS.scatter(; x = time, y = power_a + power_b + power_c, mode="lines", name = "Power_Load_$id"*"_Ep_$episode"))

            elseif id <= env.nc.num_loads_RLC+ env.nc.num_loads_LC+ env.nc.num_loads_RL+ env.nc.num_loads_L+ env.nc.num_loads_RC+ env.nc.num_loads_C
                power_a= hook.df[stepinteval, Symbol("u_load$id"*"_a")] .* ((env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+(env.nc.num_loads_RLC + env.nc.num_loads_LC + env.nc.num_loads_RL + env.nc.num_loads_L)+id])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id"*"_a")])
                power_b= hook.df[stepinteval, Symbol("u_load$id"*"_b")] .* ((env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+(env.nc.num_loads_RLC + env.nc.num_loads_LC + env.nc.num_loads_RL + env.nc.num_loads_L)+id])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id"*"_b")])
                power_c= hook.df[stepinteval, Symbol("u_load$id"*"_c")] .* ((env.nc.parameters["load"][id]["C"])*(hook.collect_state_paras[env.nc.num_fltr+env.nc.num_connections+(env.nc.num_loads_RLC + env.nc.num_loads_LC + env.nc.num_loads_RL + env.nc.num_loads_L)+id])^(-1) * hook.df[stepinteval, Symbol("op_u_load$id"*"_c")])
                push!(powers, PlotlyJS.scatter(; x = time, y = power_a + power_b + power_c, mode="lines", name = "Power_Load_$id"*"_Ep_$episode"))

            elseif id <= env.nc.num_loads_RLC+ env.nc.num_loads_LC+ env.nc.num_loads_RL+ env.nc.num_loads_L+ env.nc.num_loads_RC+ env.nc.num_loads_C+ env.nc.num_loads_R
                power_a= hook.df[stepinteval, Symbol("u_load$id"*"_a")] .* (hook.df[stepinteval, Symbol("u_load$id"*"_a")] *(env.nc.parameters["load"][id]["R"])^(-1))
                power_b= hook.df[stepinteval, Symbol("u_load$id"*"_b")] .* (hook.df[stepinteval, Symbol("u_load$id"*"_b")] *(env.nc.parameters["load"][id]["R"])^(-1))
                power_c= hook.df[stepinteval, Symbol("u_load$id"*"_c")] .* (hook.df[stepinteval, Symbol("u_load$id"*"_c")] *(env.nc.parameters["load"][id]["R"])^(-1))
                push!(powers, PlotlyJS.scatter(; x = time, y = power_a + power_b + power_c, mode="lines", name = "Power_Load_$id"*"_Ep_$episode"))
                
            else
                throw("Expect laod_ids to correspond to the amount of loads, not $id")

            end
        end

        
    end

    powers = Array{GenericTrace}(powers)
    
    p = PlotlyJS.plot(powers, layout, config = PlotConfig(scrollZoom=true))
    display(p)

end