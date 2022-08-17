using DataFrames
using PlotlyJS

include(srcdir("data_hook.jl"))

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

function plot_best_results(;agent, env, hook, state_ids_to_plot)
    reset!(env)

    act_noise_old = agent.policy.act_noise
    agent.policy.act_noise = 0.0
    copyto!(agent.policy.behavior_actor, hook.bestNNA)
    
    temphook = DataHook(collect_state_ids = state_ids_to_plot)
    
    run(agent.policy, env, StopAfterEpisode(1), temphook)
    
    layout = Layout(
            plot_bgcolor="#f1f3f7",
            title = "Results<br><sub>Run with Behavior-Actor-NNA from Episode " * string(hook.bestepisode) * "</sub>",
            xaxis_title = "Time in Seconds",
            yaxis_title = "State values",
            yaxis2 = attr(
                title="Reward",
                overlaying="y",
                side="right",
                titlefont_color="orange",
                range=[-1, 1]
            ),
            legend = attr(
                x=1,
                y=1.02,
                yanchor="bottom",
                xanchor="right",
                orientation="h"
            ),
            width = 1000,
            height = 650,
            margin=attr(l=100, r=80, b=80, t=100, pad=10)
        )
    
    trace1 = scatter(temphook.df, x = :time, y = :u_f1, mode="lines", name = "U f1")
    trace2 = scatter(temphook.df, x = :time, y = :u_1, mode="lines", name = "U 1")
    trace3 = scatter(temphook.df, x = :time, y = :u_2, mode="lines", name = "U 2")
    trace4 = scatter(temphook.df, x = :time, y = :u_l1, mode="lines", name = "U l1")
    trace5 = scatter(temphook.df, x = :time, y = :reward, yaxis = "y2", mode="lines", name = "Reward")
    
    p = plot([trace1, trace2, trace3, trace4, trace5], layout, config = PlotConfig(scrollZoom=true))
    display(p)
    
    copyto!(agent.policy.behavior_actor, hook.currentNNA)
    agent.policy.act_noise = act_noise_old
    
    reset!(env)
end

function get_data(episode)

    data = CSV.read("episode_data/$(episode).csv", DataFrame)
    data.state = eval.(Meta.parse.(data.state))
    data.next_state = eval.(Meta.parse.(data.next_state))
    data.action = eval.(Meta.parse.(data.action))
    
    data

end

function plotting_state(episodes, states, NodeConstructor) # also path to data
    
    plots = Array{Plots.Plot}(undef, length(states), length(episodes))

    phase = NodeConstructor.parameters["grid"]["phase"]
    fs = NodeConstructor.parameters["grid"]["fs"]

    state_list = get_states(NodeConstructor)

    for (idx_ep, episode) in enumerate(episodes)

        data = get_data(episode)
        t = range(0, length(data.state)-1)/fs |> collect

        for (idx_state, state) in enumerate(states)
                
            if typeof(state) === String
                state = findall(x->x==state, state_list)
                state = state[1]
            end
            
            data_state = mapreduce(permutedims, vcat, data.state)
            
            ns = NodeConstructor.num_spp
            
            plots[idx_state,idx_ep] = plot(xlabel="time", ylabel="$(state_list[state])")
            
            if phase == 1
                plot!(plots[idx_state,idx_ep], t, data_state[:, state+ns*0], legend=false)
            
            elseif phase == 3
                plot!(plots[idx_state,idx_ep], t, data_state[:, state+ns*0], label="Phase 1")
                plot!(plots[idx_state,idx_ep], t, data_state[:, state+ns*1], label="Phase 2")
                plot!(plots[idx_state,idx_ep], t, data_state[:, state+ns*2], label="Phase 3")
                    
            end
        end
    end

    plot(plots..., layout=(length(episodes), length(states)))

end

function plot_p_source(episode, source_id, NodeConstructor)

    data = get_data(episode)
    phase = NodeConstructor.parameters["grid"]["phase"]
    fs = NodeConstructor.parameters["grid"]["fs"]

    t = range(0, length(data.state)-1)/fs |> collect

    data_state = mapreduce(permutedims, vcat, data.state)
    ns = NodeConstructor.num_spp

    if source_id <= NodeConstructor.num_fltr_LCL
        LCL_source_id = source_id
        i_state = (LCL_source_id-1)*4 + 3
        u_state = (LCL_source_id-1)*4 + 4

    elseif source_id <= NodeConstructor.num_fltr_LCL+ NodeConstructor.num_fltr_LC
        LC_source_id = source_id - NodeConstructor.num_fltr_LCL
        i_state = NodeConstructor.num_fltr_LCL*4+(LC_source_id-1)*4 + 1
        u_state = NodeConstructor.num_fltr_LCL*4+(LC_source_id-1)*4 + 3
    
    elseif source_id <= NodeConstructor.num_fltr_LCL+ NodeConstructor.num_fltr_LC+ NodeConstructor.num_fltr_L
        L_source_id = source_id - NodeConstructor.num_fltr_LCL - NodeConstructor.num_fltr_LC
        i_state = NodeConstructor.num_fltr_LCL*4+NodeConstructor.num_fltr_LC*3+(L_source_id-1)*4 + 1
        u_state = NodeConstructor.num_fltr_LCL*4+NodeConstructor.num_fltr_LC*3+(L_source_id-1)*4 + 2

    else
        throw("Expect id to be suiting to the amount of sources , not $source_id")

    end

    
    state_list = get_states(NodeConstructor)
    println(state_list[i_state])
    println(state_list[u_state])
    
    # isum_computation

    CM_row = NodeConstructor.CM[source_id,:]
    indizes = CM_row[CM_row .!= 0]

    i_source_list= []
    i_source= zeros(length(data_state[:,1]))
    for idx in indizes
        idx = Int(idx)
        i_source += sign(idx)*data_state[:, NodeConstructor.num_fltr + abs(idx)]
        append!(i_source_list, (-1)*idx) 
    end

    println(i_source_list)
    
    power_i_sum = data_state[:, u_state] .* i_source

    
    p = plot(xlabel="time", ylabel="power_source_$source_id")


    if phase == 1
        
        power= data_state[:, i_state] .* data_state[:, u_state]
        plot!(p, t, power, label="u1*i1_computation")
        plot!(p, t, power_i_sum, label="u1*isum_computation")

    elseif phase == 3
            
        power1= data_state[:, i_state+ns*0] .* data_state[:, u_state+ns*0]
        power2= data_state[:, i_state+ns*1] .* data_state[:, u_state+ns*1]
        power3= data_state[:, i_state+ns*2] .* data_state[:, u_state+ns*2]

        plot!(p, t, power1, label="Phase 1")
        plot!(p, t, power2, label="Phase 2")
        plot!(p, t, power3, label="Phase 3")
            
    end
end

function plot_p_load(episode, load_id, NodeConstructor)
    
    data = get_data(episode)
    phase = NodeConstructor.parameters["grid"]["phase"]
    fs = NodeConstructor.parameters["grid"]["fs"]

    t = range(0, length(data.state)-1)/fs |> collect

    data_state = mapreduce(permutedims, vcat, data.state)
    ns = NodeConstructor.num_spp
    state_list = get_states(NodeConstructor)

    load_startpos =  NodeConstructor.num_fltr + NodeConstructor.num_connections
    
    num_twostateloads = NodeConstructor.num_loads_RLC + NodeConstructor.num_loads_LC + NodeConstructor.num_loads_RL + NodeConstructor.num_loads_L
    num_onestateloads = NodeConstructor.num_loads_RC + NodeConstructor.num_loads_C + NodeConstructor.num_loads_R

    if load_id <=  num_twostateloads
        u_state = load_startpos + (load_id-1)*2 + 1

    elseif load_id <=  num_twostateloads+num_onestateloads
        u_load_id = load_id - num_twostateloads
        u_state = num_twostateloads*2+u_load_id
    
    else

        throw("Expect id to be suiting to the amount of loads , not $load_id")
    end
    #println(state_list)
    println(state_list[u_state])
   
    CM_row = NodeConstructor.CM[NodeConstructor.num_sources+load_id,:]
    indizes = CM_row[CM_row .!= 0]
    # signs = [sign(x) for x in indizes] # get signs
    # indizes_ = indizes .* signs # delet signs from indices

    p = plot(xlabel="time", ylabel="power_load_$load_id")
    #println(indizes_)
    

    if phase == 1
        
        i_load_list= []
        i_load= zeros(length(data_state[:,1]))
        for idx in indizes
            idx = Int(idx)
            i_load += (-1)*sign(idx)*data_state[:, NodeConstructor.num_fltr + abs(idx)]
            append!(i_load_list, (-1)*idx) 
        end

        println(i_load_list)
        power= data_state[:, u_state] .* i_load
        # plot(t,i_load, label="i_load")
        # plot(t,data_state[:, u_state], label="$(state_list[u_state])")
        plot!(p, t, power)
    
    elseif phase == 3
            
        for idx in indizes_
            idx = Int(idx)
            i_load1 += data_state[:, 0*ns+ NodeConstructor.num_fltr + idx]
            i_load2 += data_state[:, 1*ns+NodeConstructor.num_fltr + idx]
            i_load3 += data_state[:, 2*ns+NodeConstructor.num_fltr + idx]
        end
        
        power1= i_load1 .* data_state[:, u_state+ns*0]
        power2= i_load2 .* data_state[:, u_state+ns*1]
        power3= i_load3 .* data_state[:, u_state+ns*2]

        plot!(p, t, power1, label="Phase 1")
        plot!(p, t, power2, label="Phase 2")
        plot!(p, t, power3, label="Phase 3")
            
    end
end