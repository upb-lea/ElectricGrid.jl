using DataFrames
using PlotlyJS

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

function plot_best_results(;agent, env, hook, state_ids_to_plot, plot_reward = true, plot_reference = false)
    reset!(env)

    act_noise_old = agent.policy.act_noise
    agent.policy.act_noise = 0.0
    copyto!(agent.policy.behavior_actor, hook.bestNNA)
    
    temphook = DataHook(collect_state_ids = state_ids_to_plot, collect_reference = true)
    
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


    traces = []

    for state_id in state_ids_to_plot
        push!(traces, scatter(temphook.df, x = :time, y = Symbol(state_id), mode="lines", name = state_id))
    end

    if plot_reward
        push!(traces, scatter(temphook.df, x = :time, y = :reward, yaxis = "y2", mode="lines", name = "Reward"))
    end

    if plot_reference
        push!(traces, scatter(temphook.df, x = :time, y = :reference, mode="lines", name = "Reference"))
    end

    traces = Array{GenericTrace}(traces)
    
    p = plot(traces, layout, config = PlotConfig(scrollZoom=true))
    display(p)
    
    copyto!(agent.policy.behavior_actor, hook.currentNNA)
    agent.policy.act_noise = act_noise_old
    
    reset!(env)

    return nothing
end