function reward(env)
    u_l1_index = findfirst(x -> x == "u_load1", env.state_ids)
    u_l1 = env.state[u_l1_index]

    return -(abs(230 - (env.norm_array[u_l1_index] * u_l1))/300)
end

function prepare_env()
    global env = ElectricGridEnv(CM=CM, num_sources=length(parameters["source"]), num_loads=length(parameters["load"]), parameters=parameters,
                        reward_function = reward, maxsteps=600)

    ns = length(env.sys_d.A[1,:])
    na = length(env.sys_d.B[1,:])
    
    global agent = create_agent_ddpg(na = na, ns = ns)

    global hook = data_hook(save_best_NNA = true, plot_rewards = false)

    global model.plot_agent[] = PlotlyBase.Plot()
end

function run_learning()
    layout = PlotlyBase.Layout(
            plot_bgcolor="#f1f3f7",
            title = "Total Reward per Episode",
            xaxis_title = "Episode",
            yaxis_title = "Total Reward",
            showlegend=false
        )

    for i in 1:20
        run(agent, env, StopAfterEpisode(1), hook)
        pl = PlotlyBase.Plot(scatter(x = collect(1:length(hook.rewards)), y = hook.rewards, mode="lines", name = "rewards"), layout)
        global model.plot_agent[] = pl
    end
end