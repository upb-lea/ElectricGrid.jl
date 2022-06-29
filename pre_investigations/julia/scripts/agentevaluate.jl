rewards = []
actions = []

function execute_env(env::SimEnv, agent::Agent, t_len::Int, debug::Bool)
    if debug
        output = zeros(length(env.Ad[1,:]), t_len+1)
    else
        output = 0.0
    end

    RLBase.reset!(env)

    for i = 1:t_len
        action = agent(env)
        env(action) 
        push!(rewards, reward(env))
        # normalised action - > denormalize
        push!(actions, action * 300)
        if debug 
            output[:,i+1] = env.state.*env.norm_array 
        end
    end

    return output
end

# params = No_Episodes, 

agent.policy.act_noise = 0.0
result = execute_env(env, agent, 300, true)

p =plot(actions, 
    title = "Actions in actor mode",
    ylabel = "Voltage (V)",
    xlabel = "Time steps",
    legend = false)
plot!(result[end,:])
plot!(result[2,:], label = ["v_source" "v_load" "v_capacitor"])
display(p)

plot(rewards, 
    title = "Rewards over each time step [actor mode]",
    ylabel = "normalised rewards",
    xlabel = "Time steps",
    legend = false)


p = plot(actions .* result[1,2:end], 
    title = "source power",
    xlabel = "Time steps",
    ylabel = "power in Watts")

Actor_pload = result[end, :] .^2 / 14
plot!(Actor_pload) #, title = "actual load")
display(p)

p =plot(result[1, :], 
    title = "Currents in actor mode",
    ylabel = "Current (A)",
    xlabel = "Time steps",
    legend = false)

plot!(result[3,:])
