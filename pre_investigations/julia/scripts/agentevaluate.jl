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
        push!(actions, action)
        if debug 
            output[:,i+1] = env.state.*env.norm_array 
        end
    end

    return output
end


# params = No_Episodes, 

result = execute_env(env, agent, 300, true)
display(plot(actions))
display(plot(rewards))