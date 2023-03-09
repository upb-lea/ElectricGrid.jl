


```@docs
ElectricGridEnv(;
        maxsteps = 500,
        ts = 1/10_000,
        action_space = nothing,
        state_space = nothing,
        prepare_action = nothing,
        featurize = nothing,
        reward_function = nothing,
        CM = nothing,
        num_sources = nothing,
        num_loads = nothing,
        parameters = nothing,
        x0 = nothing,
        t0 = 0.0,
        state_ids = nothing,
        convert_state_to_cpu = true,
        use_gpu = false,
        reward = nothing,
        action = nothing,
        action_ids = nothing,
        action_delay = 1,
        t_end = nothing,
        verbosity = 0,
        agent_dict = nothing
    )
```