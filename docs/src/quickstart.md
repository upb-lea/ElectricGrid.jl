Consider the following code snippet


    using JEG
    
    env =  ElectricGridEnv(num_sources = 1, num_loads = 1)
    Multi_Agent =  SetupAgents(env)
    hook =  Simulate(Multi_Agent, env)
    RenderHookResults(hook = hook)

This is a minimal example of a full JuliaElectricGrid setup. After running you should see some output in your terminal, created by Ipopt, an EPL licensed library JEG is using. There should also appear a plot that looks like this:
![output of the minimal example](./assets/output1.png)