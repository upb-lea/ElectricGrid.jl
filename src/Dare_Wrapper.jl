using ReinforcementLearning

mutable struct dare_setup
    env
    control_agent
    hook
end

function create_setup(;num_sources, num_loads, CM = nothing, parameters = nothing, ts = 100e-6, t_end = 0.5, env_cuda = false, agent_cuda = false)

    #_______________________________________________________________________________
    # State space representation

    #-------------------------------------------------------------------------------
    # Connectivity Matrix

    #= CM = [ 0. 0. 0. 1.
            0. 0. 0. 2.
            0. 0. 0. 3.
            -1. -2. -3. 0.] =#

    CM = [ 0. 0. 1.
    0. 0. 2
    -1. -2. 0.]

    #= CM = [0. 1.
    -1. 0.] =#

    #-------------------------------------------------------------------------------
    # Cables

    cable_list = []

    # Network Cable Impedances
    l = 2.5 # length in km
    cable = Dict()
    cable["R"] = 0.208*l # Ω, line resistance 0.722#
    cable["L"] = 0.00025*l # H, line inductance 0.264e-3#
    cable["C"] = 0.4e-6*l # 0.4e-6#

    #= push!(cable_list, cable, cable, cable) =#

    push!(cable_list, cable, cable)

    #-------------------------------------------------------------------------------
    # Sources
    source = Dict()

    source_list = []

    source["pwr"] = 200e3
    source["vdc"] = 800
    source["fltr"] = "LC"
    source["p_set"] = 50e3
    source["q_set"] = 10e3
    source["v_pu_set"] = 1.05

    push!(source_list, source)

    source["pwr"] = 100e3
    source["vdc"] = 800
    source["fltr"] = "LC"
    source["p_set"] = 50e3
    source["q_set"] = 10e3
    source["v_pu_set"] = 1.0

    push!(source_list, source)

    #= source = Dict()

    source["pwr"] = 200e3
    source["vdc"] = 800
    source["fltr"] = "LC"
    source["i_rip"] = 0.15
    source["v_rip"] = 0.01537

    push!(source_list, source) =#

    #-------------------------------------------------------------------------------
    # Loads

    load_list = []
    load = Dict()

    R1_load, L_load, _, _ = Load_Impedance_2(50e3, 0.6, 230)
    #R2_load, C_load, _, _ = Load_Impedance_2(150e3, -0.8, 230)

    load["impedance"] = "RL"
    load["R"] = R1_load# + R2_load # 
    load["L"] = L_load
    #load["C"] = C_load

    push!(load_list, load)

    #-------------------------------------------------------------------------------
    # Amalgamation

    parameters = Dict()

    parameters["source"] = source_list
    parameters["cable"] = cable_list
    parameters["load"] = load_list
    parameters["grid"] = Dict("fs" => fs, "phase" => 3, "v_rms" => 230)

    # Define the environment

    num_sources = length(source_list)
    num_loads = length(load_list)
        
    maxsteps = Int(t_end / ts)

    function reward(env, name = nothing)
        u_l1_index = findfirst(x -> x == "source1_v_C", env.state_ids)

        u_l1 = env.state[u_l1_index]

        ref = sqrt(2)*230 * cos.(2*pi*50*env.t)

        r = -abs(ref/600 - u_l1)/3 * (1 - 0.99)
    
        return r
    end

    env = SimEnv(reward_function = reward, num_sources = num_sources, num_loads = num_loads, CM = CM, parameters = parameters, ts = ts, maxsteps = maxsteps, use_gpu = env_cuda)

    ns = length(env.sys_d.A[1,:])
    na = length(env.sys_d.B[1,:])

    agent = create_agent_ddpg(na = na, ns = ns, use_gpu = env_cuda)

    hook = DataHook(save_best_NNA = true, plot_rewards = true, collect_state_ids = env.state_ids, collect_action_ids = [])

    dare_setup(env, agent, hook)
end

function run(dare_setup, no_episodes)
    RLBase.run(dare_setup.control_agent, dare_setup.env, StopAfterEpisode(no_episodes), dare_setup.hook)
end