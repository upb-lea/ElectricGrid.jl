using ReinforcementLearning

mutable struct dare_setup
    env
    control_agent
    hook
end

function create_setup(;num_sources, num_loads, CM=nothing, parameters=nothing, ts=0.0001, t_end=0.05, env_cuda=false, agent_cuda=false)

    CM = [ 0. 0. 1.
     0. 0. 2
     -1. -2. 0.]

    #-------------------------------------------------------------------------------
    # Cables

    cable_list = []

    # Network Cable Impedances
    l = 1.0 # length in km
    cable = Dict()
    cable["R"] = 0.208*l # Ω, line resistance 0.722#
    cable["L"] = 0.00025*l # H, line inductance 0.264e-3#
    cable["C"] = 0.4e-6*l # 0.4e-6#

    push!(cable_list, cable, cable)

    # Sources

    source_list = []
    source = Dict()
    # time step
    ts = 0.0001

    fs = 1/ts

    source["pwr"] = 200e3
    #source["vdc"] = 800
    source["fltr"] = "LC"
    Lf, Cf, _ = Filter_Design(source["pwr"], fs)
    source["R1"] = 0.4
    source["R_C"] = 0.0006
    source["L1"] = Lf
    source["C"] = Cf
    source["v_limit"]= 1500
    source["i_limit"]= 1000

    push!(source_list, source)

    source = Dict()

    source["pwr"] = 100e3
    source["vdc"] = 600
    source["fltr"] = "LC"
    Lf, Cf, _ = Filter_Design(source["pwr"], fs)
    source["R1"] = 0.4
    source["R_C"] = 0.0006
    source["L1"] = Lf
    source["C"] = Cf
    source["v_limit"]= 1500 
    source["i_limit"]= 1000

    push!(source_list, source)

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
    parameters["grid"] = Dict("fs" => fs, "phase" => 1, "v_rms" => 230)

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

function create_setup_2(;num_sources = 0, num_loads = 0, CM = nothing, parameters = nothing, 
                        t_end = 0.5, env_cuda = false, agent_cuda = false, hook = nothing)

    #-------------------------------------------------------------------------------
    # Auxiliary Functions

    function reward(env, name = nothing)
        u_l1_index = findfirst(x -> x == "source1_v_C", env.state_ids)

        u_l1 = env.state[u_l1_index]

        ref = sqrt(2)*230 * cos.(2*pi*50*env.t)

        r = -abs(ref/600 - u_l1)/3 * (1 - 0.99)
    
        return r
    end

    function featurize(x0 = nothing, t0 = nothing; env = nothing, name = nothing)

        if !isnothing(name)
            state = env.state
            if name == agentname
                global state_ids_agent
                global state_ids
                state = state[findall(x -> x in state_ids_agent, state_ids)]
                state = vcat(state, reference(env.t)/600)
            else
                global state_ids_classic
                global state_ids
                state = env.x[findall(x -> x in state_ids_classic, state_ids)]
            end
        elseif isnothing(env)
            return x0
        else
            return env.state
        end
        return state
    end

    function RLBase.action_space(env::SimEnv, name::String)
        if name == "agent"
            return Space(fill(-1.0..1.0, size(action_ids_agent)))
        else
            return Space(fill(-1.0..1.0, size(action_ids_classic)))
        end
    end

    #-------------------------------------------------------------------------------
    # Finalising the Node Constructor

    if !isempty(parameters)
        num_sources = length(source_list)
        num_loads = length(load_list)
    end

    if parameters["grid"]["fs"] === nothing
        ts = 100e-6
    else
        ts = 1/(parameters["grid"]["fs"])
    end
        
    maxsteps = Int(t_end / ts) + 1

    #-------------------------------------------------------------------------------
    # Defining the Environment

    env = SimEnv(reward_function = reward, num_sources = num_sources, num_loads = num_loads, 
                CM = CM, parameters = parameters, ts = ts, 
                maxsteps = maxsteps, use_gpu = env_cuda, action_delay = 1)
    
    state_ids = get_state_ids(env.nc)
    action_ids = get_action_ids(env.nc)

    #-------------------------------------------------------------------------------
    # Outputs
    
    if hook === nothing
        hook = DataHook(save_best_NNA = true, plot_rewards = true, 
        collect_state_ids = env.state_ids, collect_action_ids = [])
    end

    #-------------------------------------------------------------------------------
    # Setting up the controls for the Classical Sources

    Animo = NamedPolicy("classic", Classical_Policy(env))

    num_clas_sources = length(Animo.policy.Source_Indices) # number of classically controlled sources
    state_ids_classic = Animo.policy.state_ids
    action_ids_classic = Animo.policy.action_ids

    #-------------------------------------------------------------------------------
    # Setting up the controls for the Reinforcement Learning Sources

    num_RL_sources = env.nc.num_sources - num_clas_sources # number of reinforcement learning controlled sources

    if num_RL_sources == 1

        ns = length(env.sys_d.A[1,:])
        na = length(env.sys_d.B[1,:])

        agent = create_agent_ddpg(na = na, ns = ns, use_gpu = env_cuda)
        agent = Agent(policy = NamedPolicy("agent", agent.policy), trajectory = agent.trajectory)
    end

    #-------------------------------------------------------------------------------
    # Finalising the Control

    num_policies = num_RL_sources + convert(Int64, num_clas_sources > 0)

    Multi_Agents = Dict()

    if num_clas_sources > 0

        polc = Dict()

        polc["policy"] = Animo
        polc["state_ids"] = state_ids_classic
        polc["action_ids"] = action_ids_classic

        Multi_Agents[nameof(Animo)] = polc
    end
    
    Multi_Agent = MultiAgentGridController(Multi_Agents, action_ids)

    return dare_setup(env, Multi_Agent, hook)
end

function run(dare_setup, no_episodes)

    RLBase.run(dare_setup.control_agent, dare_setup.env, StopAfterEpisode(no_episodes), dare_setup.hook)

    return nothing
end