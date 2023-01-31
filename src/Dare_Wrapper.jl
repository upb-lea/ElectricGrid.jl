
"""
    Power_System_Dynamics(env, hook; num_episodes = 1, return_Agents = false)

# Description
Runs the simulation with the specified parameters.

# Arguments
- `env::SimEnv`: mutable struct containing the environment.
- `hook::DataHook`: mutable struct definining the signals to be saved.

# Keyword Arguments
- `num_episodes::Int`: the number of times that the simulation is run.
- `return_Agents::Float`: returns the mutable struct containing the agents.

# Return Values
- `Multi_Agent::MultiAgentGridController`: (optional)

"""


function Power_System_Dynamics(env, hook; num_episodes = 1, return_Agents = false)


    #_______________________________________________________________________________
    # Setting up the Agents

    Agents = Dict()

    #-------------------------------------------------------------------------------
    # Setting up the controls for the Classical Sources

    Animo = NamedPolicy("classic", Classical_Policy(env))

    num_clas_sources = Animo.policy.Source.num_sources # number of classically controlled sources

    #-------------------------------------------------------------------------------
    # Setting up the controls for the Reinforcement Learning Sources

    num_RL_sources = env.nc.num_sources - num_clas_sources # number of reinforcement learning controlled sources

    

    #= for s in axes(Source_Indices, 1)

        s_idx = string(Source_Indices[s])

        Source.V_cable_loc[:, s]  = findall(contains(s_idx*"_v_C_cable"), state_ids_classic)
        Source.I_inv_loc[:, s] = findall(contains(s_idx*"_i_L1"), state_ids_classic)

        if Source.filter_type[s] == "LC"

            Source.V_cap_loc[:, s]  = findall(contains(s_idx*"_v_C_filt"), state_ids_classic)
            Source.I_poc_loc[:, s] = findall(contains(s_idx*"_v_C_filt"), state_ids)

        elseif Source.filter_type[s] == "LCL"

            Source.V_cap_loc[:, s]  = findall(contains(s_idx*"_v_C_filt"), state_ids_classic)
            Source.I_poc_loc[:, s] = findall(contains(s_idx*"_i_L2"), state_ids_classic)
        
        elseif Source.filter_type[s] == "L"

            Source.I_poc_loc[:, s] = findall(contains(s_idx*"_v_C_cable"), state_ids)   # TODO check _v_C_??
        end
    end

    letterdict = Dict("a" => 1, "b" => 2, "c" => 3)

    Source.Action_loc = [[findfirst(y -> y == parse(Int64, SubString(split(x, "_")[1], 7)), 
    Source_Indices), letterdict[split(x, "_")[3]]] for x in action_ids_classic] =#

    RL_Source_Indices = Array{Int64, 1}(undef, 0)
    for ns in 1:env.nc.num_sources

        if env.nc.parameters["source"][ns]["control_type"] == "RL"
            RL_Source_Indices = [RL_Source_Indices; Int(ns)]
        end
    end

    ssa = "source".*string.(RL_Source_Indices)
    env.state_ids_RL = filter(x -> !isempty(findall(y -> y == split(x, "_")[1], ssa)), env.state_ids)     
    env.action_ids_RL = filter(x -> !isempty(findall(y -> y == split(x, "_")[1], ssa)), env.action_ids)

    if num_RL_sources > 0

        ns = length(env.state_ids) - length(Animo.policy.state_ids)  # all RL_type sources are controlled by 1 RL_agent (#TODO: spilt up to more then 1 agent?)
        na = length(env.action_ids) - length(Animo.policy.action_ids)

        agent = create_agent_ddpg(na = na, ns = length(env.state_ids_RL), use_gpu = false)
        agent = Agent(policy = NamedPolicy("agent", agent.policy), trajectory = agent.trajectory)

        RL_policy = Dict()
        RL_policy["policy"] = agent
        RL_policy["state_ids"] = env.state_ids_RL
        RL_policy["action_ids"] = env.action_ids_RL

        Agents[nameof(agent)] = RL_policy
    end
    
    #-------------------------------------------------------------------------------
    # Finalising the Control

    num_policies = num_RL_sources + convert(Int64, num_clas_sources > 0)

    if num_clas_sources > 0

        polc = Dict()

        polc["policy"] = Animo
        polc["state_ids"] = Animo.policy.state_ids
        polc["action_ids"] = Animo.policy.action_ids

        Agents[nameof(Animo)] = polc
    end
    
    Multi_Agent = MultiAgentGridController(Agents, env.action_ids)

    #_______________________________________________________________________________
    # Logging

    if env.verbosity > 1

        @info "Time simulation run time: $((env.maxsteps - 1)*env.ts) [s] ~ $(Int(env.maxsteps)) steps"
        @info "Number of episodes: $(num_episodes)"
    end
    
    #_______________________________________________________________________________
    # Running the time simulation
    
    RLBase.run(Multi_Agent, env, StopAfterEpisode(num_episodes), hook);    

    if return_Agents

        return Multi_Agent
    else
        
        return nothing
    end
end


#= 
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

    function run(dare_setup, no_episodes)

        RLBase.run(dare_setup.control_agent, dare_setup.env, StopAfterEpisode(no_episodes), dare_setup.hook)

        return nothing
    end 
=#