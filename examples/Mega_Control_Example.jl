using Dare
using ReinforcementLearning

function reference(t)
    
    u = [sqrt(2)*230 * cos.(2*pi*50*t .- 2/3*pi*(i-1)) for i = 1:3]
    #return vcat(u,u)  # to control 2 sources
    return u
end

function reward(env, name = nothing)
    r = 0.0
    
    if !isnothing(name)
        if name == "agent"
            u_l1_index = findfirst(x -> x == "source1_v_C_cables_a", env.state_ids)
            u_l2_index = findfirst(x -> x == "source1_v_C_cables_b", env.state_ids)
            u_l3_index = findfirst(x -> x == "source1_v_C_cables_c", env.state_ids)
        else
            u_l1_index = findfirst(x -> x == "source1_v_C_cables_a", env.state_ids)
            u_l2_index = findfirst(x -> x == "source1_v_C_cables_b", env.state_ids)
            u_l3_index = findfirst(x -> x == "source1_v_C_cables_c", env.state_ids)
        end

        u_l1 = env.state[u_l1_index]
        u_l2 = env.state[u_l2_index]
        u_l3 = env.state[u_l3_index]

        u = [u_l1, u_l2, u_l3]
        refs = reference(env.t)

        r = -(sum(abs.(refs/600 - u)/3))
    end

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

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Parameters - Time simulation
Timestep = 100 #time step in μs ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_final = 1.0 #time in seconds, total simulation run time

ts = Timestep*1e-6
t = 0:ts:t_final # time

fs = 1/ts # Hz, Sampling frequency of controller ~ 15 kHz < fs < 50kHz

#_______________________________________________________________________________
# State space representation

#-------------------------------------------------------------------------------
# Connectivity Matrix

#= CM = [ 0. 0. 0. 1.
        0. 0. 0. 2.
        0. 0. 0. 3.
        -1. -2. -3. 0.] =#

CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.]

#= CM = [0. 1.
   -1. 0.] =#

#-------------------------------------------------------------------------------
# Cables

cable_list = []

# Network Cable Impedances
l = 1# length in km
cable = Dict()
cable["R"] = 0.208*l # Ω, line resistance 0.722#
cable["L"] = 0.00025*l # H, line inductance 0.264e-3#
cable["C"] = 0.4e-3*l # 0.4e-6#

#push!(cable_list, cable, cable, cable)

push!(cable_list, cable, cable)

#push!(cable_list, cable)

#-------------------------------------------------------------------------------
# Sources

#= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "Voltage Control" - voltage source with controller dynamics

    3 -> "PQ Control" - grid following controllable source/load (active and reactive Power)
    4 -> "PV Control" - grid following controllable source (active power and voltage magnitude)

    5 -> "Droop Control" - simple grid forming with power balancing (i.e. load partitioning)
    6 -> "Full-Synchronverter" - droop control on real and imaginary powers
    7 -> "Semi-Synchronverter" - droop characteristic on real power, and active control on voltage
=#

source = Dict()

source_list = []

source["pwr"] = 200e3 #VA
source["vdc"] = 800 #V
source["fltr"] = "LCL"
source["p_set"] = 150e3 #Watt
source["q_set"] = 10e3 #VAr
source["v_pu_set"] = 1.0 #p.u.
source["v_δ_set"] = 0 # degrees
source["mode"] = 7
source["control_type"] = "classic"
source["v_rip"] = 0.01537
source["i_rip"] = 0.15
source["τv"] = 0.002
source["τf"] = 0.002
source["std_asy"] = 50 # asymptotic standard deviation
#source["κ"] = 3 # mean reversion parameter
source["σ"] = 1000.0 # Brownian motion scale i.e. ∝ diffusion, volatility parameter
source["γ"] = 0 # asymptotoic mean
source["X₀"] = 0 # initial value
source["Δt"] = ts # time step

#= source["L1"] = 0.002
source["R1"] = 0.04
source["L2"] = 0.002
source["R2"] = 0.05
source["R_C"] = 0.09
source["C"] = 0.003 =#

push!(source_list, source)

source = Dict()

source["pwr"] = 100e3
source["vdc"] = 800
source["fltr"] = "LC"
source["p_set"] = 50e3
source["q_set"] = 10e3
source["v_pu_set"] = 1.0
source["v_δ_set"] = 0 # degrees
source["mode"] = 7
source["control_type"] = "classic"
source["v_rip"] = 0.01537
source["i_rip"] = 0.15
source["τv"] = 0.002
source["τf"] = 0.002
source["std_asy"] = 50 # asymptotic standard deviation
source["σ"] = 1000.0 # Brownian motion scale i.e. ∝ diffusion parameter
source["γ"] = 0 # asymptotoic mean
source["X₀"] = 0 # initial values
source["Δt"] = 4

#= source["L1"] = 1e-3
source["R1"] = 0.05 =#

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

R1_load, L_load, _, _ = Parallel_Load_Impedance(100e3, 0.6, 230)
#R2_load, C_load, _, _ = Parallel_Load_Impedance(150e3, -0.8, 230)

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

#setup = create_setup(parameters)

# Define the environment

num_sources = length(source_list)
num_loads = length(load_list)

env = SimEnv(reward_function = reward, featurize = featurize, 
ts = ts, use_gpu = false, CM = CM, num_sources = num_sources, num_loads = num_loads, 
parameters = parameters, maxsteps = length(t), action_delay = 1)

state_ids = get_state_ids(env.nc)
action_ids = get_action_ids(env.nc)

#_______________________________________________________________________________
# Setting up the Classical Sources

Animo = NamedPolicy("classic", Classical_Policy(env))

state_ids_classic = Animo.policy.state_ids
action_ids_classic = Animo.policy.action_ids

Multi_Agents = Dict()
Multi_Agent_list = []

polc = Dict()

polc["policy"] = Animo
polc["state_ids"] = state_ids_classic
polc["action_ids"] = action_ids_classic

Multi_Agents[nameof(Animo)] = polc

ma = MultiAgentGridController(Multi_Agents, action_ids)

agentname = "agent"

#_______________________________________________________________________________
#%% Setting up data hooks

#= plt_state_ids = ["source1_v_C_a", "source1_v_C_b", "source1_v_C_c",
                "source2_v_C_a", "source2_v_C_b", "source2_v_C_c", 
                "source1_i_L1_a", "source1_i_L1_b", "source1_i_L1_c", 
                "source2_i_L1_a", "source2_i_L1_b", "source2_i_L1_c"] =#
plt_state_ids = []               
plt_action_ids = []#"source1_u_a", "u_v1_b", "u_v1_c"]
hook = DataHook(collect_state_ids = plt_state_ids, collect_action_ids = plt_action_ids,  collect_sources = [1 2],
collect_cables = [1], collect_vrms_ids = [1 2], collect_irms_ids = [1 2], collect_pq_ids = [1 2], collect_vdq_ids = [1 2], collect_idq_ids = [1 2],
save_best_NNA = false, collect_reference = false, plot_rewards = false, collect_debug = [3 4])

#_______________________________________________________________________________
# Starting time simulation
num_eps = 1
RLBase.run(ma, env, StopAfterEpisode(num_eps), hook);

#_______________________________________________________________________________
# Plotting
#= st = ["source1_i_L1_a", "source1_i_L1_b", "source1_i_L1_c", 
"source1_v_C_filt_a", "source1_v_C_filt_b", "source1_v_C_filt_c",
"source1_i_C_filt_a", "source1_i_C_filt_b", "source1_i_C_filt_c",
"cable1_i_L_a", "cable1_i_L_b", "cable1_i_L_c",
"source1_v_C_cables_a", "source1_v_C_cables_b", "source1_v_C_cables_c",
"source1_i_C_cables_a", "source1_i_C_cables_b", "source1_i_C_cables_c"] =#
#= plot_hook_results(; hook = hook, states_to_plot = ["source1_i_L2_a", "source1_v_C_filt_a", "source1_v_C_cables_a" ], actions_to_plot = [], episode = 1, 
pq_to_plot = [], vrms_to_plot = [], irms_to_plot = [], vdq_to_plot = []) =#

for eps in 1:num_eps

    plot_hook_results(; hook = hook, states_to_plot = [], actions_to_plot = [], episode = eps, 
    pq_to_plot = [1 2], vrms_to_plot = [1 2], irms_to_plot = [], vdq_to_plot = [], idq_to_plot = [])
end

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
