using Dare;

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Parameters - Time simulation
Timestep = 100 #time step in μs ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_final = 0.5 #time in seconds, total simulation run time

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

#= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "Voltage Control" - voltage source with controller dynamics

    3 -> "PQ Control" - grid following controllable source/load

    4 -> "Droop Control" - simple grid forming with power balancing
    5 -> "Full-Synchronverter" - droop control on real and imaginary powers
    6 -> "Semi-Synchronverter" - droop characteristic on real power, and active control on voltage
=#

source = Dict()

source_list = []

source["pwr"] = 200e3
source["vdc"] = 800
source["fltr"] = "LC"
source["p_set"] = 50e3
source["q_set"] = 10e3
source["v_pu_set"] = 1.05
source["mode"] = 6
source["control_type"] = "classic"

push!(source_list, source)

source = Dict()

source["pwr"] = 100e3
source["vdc"] = 800
source["fltr"] = "LC"
source["p_set"] = 50e3
source["q_set"] = 10e3
source["v_pu_set"] = 1.0
source["mode"] = 3
source["control_type"] = "classic"

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

#setup = create_setup(parameters = parameters, CM = CM)

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
#%% Starting time simulation

plt_state_ids = []#"source1_v_C_a", "source2_v_C_a", "source1_i_L1_a", "source2_i_L1_a"]
plt_action_ids = []#"u_v1_a", "u_v1_b", "u_v1_c"]
hook = DataHook(collect_state_ids = plt_state_ids, collect_action_ids = plt_action_ids, 
collect_vrms_ids = [1 2], collect_irms_ids = [1 2], collect_pq_ids = [1 2],
save_best_NNA = true, collect_reference = false, plot_rewards = false)

RLBase.run(ma, env, StopAfterEpisode(1), hook);

#_______________________________________________________________________________
# Plotting

plot_hook_results(; hook = hook, actions_to_plot = [], episode = 1, 
pq_to_plot = [1 2], vrms_to_plot = [1 2], irms_to_plot = [1 2])

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
