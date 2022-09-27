using DrWatson
@quickactivate "dare"

using ReinforcementLearning
using PlotlyJS

include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent_ddpg.jl"))
include(srcdir("data_hook.jl"))
include(srcdir("Classical_Control.jl"))
include(srcdir("Power_System_Theory.jl"))
include(srcdir("MultiAgentGridController.jl"))

include(srcdir("Classical_Control_Plots.jl"))

function reference(t)
    
    u = [sqrt(2)*230 * cos.(2*pi*50*t .- 2/3*pi*(i-1)) for i = 1:3]
    #return vcat(u,u)  # to control 2 sources
    return u
end

function reward(env, name = nothing)
    r = 0.0
    
    if !isnothing(name)
        if name == "agent"
            u_l1_index = findfirst(x -> x == "source1_u_C_a", env.state_ids)
            u_l2_index = findfirst(x -> x == "source1_u_C_b", env.state_ids)
            u_l3_index = findfirst(x -> x == "source1_u_C_c", env.state_ids)
        else
            u_l1_index = findfirst(x -> x == "source2_u_C_a", env.state_ids)
            u_l2_index = findfirst(x -> x == "source2_u_C_b", env.state_ids)
            u_l3_index = findfirst(x -> x == "source2_u_C_c", env.state_ids)
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

#_______________________________________________________________________________
# Parameters - Time simulation
Timestep = 100 #time step in μs ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_final = 3.4 #time in seconds, total simulation run time

ts = Timestep*1e-6
t = 0:ts:t_final # time

fs = 1/ts # Hz, Sampling frequency of controller ~ 15 kHz < fs < 50kHz

#_______________________________________________________________________________
# State space representation

#-------------------------------------------------------------------------------
# Connectivity Matrix

CM = [ 0. 0. 1.
     0. 0. 2
     -1. -2. 0.]

#-------------------------------------------------------------------------------
# Cables

cable_list = []

# Network Cable Impedances
l = 0.50 # length in km
cable = Dict()
cable["R"] = 0.208*l # Ω, line resistance 0.722#
cable["L"] = 0.00025*l # H, line inductance 0.264e-3#
cable["C"] = 0.4e-6*l # 0.4e-6#

push!(cable_list, cable, cable)

# Sources

source_list = []
source = Dict()

source["pwr"] = 200e3
source["vdc"] = 800
source["fltr"] = "LC"
Lf, Cf, _ = Filter_Design(source["pwr"], fs)
source["R1"] = 0.4
source["R_C"] = 0.0006
source["L1"] = Lf
source["C"] = Cf

push!(source_list, source)

source = Dict()

source["pwr"] = 200e3
source["vdc"] = 800
source["fltr"] = "LC"
Lf, Cf, _ = Filter_Design(source["pwr"], fs)
source["R1"] = 0.4
source["R_C"] = 0.0006
source["L1"] = Lf
source["C"] = Cf

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
parameters["grid"] = Dict("fs" => fs, "phase" => 3, "v_rms" => 230)

# Define the environment

env = SimEnv(reward_function = reward, featurize = featurize, 
v_dc = 800, ts = ts, use_gpu = false, CM = CM, num_sources = 2, num_loads = 1, parameters = parameters,
maxsteps = length(t), action_delay = 0)

#_______________________________________________________________________________
# Setting up the Classical Sources

state_ids = get_state_ids(env.nc)
action_ids = get_action_ids(env.nc)

state_ids_classic = filter(x -> (split(x, "_")[1] == "source1" || split(x, "_")[1] == "source2"), state_ids)
action_ids_classic = filter(x -> (split(x, "_")[1] == "source1" || split(x, "_")[1] == "source2"), action_ids)

Animo = NamedPolicy("classic", Classical_Policy(action_space = Space([-1.0..1.0 for i in 1:length(action_ids_classic)]), t_final = t_final, 
fs = fs, num_sources = 2, state_ids = state_ids_classic, action_ids = action_ids_classic))

#= Modes:
    1 -> "Swing Mode" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "Voltage Control Mode" - voltage source with controller dynamics
    3 -> "PQ Control Mode" - grid following controllable source/load
    4 -> "Droop Control Mode" - simple grid forming with power balancing
    
        "Synchronverter Modes" - grid forming with power balancing via virtual motor (advanced controllable source/load)
    5 -> "Full-Synchronverter Mode" - droop control on real and imaginary powers
    6 -> "Self-Synchronverter Mode" - active control on real and imaginary powers
    7 -> "Infinite-Synchronverter Mode" - droop characteristic on real power, and active control on voltage
    8 -> "Semi-Synchronverter Mode" - droop characteristic on imaginary power, and active control on real power
    9 -> "Null-Synchronverter Mode" - active control on real power and voltage
=#

Source_Initialiser(env, Animo, [6 5])

Animo.policy.Source.τv = 0.02 # time constant of the voltage loop # 0.02
Animo.policy.Source.τf = 0.002 # time constant of the frequency loop # 0.002

nm_src = 1 # changing the power set points of the source
Animo.policy.Source.pq0_set[nm_src, 1] = 160e3 # W, Real Power
Animo.policy.Source.pq0_set[nm_src, 2] = -125e3 # VAi, Imaginary Power

Animo.policy.Source.V_pu_set[nm_src, 1] = 1.02
Animo.policy.Source.V_δ_set[nm_src, 1] = 0 #-90*π/180

nm_src = 2 # changing the power set points of the source
Animo.policy.Source.pq0_set[nm_src, 1] = 35e3 # W, Real Power
Animo.policy.Source.pq0_set[nm_src, 2] = 55e3 # VAi, Imaginary Power

Animo.policy.Source.V_pu_set[nm_src, 1] = 0.98
Animo.policy.Source.V_δ_set[nm_src, 1] = 0 #-90*π/180
ma_agents = Dict(nameof(Animo) => Dict("policy" => Animo,
                            "state_ids" => state_ids_classic,
                            "action_ids" => action_ids_classic))
                            
ma = MultiAgentGridController(ma_agents, action_ids)

agentname = "agent"

#_______________________________________________________________________________
#%% Starting time simulation

plt_state_ids = ["source1_u_C_a", "source1_u_C_b", "source1_u_C_c", "source2_u_C_a", "source2_u_C_b", "source2_u_C_c", "source1_i_L1_a", "source2_i_L1_a"]
plt_action_ids = []#"u_v1_a", "u_v1_b", "u_v1_c"]
hook = DataHook(collect_state_ids = plt_state_ids, collect_action_ids = plt_action_ids, save_best_NNA = true, collect_reference = false, plot_rewards=false)

run(ma, env, StopAfterEpisode(1), hook);

#_______________________________________________________________________________
# Plotting

plot_hook_results(; hook = hook, actions_to_plot = [] ,plot_reward = false, plot_reference = false, episode = 1)

#= Plot_Vrms(5, 5000, Animo.policy.Source, num_source = 1)
Plot_Vrms(5, 5000, Animo.policy.Source, num_source = 2)

Plot_Real_Imag_Active_Reactive(0, 5000, Animo.policy.Source, num_source = 1)
Plot_Real_Imag_Active_Reactive(0, 5000, Animo.policy.Source, num_source = 2) =#

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
