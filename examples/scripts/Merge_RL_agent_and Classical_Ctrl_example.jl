using DrWatson
@quickactivate "JEG"

using ReinforcementLearning
using PlotlyJS

include(srcdir("node_constructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent_ddpg.jl"))
include(srcdir("data_hook.jl"))
include(srcdir("Classical_Control.jl"))
include(srcdir("power_system_theory.jl"))
include(srcdir("MultiAgentGridController.jl"))


agentname = "agent"
classicname = "classic"

function reference(t)
    
    u = [sqrt(2)*230 * cos.(2*pi*50*t .- 2/3*pi*(i-1)) for i = 1:3]
    #return vcat(u,u)  # to control 2 sources
    return u
end

function reward(env, name = nothing)
    r = 0.0
    
    if !isnothing(name)
        if name == "agent"
            u_l1_index = findfirst(x -> x == "source1_v_C_a", env.state_ids)
            u_l2_index = findfirst(x -> x == "source1_v_C_b", env.state_ids)
            u_l3_index = findfirst(x -> x == "source1_v_C_c", env.state_ids)
        else
            u_l1_index = findfirst(x -> x == "source2_v_C_a", env.state_ids)
            u_l2_index = findfirst(x -> x == "source2_v_C_b", env.state_ids)
            u_l3_index = findfirst(x -> x == "source2_v_C_c", env.state_ids)
        end

        u_l1 = env.state[u_l1_index]
        u_l2 = env.state[u_l2_index]
        u_l3 = env.state[u_l3_index]

        u = [u_l1, u_l2, u_l3]
        refs = reference(env.t)

        r = -(sum(abs.(refs/600 - u)/3)) * (1 - 0.99)
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

env_cuda = false
agent_cuda = false

#CM = [ 0. 1.
#     -1.  0.]

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
Lf, Cf, _ = FilterDesign(source["pwr"], fs)
source["R1"] = 0.4
source["R_C"] = 0.0006
source["L1"] = Lf
source["C"] = Cf
source["v_limit"]= 220 
source["i_limit"]= 1000

push!(source_list, source)

source = Dict()

source["pwr"] = 100e3
source["vdc"] = 600
source["fltr"] = "LC"
Lf, Cf, _ = FilterDesign(source["pwr"], fs)
source["R1"] = 0.4
source["R_C"] = 0.0006
source["L1"] = Lf
source["C"] = Cf
source["v_limit"]= 230 
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
parameters["grid"] = Dict("fs" => fs, "phase" => 3, "v_rms" => 230)

# Define the environment
V_source = 800

env = SimEnv(reward_function = reward, featurize = featurize, ts = ts, use_gpu = env_cuda, CM = CM, num_sources = 2, num_loads = 1, parameters = parameters,
maxsteps = 1000, action_delay = 0)

state_ids = get_state_ids(env.nc)
action_ids = get_action_ids(env.nc)

state_ids_agent = filter(x -> split(x, "_")[1] == "source1", state_ids)
action_ids_agent = filter(x -> split(x, "_")[1] == "source1", action_ids)
state_ids_classic = filter(x -> split(x, "_")[1] == "source2", state_ids)
action_ids_classic = filter(x -> split(x, "_")[1] == "source2", action_ids)

function RLBase.action_space(env::SimEnv, name::String)
    if name == "agent"
        return Space(fill(-1.0..1.0, size(action_ids_agent)))
    else
        return Space(fill(-1.0..1.0, size(action_ids_classic)))
    end
end

na = length(env.action_space)
agent = create_agent_ddpg(na = length(action_ids_agent), ns = length(state(env,agentname)), use_gpu = agent_cuda)
agent = Agent(policy = NamedPolicy(agentname, agent.policy), trajectory = agent.trajectory)

#agent = NamedPolicy(agentname, Classical_Policy(action_space = Space([-1.0..1.0 for i in 1:length(action_ids_agent)]), t_final = ts*1001, fs = fs, num_sources = 1, state_ids = state_ids_agent, action_ids = action_ids_agent))

#_______________________________________________________________________________
# Setting up the Sources

Animo = NamedPolicy(classicname, Classical_Policy(action_space = Space([-1.0..1.0 for i in 1:length(action_ids_classic)]), t_final = ts*1001, 
fs = fs, num_sources = 1, state_ids = state_ids_classic, action_ids = action_ids_classic))

#=
    1 -> "Swing Mode" - voltage source without dynamics
    2 -> "Voltage Control Mode" - voltage source with controller dynamics
    3 -> "PQ Control Mode" - grid following controllable load/source
    4 -> "Droop Control Mode" - grid forming with power balancing
    5 -> "Synchronverter Mode" - grid forming with power balancing via virtual motor
    6 -> "Self-Synchronverter Mode" - grid forming with power balancing via virtual motor
=#

Source_Initialiser(env, Animo, [2])

#TODO: entfernen!!!
#Animo = create_agent_ddpg(na = length(action_ids_classic), ns = length(state(env,classicname)), use_gpu = agent_cuda)

#Animo = Agent(policy = NamedPolicy(classicname, Animo.policy), trajectory = Animo.trajectory)

ma_agents = Dict(nameof(agent) => Dict("policy" => agent,
                            "state_ids" => state_ids_agent,
                            "action_ids" => action_ids_agent),
                nameof(Animo) => Dict("policy" => Animo,
                            "state_ids" => state_ids_classic,
                            "action_ids" => action_ids_classic))
                            
ma = MultiAgentGridController(ma_agents, action_ids)

plt_state_ids = ["source1_v_C_a", "source1_v_C_b", "source1_v_C_c", "source2_v_C_a", "source2_v_C_b", "source2_v_C_c", "source1_i_L1_a", "source2_i_L1_a"]
plt_action_ids = []#"u_v1_a", "u_v1_b", "u_v1_c"]
hook = data_hook(collect_state_ids = plt_state_ids, collect_action_ids = plt_action_ids, save_best_NNA = true, collect_reference = true, plot_rewards=true, collect_vdc_idx = [2])

ma["agent"]["policy"].policy.policy.act_noise = 0.1
run(ma, env, StopAfterEpisode(1), hook);


###############################
# Plotting
#plot_hook_results(; hook = hook, actions_to_plot = [] ,plot_reward = false, plot_reference = true, episode = 3000)

plot_best_results(;agent = ma, env = env, hook = hook, states_to_plot = plt_state_ids, plot_reward = false, plot_reference = true, use_best = false)