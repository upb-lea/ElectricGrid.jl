#using Dare
#using ReinforcementLearning
#using PlotlyJS

using DrWatson
@quickactivate "dare"

using PlotlyJS
using ReinforcementLearning

include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent_ddpg.jl"))
include(srcdir("data_hook.jl"))
include(srcdir("Classical_Control.jl"))
include(srcdir("Power_System_Theory.jl"))
include(srcdir("MultiAgentGridController.jl"))


function reference(t)
    
    u = [sqrt(2)*30 * cos.(2*pi*50*t .- 2/3*pi*(i-1)) for i = 1:3]
    #u = [230, 230, 230]
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
t_final = 0.1 #time in seconds, total simulation run time
ts = 1e-4
t = 0:ts:t_final # time

fs = 1/ts # Hz, Sampling frequency of controller ~ 15 kHz < fs < 50kHz


#=CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.] =#
 CM = [0. 1.
   -1. 0.] 
#-------------------------------------------------------------------------------
#= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "Voltage Control" - voltage source with controller dynamics

    3 -> "PQ Control" - grid following controllable source/load (active and reactive Power)
    4 -> "PV Control" - grid following controllable source (active power and voltage magnitude)

    5 -> "Droop Control" - simple grid forming with power balancing (i.e. load partitioning)
    6 -> "Full-Synchronverter" - droop control on real and imaginary powers
    7 -> "Semi-Synchronverter" - droop characteristic on real power, and active control on voltage
=#
#-------------------------------------------------------------------------------


#R_load, L_load, _, Z = Parallel_Load_Impedance(1e3, 0.6, 230)

parameters = Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}("fltr"=>"LC", "pwr"=>10000e3, "control_type" =>"classic", "source_type"=>"ideal", "mode" => 8,  "R1"=>1.1e-3, "L1"=>1e-3, "C"=>1e-3, "R_C"=>7e-3, "vdc"=>800, "v_limit"=>10000, "i_limit"=>10000)
                    ],
    "load"   => Any[
                    Dict{Any, Any}("impedance"=>"RLC", "R"=>100, "L"=>1e-2, "C"=>1e-2, "pf"=>0.8, "v_limit"=>10000, "i_limit"=>10000) 
                    ],
    "cable" => Any[
                    Dict{Any, Any}("len"=>1, "R"=>1e-3, "L"=>1e-4, "C"=>1e-4, "i_limit"=>10000),
                    ],
    "grid"   => Dict{Any, Any}("fs"=>50.0, "phase"=>3, "v_rms"=>230, "f_grid" => 50, "ramp_end"=>0.0)
                
)


env = SimEnv(ts = ts, use_gpu = false, CM = CM, num_sources = 1, num_loads = 1, parameters = parameters, maxsteps = length(t), action_delay = 1)

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
plt_state_ids = ["source1_v_C_filt_a", "source1_i_L1_a"]               
plt_action_ids = ["source1_u_a"]#"source1_u_a", "u_v1_b", "u_v1_c"]
hook = DataHook(collect_state_ids = plt_state_ids, collect_action_ids = plt_action_ids,  collect_sources = [],
collect_cables = [], collect_vrms_ids = [1], collect_irms_ids = [], collect_pq_ids = [], collect_vdq_ids = [], collect_idq_ids = [],
save_best_NNA = false, collect_reference = false, plot_rewards = false)

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

    plot_hook_results(; hook = hook, states_to_plot = plt_state_ids, actions_to_plot = plt_action_ids, episode = eps, 
    pq_to_plot = [], vrms_to_plot = [], irms_to_plot = [], vdq_to_plot = [], idq_to_plot = [])
end

#===
for eps in 1:num_eps

    plot_hook_results(; hook = hook, states_to_plot = [], actions_to_plot = [], episode = eps, 
    pq_to_plot = [1 2], vrms_to_plot = [1 2], irms_to_plot = [], vdq_to_plot = [], idq_to_plot = [])
end
===#
print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
