using Dare
using ReinforcementLearning
using IntervalSets

print("\n...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Network Configuration 

#-------------------------------------------------------------------------------
# Time simulation

Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_end    = 0.2     # total run time, seconds
num_eps  = 100       # number of episodes to run

#-------------------------------------------------------------------------------
# Connectivity Matrix
#=
CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.]
=#

CM = [0. 1.
      -1. 0.]

#-------------------------------------------------------------------------------
# Parameters

#= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "PQ" - grid following controllable source/load (active and reactive Power)
    3 -> "Droop" - simple grid forming with power balancing
    4 -> "Synchronverter" - enhanced droop control
=#
# Filter
R_1 = 1.1e-3;
L_1 = 1e-4;
R_c = 7e-3;
C_1 = 1e-5;

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => 200e3, "control_type" => "RL", "fltr" => "L", "L1" => L_1, "R1" => R_1, "i_limit"=>100000),
                        #Dict{Any, Any}("pwr" => 200e3, "fltr" => "LC", "control_type" => "classic", "mode" => 1, "R1"=>R_1, "L1"=>L_1, "C"=>C_1, "R_C"=>R_c, "vdc"=>800, "v_limit"=>10000, "i_limit"=>10e8),
                        ],
         "load"   => Any[
                        Dict{Any, Any}("impedance" => "RL", "R" => 2.64, "L" => 0.006, "v_limit"=>10000, "i_limit"=>10e8),
                        ],
        "cable"   => Any[
                        Dict{Any, Any}("R" => 1e-3, "L" => 1e-4, "C" => 1e-4, "i_limit" => 10e8,),
                        ],
        "grid" => Dict{Any, Any}("ramp_end" => 0.0)
    )
#_______________________________________________________________________________
# Defining the environment

function reference(t)
    θ = 2*pi*50*t
    θph = [θ; θ - 120π/180; θ + 120π/180]
    #i = [10 * cos.(2*pi*50*t .- 2/3*pi*(i-1)) for i = 1:3]
    return 1# * cos.(θph)
end

function reward(env)
    
    index_1 = findfirst(x -> x == "source1_i_L1_a", env.state_ids)
    index_2 = findfirst(x -> x == "source1_i_L1_b", env.state_ids)
    index_3 = findfirst(x -> x == "source1_i_L1_c", env.state_ids)
    
    u_l1 = env.state[index_1]
    u_l2 = env.state[index_2]
    u_l3 = env.state[index_3]

    # better:
    #u = env.state[findall(x -> x in env.state_ids, env.state_ids_RL)]
    # problem, u_cable is state as well

    u = [u_l1, u_l2, u_l3]
    any(abs.(env.x./env.norm_array) .> 1)
    if any(u.>1)
        return -1
    else

        refs = reference(env.t)
        norm_ref = env.nc.parameters["source"][1]["i_limit"]
        r = 1-(sum(abs.(refs/norm_ref - u)/2))/3  # TODO: replace by correct entry of norm array
        return r
    end
end


function featurize(x0 = nothing, t0 = nothing; env = nothing, name = nothing)
        if !isnothing(name)
            state = env.state
            if name == "agent"
                state = state[findall(x -> x in env.state_ids_RL, env.state_ids)]
                #state = vcat(state, reference(env.t)/600)
            #else
            #    global state_ids_classic
            #    global state_ids
            #    state = env.x[findall(x -> x in state_ids_classic, state_ids)]
            end
        elseif isnothing(env)
            return x0
        else
            return env.state
        end
        return state
end


env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 2, featurize = featurize)

#_______________________________________________________________________________
# Setting up data hooks

hook = DataHook(collect_vrms_ids = [], 
                collect_irms_ids = [], 
                collect_pq_ids   = [], 
                collect_freq     = [],
                collect_sources  = [1, 2])

#_______________________________________________________________________________
# Running the Time Simulation

function RLBase.action_space(env::SimEnv, name::String)
        if name == "agent"
                return Space(fill(-1.0..1.0, 3))#size(env.action_ids_agent)))
        end
end

ma = Power_System_Dynamics(env, hook, num_episodes = num_eps, return_Agents = true)

#_______________________________________________________________________________
# Plotting

plot_hook_results(hook = hook, 
                    episode = num_eps,
                    #states_to_plot  = ["source1_i_L1_a", "source2_i_L1_a", "source2_v_C_filt_a"],
                    states_to_plot  = env.state_ids,  
                    actions_to_plot = env.action_ids,  
                    p_to_plot       = [], 
                    q_to_plot       = [], 
                    vrms_to_plot    = [], 
                    irms_to_plot    = [],
                    freq_to_plot    = [])

print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")