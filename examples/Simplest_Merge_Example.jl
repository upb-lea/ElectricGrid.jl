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

CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.]


#-------------------------------------------------------------------------------
# Parameters

#= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "PQ" - grid following controllable source/load (active and reactive Power)
    3 -> "Droop" - simple grid forming with power balancing
    4 -> "Synchronverter" - enhanced droop control
=#

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => 200e3, "control_type" => "RL"),
                        Dict{Any, Any}("pwr" => 200e3),
                        ],
         "load"   => Any[
                        Dict{Any, Any}("impedance" => "RL", "R" => 2.64, "L" => 0.006),
                        ],
        "cable"   => Any[
                        Dict{Any, Any}("R" => 1e-3, "L" => 1e-4, "C" => 1e-4, "i_limit" => 10e4,),
                        ],
        "grid" => Dict{Any, Any}("ramp_end" => 0.0)
    )
#_______________________________________________________________________________
# Defining the environment
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

hook = DataHook(collect_vrms_ids = [1], 
                collect_irms_ids = [1], 
                collect_pq_ids   = [1], 
                collect_freq     = [1],
                collect_sources  = [1])

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
                    states_to_plot  = [], 
                    actions_to_plot = [],  
                    p_to_plot       = [1], 
                    q_to_plot       = [1], 
                    vrms_to_plot    = [1], 
                    irms_to_plot    = [],
                    freq_to_plot    = [])

print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")