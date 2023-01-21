#using Dare

using DrWatson
@quickactivate "dare"

include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent_ddpg.jl"))
include(srcdir("data_hook.jl"))
include(srcdir("Dare_Wrapper.jl"))
include(srcdir("Classical_Control.jl"))
include(srcdir("Power_System_Theory.jl"))
include(srcdir("MultiAgentGridController.jl"))

print("\n...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Network Parameters 

#-------------------------------------------------------------------------------
# Time simulation

Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_end    = 0.2     # total run time, seconds

#-------------------------------------------------------------------------------
# Connectivity Matrix

CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.]

#-------------------------------------------------------------------------------
# Sources

source_list = []

source = Dict()

source["pwr"]      = 200e3  # Rated Apparent Power, VA

push!(source_list, source)

source = Dict()

source["pwr"]      = 100e3  # Rated Apparent Power, VA

push!(source_list, source)

#-------------------------------------------------------------------------------
# Loads

load_list = []
load = Dict()

load["impedance"] = "RL"
load["R"] = 2.64
load["L"] = 0.006

push!(load_list, load)

#-------------------------------------------------------------------------------
# Amalgamation

parameters = Dict()

parameters["source"] = source_list
parameters["cable"] = cable_list
parameters["load"] = load_list

#= parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr"=>10000e3)
                        ],
        "load"   => Any[
                        Dict{Any, Any}("impedance"=>"RLC", "R"=>100, "L"=>1e-2, "C"=>1e-2, "pf"=>0.8, "v_limit"=>10000, "i_limit"=>10000) 
                        ],
        "grid"   => Dict{Any, Any}("fs"=>50.0, "phase"=>3, "v_rms"=>230, "f_grid" => 50, "ramp_end"=>0.0)
    ) =#

#bugs ensue when one of these is removed, but not both
parameters["grid"] = Dict("v_rms" => 230, "Δfmax" => 0.005) 

# Define the environment

env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end)

#_______________________________________________________________________________
# Setting up data hooks

hook = DataHook(collect_vrms_ids = [1 2], 
                collect_irms_ids = [1 2], 
                collect_pq_ids   = [1 2],
                collect_freq     = [1 2])

#_______________________________________________________________________________
# Running the Time Simulation

Power_System_Dynamics(env, hook, num_episodes = 1)

#_______________________________________________________________________________
# Plotting

plot_hook_results(hook = hook, 
                    states_to_plot  = [], 
                    actions_to_plot = [],  
                    p_to_plot       = [1 2], 
                    q_to_plot       = [], 
                    vrms_to_plot    = [1 2], 
                    irms_to_plot    = [1 2],
                    freq_to_plot    = [1 2])

print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")
