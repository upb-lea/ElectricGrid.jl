using Dare

#= using DrWatson
@quickactivate "dare"

include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("agent_ddpg.jl"))
include(srcdir("data_hook.jl"))
include(srcdir("Dare_Wrapper.jl"))
include(srcdir("Classical_Control.jl"))
include(srcdir("Power_System_Theory.jl"))
include(srcdir("MultiAgentGridController.jl")) =#

print("\n...........o0o----ooo0§0ooo~~~  START  ~~~ooo0§0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Network Parameters 

#-------------------------------------------------------------------------------
# Time simulation

Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_end    = 0.2     # total run time, seconds
num_eps  = 4       # number of episodes to run

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
# Cable Impedances

cable_list = []

cable = Dict()
cable["R"]       = 0.208   # Ω, line resistance
cable["L"]       = 0.00025 # H, line inductance
cable["C"]       = 0.4e-3  # F, line capacitance
cable["i_limit"] = 10e12   # A, line current limit

#push!(cable_list, cable, cable, cable)

push!(cable_list, cable, cable)

#push!(cable_list, cable)

#-------------------------------------------------------------------------------
# Sources

#= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "PQ" - grid following controllable source/load (active and reactive Power)
    3 -> "Droop" - simple grid forming with power balancing
    4 -> "Synchronverter" - enhanced droop control
=#

source_list = []

source = Dict()

source["mode"]     = 4

source["fltr"]     = "LC"  # Filter type

source["pwr"]      = 200e3  # Rated Apparent Power, VA
source["p_set"]    = 50e3   # Real Power Set Point, Watt
source["q_set"]    = 10e3   # Imaginary Power Set Point, VAi

source["v_pu_set"] = 1.00   # Voltage Set Point, p.u.
source["v_δ_set"]  = 0      # Voltage Angle, degrees

source["std_asy"]  = 50e3   # Asymptotic Standard Deviation
source["σ"]        = 0.0    # Brownian motion scale i.e. ∝ diffusion, volatility parameter
source["Δt"]       = 0.02   # Time Step, seconds
source["X₀"]       = 0      # Initial Process Values, Watt
source["k"]        = 0      # Interpolation degree

source["τv"]       = 0.002  # Time constant of the voltage loop, seconds
source["τf"]       = 0.002  # Time constant of the frequency loop, seconds

source["Observer"] = true   # Discrete Luenberger Observer

push!(source_list, source)

source = Dict()

source["mode"]     = 2

source["fltr"]     = "L"   # Filter type

source["pwr"]      = 100e3  # Rated Apparent Power, VA
source["p_set"]    = 50e3   # Real Power Set Point, Watt
source["q_set"]    = 10e3   # Imaginary Power Set Point, VAi

source["v_pu_set"] = 1.00   # Voltage Set Point, p.u.
source["v_δ_set"]  = 0      # Voltage Angle, degrees

source["std_asy"]  = 50e3   # Asymptotic Standard Deviation
source["σ"]        = 50e3   # Brownian motion scale i.e. ∝ diffusion, volatility parameter
source["Δt"]       = 0.01   # Time Step, seconds
source["X₀"]       = 0      # Initial Process Values, Watt
source["k"]        = 0      # Interpolation degree

source["τv"]       = 0.002  # Time constant of the voltage loop, seconds
source["τf"]       = 0.002  # Time constant of the frequency loop, seconds

source["Observer"] = true   # Discrete Luenberger Observer

#= 
source["Dp"]       = 202 # frequency droop coefficient
source["Dq"]       = 6148 # voltage droop coefficient =#

#= 
source["I_kp"]     = 0.0032 # A/V
source["I_ki"]     = 0.3497 # A/Vs

source["V_kp"]     = 0.2964# A/V
source["V_ki"]     = 5.856 # A/Vs =#

push!(source_list, source)

#= 
source["Dp"]           = 202 # frequency droop coefficient
source["Dq"]           = 6148 # voltage droop coefficient
source["I_kp"]         = 0.0032 # A/V
source["I_ki"]         = 0.3497 # A/Vs
source["V_kp"]         = 0.2964# A/V
source["V_ki"]         = 5.856 # A/Vs
source["fltr"]         = "LCL"
source["control_type"] = "classic"
source["v_δ_set"]      = 0 # degrees
source["v_rip"]        = 0.01537
source["i_rip"]        = 0.15
source["vdc"]          = 800 #V
source["τv"]           = 0.002
source["τf"]           = 0.002
source["pf"]           = 0.8    # Power Factor
source["κ"]            = 3 # mean reversion parameter
source["γ"]            = 50e3   # Asymptotoic Mean
source["X₀"]           = 25 # initial values
source["L1"]           = 0.002
source["R1"]           = 0.04
source["L2"]           = 0.002
source["R2"]           = 0.05
source["R_C"]          = 0.09
source["C"]            = 0.003 =#

#= source = Dict()

source["pwr"] = 200e3
source["vdc"] = 800
source["fltr"] = "LC"
source["i_rip"] = 0.15
source["v_rip"] = 0.01537

push!(source_list, source) =#

#-------------------------------------------------------------------------------
# Loads

R_load, L_load, _, _ = Parallel_Load_Impedance(10e3, 0.6, 230)

load_list = []
load = Dict()

load["impedance"] = "RL"
load["R"] = R_load
load["L"] = L_load
load["i_limit"] = 10e12
load["v_limit"] = 10e12

push!(load_list, load)

#-------------------------------------------------------------------------------
# Network

grid = Dict()

grid["v_rms"] = 230
load["ramp_end"] = 0.04
load["process_start"] = 0.04
load["f_grid"] = 50 
load["Δfmax"] = 0.005 # The drop (increase) in frequency that causes a 100% increase (decrease) in power
load["ΔEmax"] = 0.05 # The drop (increase) in rms voltage that causes a 100% increase (decrease) in reactive power (from nominal)

#-------------------------------------------------------------------------------
# Amalgamation

parameters = Dict()

parameters["source"] = source_list
parameters["cable"] = cable_list
parameters["load"] = load_list
parameters["grid"] = grid

#_______________________________________________________________________________
# Defining the environment

env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 2)

#_______________________________________________________________________________
# Setting up data hooks

hook = DataHook(collect_vrms_ids = [1 2], 
                collect_irms_ids = [1 2], 
                collect_pq_ids   = [1 2],
                collect_freq     = [1 2])

#_______________________________________________________________________________
# Running the Time Simulation

Multi_Agent = Power_System_Dynamics(env, hook, num_episodes = num_eps, return_Agents = true)
Source = Multi_Agent.agents["classic"]["policy"].policy.Source

#_______________________________________________________________________________
# Plotting

for eps in 1:num_eps

    plot_hook_results(hook = hook, 
                      episode = eps,
                      states_to_plot  = [], 
                      actions_to_plot = [],  
                      p_to_plot       = [1 2], 
                      q_to_plot       = [], 
                      vrms_to_plot    = [1 2], 
                      irms_to_plot    = [1 2],
                      freq_to_plot    = [1 2])
end

print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")
