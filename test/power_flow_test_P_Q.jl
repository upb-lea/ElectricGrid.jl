using JEG
using Test

import JEG


Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_end = 0.8     # total run time, seconds
num_eps = 1       # number of episodes to run

#-------------------------------------------------------------------------------
# Connectivity Matrix

# CM = [ 0. 0. 1.
#         0. 0. 2.
#         -1. -2. 0.]

CM = [0.0 1.0
    -1.0 0.0]

#-------------------------------------------------------------------------------
# Parameters

#= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "PQ" - grid following controllable source/load (active and reactive Power)
    3 -> "Droop" - simple grid forming with power balancing
    4 -> "Synchronverter" - enhanced droop control
=#


R_load, L_load, _, _ = ParallelLoadImpedance(100e3, 0.95, 230)

parameters = Dict{Any,Any}(
    "source" => Any[
        Dict{Any,Any}(
            "pwr" => 200e3,
            "mode" => 1,
            "fltr" => "L",
            "L1" => 1e-6,
            "R1" => 1e-6,
            "v_pu_set" => 1.0,
            "v_δ_set" => 0.0,
        ),
    ],

    "load" => Any[
        Dict{Any,Any}(
            "impedance" => "RL", 
            "R" => R_load, 
            "L" => L_load
        ),
    ],

    "grid" => Dict{Any,Any}(
                "ramp_end" => 0.04,
                "phase" => 3)
)
#_______________________________________________________________________________
# Defining the environment

env = ElectricGridEnv(ts=Timestep, CM=CM, parameters=parameters, t_end=t_end, verbosity=2, action_delay=1)

#_______________________________________________________________________________
# initialising the agents 

Multi_Agent = SetupAgents(env)
Source = Multi_Agent.agents["classic"]["policy"].policy.Source

#_______________________________________________________________________________
# running the time simulation 

hook = Simulate(Multi_Agent, env)

#_______________________________________________________________________________
# Plotting

RenderHookResults(hook=hook,
    states_to_plot=[],
    actions_to_plot=[],
    power_p_inv=[1],
    power_q_inv=[1],
    v_mag_inv=[1],
    i_mag_inv=[],
    freq=[1],
    angles=[1])

@testset "Power flow equations: power balance" begin

    # TODO: Compare P Q values 
    # Time domain : p, q
    # Freq domain : P, Q
    # Collect steady state values from DataHook
    # p q values from DataHook
    
    num_sources = env.nc.num_sources
    num_loads = env.nc.num_loads
    parameters = env.nc.parameters

    p_test = 0
    q_test = 0

    for i in 1:env.nc.num_sources
        p_test += hook.df[!, "source$(i)_p_inv"][end]
        q_test += hook.df[!, "source$(i)_q_inv"][end]
    end

    @show p_test, q_test

    # total_P_load, total_Q_load, s_load_total, total_S_source = CheckPowerBalance(parameters, num_sources, num_loads, CM)

    @show 3*[JEG.P_source, JEG.Q_source]
    
    @test (p_test - 3 * JEG.P_source[1]) < 1e1 
    # @test q_test ≈ 3 * JEG.Q_source 
end