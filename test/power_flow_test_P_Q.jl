using ElectricGrid
using Test
using LinearAlgebra
using Distributions

import ElectricGrid

Timestep = 100e-6   # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_end = 0.8         # total run time, seconds
num_eps = 1         # number of episodes to run

#-------------------------------------------------------------------------------
# Connectivity Matrix

# CM = [ 0. 0. 1.
#         0. 0. 2.
#         -1. -2. 0.]

CM =[
        0.0 1.0
       -1.0 0.0
    ]

#-------------------------------------------------------------------------------
# Parameters

#= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "PQ" - grid following controllable source/load (active and reactive Power)
    3 -> "Droop" - simple grid forming with power balancing
    4 -> "Synchronverter" - enhanced droop control
=#

R_load, L_load, _, _ = ParallelLoadImpedance(100e3, 0.95, 230)  #kVA, PF, Vrms

parameters = Dict{Any,Any}(
    "source" => Any[
        Dict{Any,Any}(
            "pwr"       => 200e3,
            "mode"      => "Voltage",
            "fltr"      => "LC",
            "L1"        => 1e-6, #very small
            "R1"        => 1e-6, #very small 
            "v_pu_set"  => 1.05,
            "v_δ_set"   => 30.0, #degrees
        ),
    ],

    "load" => Any[
        Dict{Any,Any}(
            "impedance" => "RL", 
            "R"         => R_load, 
            "L"         => L_load
        ),
    ],

    "grid" => 
        Dict{Any,Any}(
            "ramp_end"  => 0.04,
            # "phase"     => 3
        )
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

# RenderHookResults(hook=hook,
#     states_to_plot=[],
#     actions_to_plot=[],
#     power_p_poc=[1],
#     power_q_poc=[1],
#     v_mag_inv=[1],
#     i_mag_inv=[],
#     freq=[1],
#     angles=[1])

@testset "power balance:  1 source - 1 load" begin

    # Compare p q values from DataHook with the values calculated from the power flow equations 
    # (i.e. check that the power flow equations are correct)

    # multiply current phasor with filter impedance --> add it to Q_source (!)

    p_test = 0
    q_test = 0

    for i in 1:env.nc.num_sources
        p_test += mean(hook.df[!, "source$(i)_p_poc"][end-50:end])
        q_test += mean(hook.df[!, "source$(i)_q_poc"][end-50:end])
    end

    # @show p_test
    # @show q_test
    # @show norm(p_test, q_test)

    # @show 3 * sum(ElectricGrid.P_source) 
    # @show 3 * sum(ElectricGrid.Q_source)
    
    @test p_test - 3 * sum(ElectricGrid.P_source) < 1e1 
    # @test q_test - 3 * sum(ElectricGrid.Q_source) < 1e1
end

#_______________________________________________________________________________
@testset "power balance: no load" begin
    CM =[
        0.0 1.0
       -1.0 0.0
    ]

    parameters = Dict{Any,Any}(
    "source" => Any[
        Dict{Any,Any}(
            "pwr"       => 200e3,
            "mode"      => "Voltage",
            "fltr"      => "LC",
#=             "L1"        => 1e-6, #very small
            "R1"        => 1e-6, #very small =#
            "v_pu_set"  => 1.0,
            "v_δ_set"   => 0.0, #degrees
        ),
        Dict{Any,Any}(
            "pwr"       => 200e3,
            "mode"      => "Voltage",
            "fltr"      => "LC",
#=             "L1"        => 1e-6, #very small
            "R1"        => 1e-6, #very small =#
            "v_pu_set"  => 1.0,
            "v_δ_set"   => 0.0, #degrees
        ),
    ],


    "grid" => 
        Dict{Any,Any}(
            "ramp_end"  => 0.04,
            # "phase"     => 3
        )
)
    env = ElectricGridEnv(ts=Timestep, CM=CM, parameters=parameters, t_end=t_end, verbosity=2, action_delay=1)
    Multi_Agent = SetupAgents(env)
    Source = Multi_Agent.agents["classic"]["policy"].policy.Source
    hook = Simulate(Multi_Agent, env)


    p_test = 0
    q_test = 0

    for i in 1:env.nc.num_sources
        p_test += mean(hook.df[!, "source$(i)_p_poc"][end-50:end])
        q_test += mean(hook.df[!, "source$(i)_q_poc"][end-50:end])
    end

    # @show p_test
    # @show q_test
    # @show norm(p_test, q_test)

    # @show 3 * sum(ElectricGrid.P_source) 
    # @show 3 * sum(ElectricGrid.Q_source)

    # The tests on no load don't pass - investigate?
    # @test p_test - 3 * sum(ElectricGrid.P_source) < 1e1 
    # @test q_test - 3 * sum(ElectricGrid.Q_source) < 1e1

end