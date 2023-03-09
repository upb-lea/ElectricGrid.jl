using JEG
using Test


Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_end    = 0.8     # total run time, seconds
num_eps  = 1       # number of episodes to run

#-------------------------------------------------------------------------------
# Connectivity Matrix

# CM = [ 0. 0. 1.
#         0. 0. 2.
#         -1. -2. 0.]

 CM = [ 0. 1.
        -1. 0.] 

#-------------------------------------------------------------------------------
# Parameters

#= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "PQ" - grid following controllable source/load (active and reactive Power)
    3 -> "Droop" - simple grid forming with power balancing
    4 -> "Synchronverter" - enhanced droop control
=#

R_load, L_load, _, _ = Parallel_Load_Impedance(100e3, 1.0, 230)

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}( "pwr"   =>  200e3, 
                                        "mode"  =>  1,
                                        "fltr"  =>  "L",
                                        "L  "   =>  1e-3,
                                        "R1"    =>  0.1,
                                        "v_pu_set"=>    1.0,
                                        "v_δ_set"=>    0.0,
                                        ),
                        # Dict{Any, Any}("pwr" => 100e3, "mode" => 4),
                        ],
        "load"   => Any[
                        Dict{Any, Any}("impedance" => "RL", "R" => R_load, "L" => L_load),
                        ],
        # "cable"   => Any[
        #                 Dict{Any, Any}("R" => 0.208, "L" => 0.00025, "C" => 0.4e-3, "i_limit" => 10e4,),
                        # ],
        "grid" => Dict{Any, Any}("ramp_end" => 0.04)
    )
#_______________________________________________________________________________
# Defining the environment

env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 2, action_delay = 1)

#_______________________________________________________________________________
# initialising the agents 

Multi_Agent = setup_agents(env)
Source = Multi_Agent.agents["classic"]["policy"].policy.Source

#_______________________________________________________________________________
# running the time simulation 

hook = simulate(Multi_Agent, env)

#_______________________________________________________________________________
# Plotting

plot_hook_results(hook = hook, 
                    states_to_plot  = [], 
                    actions_to_plot = [],  
                    power_p         = [1 2], 
                    power_q         = [1 2], 
                    v_mag           = [1 2], 
                    i_mag           = [],
                    freq            = [1 2],
                    angles          = [1 2])

@testset "Power flow equations: power balance" begin

        # TODO: Compare P Q values 
            # Time domain : p, q
            # Freq domain : P, Q
                    # Collect steady state values from DataHook
            # p q values from DataHook
            p = 0;
            q = 0;
            for i in 1:2
                p += hook.df[!,"source$(i)_i"][end]
                q += hook.df[!,"source$(i)_p"][end]
            end

            @show p, q

            total_P_load, total_Q_load, s_load_total, total_S_source = CheckPowerBalance(parameters, num_source, num_load, CM)

            @test 
end