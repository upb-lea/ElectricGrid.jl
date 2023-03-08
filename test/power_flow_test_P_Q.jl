using JEG
using Test

@testset "Layout cables with Power flow equations" begin
 
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

            parameters = Dict{Any, Any}(
                    "source" => Any[
                                    Dict{Any, Any}("pwr" => 200e3),
                                    Dict{Any, Any}("pwr" => 200e3)
                                    ],
                    "load"   => Any[
                                    Dict{Any, Any}("impedance" => "RL", "R" => 2.64, "L" => 0.006) 
                                    ]
                )
            #_______________________________________________________________________________
        #     TODO cable layout


            # Defining the environment
            env = ElectricGridEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end)

            #_______________________________________________________________________________
            # Setting up data hooks

            hook = data_hook(collect_vrms_ids = [1 2], 
                            collect_irms_ids = [1 2], 
                            collect_pq_ids   = [1 2], #collecting p and q for sources 1, 2
                            collect_freq     = [1 2],
                            collect_sources  = [1 2])

            #_______________________________________________________________________________
            # Running the Time Simulation

            Power_System_Dynamics(env, hook)

            #_______________________________________________________________________________
            # Plotting

            RenderHookResults(hook = hook, 
                                states_to_plot  = [], 
                                actions_to_plot = [],  
                                p_to_plot       = [1 2], 
                                q_to_plot       = [],  #no plotting q
                                vrms_to_plot    = [1 2], 
                                irms_to_plot    = [1 2],
                                freq_to_plot    = [1 2])

            print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")

        # TODO: Compare P Q values 
            # Time domain : p, q
            # Freq domain : P, Q
end