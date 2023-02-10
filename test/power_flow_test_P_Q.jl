using Dare
using Test

#= 
    Rationale for this unit test:
    - The goal is to compare the P Q values from the time domain simulation with the P Q values from the frequency domain simulation
    - The frequency domain simulation is run with the Power_System_Dynamics() function
    - The time domain simulation is run with the SimEnv() function
=#



# with basic commands without SmallWorld
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
            # invoke layout_cabels() by not specifying the cables in parameters
            parameters = Dict{Any, Any}(
                    "source" => Any[
                                    Dict{Any, Any}("pwr" => 200e3, "mode" => 1),
                                    Dict{Any, Any}("pwr" => 200e3, "mode" => 1)
                                    ],
                    "load"   => Any[
                                    Dict{Any, Any}("impedance" => "RL", "R" => 2.64, "L" => 0.006) 
                                    ],
                #     "cable"   => Any[
                #                     Dict{Any, Any}("R" => 1e-3, "L" => 1e-4, "C" => 1e-4, "i_limit" => 10e4,) #covered in hotfix
                #                     ]
                )
            #_______________________________________________________________________________
        #     TODO cable layout

            
            # Defining the environment
            env = SimEnv(   ts = Timestep, 
                            CM = CM, 
                            parameters = parameters, 
                            t_end = t_end)

            #_______________________________________________________________________________
            # Setting up data hooks

            hook = DataHook(collect_vrms_ids = [1 2], 
                            collect_irms_ids = [1 2], 
                            collect_pq_ids   = [1 2], #collecting p and q for sources 1, 2
                            collect_freq     = [1 2],
                            collect_sources  = [1 2])

            #_______________________________________________________________________________
            # Running the Time Simulation

            Power_System_Dynamics(env, hook)

            #_______________________________________________________________________________
            # Plotting

            plot_hook_results(hook = hook, 
                                states_to_plot  = [], 
                                actions_to_plot = [],  
                                p_to_plot       = [1 2], 
                                q_to_plot       = [],  #no plotting q
                                vrms_to_plot    = [1 2], 
                                irms_to_plot    = [1 2],
                                freq_to_plot    = [1 2])

            print("\n...........o0o----ooo0§0ooo~~~   END   ~~~ooo0§0ooo----o0o...........\n")

        # TODO: Compare P Q values 
            # Time domain : p, q from DataHook
            # Freq domain : P, Q from optimizer

        # Collect steady state values from DataHook
            # p q values from DataHook
            p = 0;
            q = 0;
            for i in 1:2
                p += hook.df[!,"source$(i)_i"][end]
                q += hook.pq[!,"source$(i)_p"][end]
            end

            @show p, q

            total_P_load, total_Q_load, s_load_total, total_S_source = CheckPowerBalance(parameters, num_source, num_load, CM)


end