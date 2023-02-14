using Dare
using Test
using Logging
using CSV
using DataFrames
using Distributions

#= @testset "Classical_Controllers" begin

        #_______________________________________________________________________________
        # Network Parameters 

        #-------------------------------------------------------------------------------
        # Time simulation

        Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
        t_end    = 0.045     # total run time, seconds
        num_eps  = 3       # number of episodes to run

        #-------------------------------------------------------------------------------
        # Connectivity Matrix

        CM = [ 0. 0. 0. 1.
                0. 0. 0. 2.
                0. 0. 0. 3.
                -1. -2. -3. 0.]

        #-------------------------------------------------------------------------------
        # Cable Impedances

        cable_list = []

        cable = Dict()
        cable["R"]       = 0.208   # Ω, line resistance
        cable["L"]       = 0.00025 # H, line inductance
        cable["C"]       = 0.4e-3  # F, line capacitance
        cable["i_limit"] = 10e12   # A, line current limit

        push!(cable_list, cable, cable, cable)

        #-------------------------------------------------------------------------------
        # Sources

        source_list = []

        source = Dict()

        source["mode"]         = 2
        source["control_type"] = "classic"

        source["fltr"]         = "LCL"  # Filter type

        source["pwr"]          = 100e3  # Rated Apparent Power, VA
        source["p_set"]        = 50e3   # Real Power Set Point, Watt
        source["q_set"]        = 10e3   # Imaginary Power Set Point, VAi

        source["v_pu_set"]     = 1.00   # Voltage Set Point, p.u.
        source["v_δ_set"]      = 0      # Voltage Angle, degrees
        source["vdc"]          = 800    # V

        source["std_asy"]      = 50e3   # Asymptotic Standard Deviation
        source["σ"]            = 0.0    # Brownian motion scale i.e. ∝ diffusion, volatility parameter
        source["Δt"]           = 0.02   # Time Step, seconds
        source["X₀"]           = 0      # Initial Process Values, Watt
        source["k"]            = 0      # Interpolation degree

        source["τv"]           = 0.002  # Time constant of the voltage loop, seconds
        source["τf"]           = 0.002  # Time constant of the frequency loop, seconds

        source["Observer"]     = true   # Discrete Luenberger Observer

        source["L1"]           = 0.0003415325753131024
        source["R1"]           = 0.06830651506262048
        source["L2"]           = 4.670680436888136e-5 
        source["R2"]           = 0.009341360873776272
        source["R_C"]          = 0.17210810076926025
        source["C"]            = 0.00015412304086843381

        push!(source_list, source)

        source = Dict()

        source["mode"]         = 3
        source["control_type"] = "classic"

        source["fltr"]         = "LC"   # Filter type

        source["pwr"]          = 150e3  # Rated Apparent Power, VA
        source["p_set"]        = 50e3   # Real Power Set Point, Watt
        source["q_set"]        = 10e3   # Imaginary Power Set Point, VAi

        source["v_pu_set"]     = 1.00   # Voltage Set Point, p.u.
        source["v_δ_set"]      = 0      # Voltage Angle, degrees
        source["vdc"]          = 800    # V

        source["std_asy"]      = 50e3   # Asymptotic Standard Deviation
        source["σ"]            = 0.0    # Brownian motion scale i.e. ∝ diffusion, volatility parameter
        source["Δt"]           = 0.02   # Time Step, seconds
        source["X₀"]           = 0      # Initial Process Values, Watt
        source["k"]            = 0      # Interpolation degree

        source["τv"]           = 0.002  # Time constant of the voltage loop, seconds
        source["τf"]           = 0.002  # Time constant of the frequency loop, seconds

        source["Observer"]     = true   # Discrete Luenberger Observer

        source["L1"]           = 0.0002276883835420683
        source["R1"]           = 0.04553767670841366
        source["R_C"]          = 0.006473167716474219
        source["C"]            = 0.00023118456130265068

        push!(source_list, source)

        source = Dict()

        source["mode"]         = 4
        source["control_type"] = "classic"

        source["fltr"]         = "LCL"  # Filter type

        source["pwr"]          = 100e3  # Rated Apparent Power, VA
        source["p_set"]        = 50e3   # Real Power Set Point, Watt
        source["q_set"]        = 10e3   # Imaginary Power Set Point, VAi

        source["v_pu_set"]     = 1.00   # Voltage Set Point, p.u.
        source["v_δ_set"]      = 0      # Voltage Angle, degrees
        source["vdc"]          = 800    # V

        source["std_asy"]      = 50e3   # Asymptotic Standard Deviation
        source["σ"]            = 0.0    # Brownian motion scale i.e. ∝ diffusion, volatility parameter
        source["Δt"]           = 0.02   # Time Step, seconds
        source["X₀"]           = 0      # Initial Process Values, Watt
        source["k"]            = 0      # Interpolation degree

        source["τv"]           = 0.002  # Time constant of the voltage loop, seconds
        source["τf"]           = 0.002  # Time constant of the frequency loop, seconds

        source["Observer"]     = true   # Discrete Luenberger Observer

        source["L1"]           = 0.0003415325753131024
        source["R1"]           = 0.06830651506262048
        source["L2"]           = 4.670680436888136e-5
        source["R2"]           = 0.009341360873776272
        source["R_C"]          = 0.17210810076926025
        source["C"]            = 0.00015412304086843381

        push!(source_list, source)

        #-------------------------------------------------------------------------------
        # Loads

        R_load, L_load, _, _ = Parallel_Load_Impedance(100e3, 0.6, 230)

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
        grid["ramp_end"] = 0.04
        grid["process_start"] = 0.04
        grid["f_grid"] = 50 
        grid["Δfmax"] = 0.005 # The drop (increase) in frequency that causes a 100% increase (decrease) in power
        grid["ΔEmax"] = 0.05 # The drop (increase) in rms voltage that causes a 100% increase (decrease) in reactive power (from nominal)

        #-------------------------------------------------------------------------------
        # Amalgamation

        parameters = Dict()

        parameters["source"] = source_list
        parameters["cable"] = cable_list
        parameters["load"] = load_list
        parameters["grid"] = grid

        #_______________________________________________________________________________
        # Defining the environment

        env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 0)

        #_______________________________________________________________________________
        # Setting up data hooks

        hook = DataHook(collect_sources  = [1 2 3],
                        collect_vrms_ids = [1 2 3], 
                        collect_irms_ids = [1 2 3], 
                        collect_pq_ids   = [1 2 3],
                        collect_freq     = [1 2 3],
                        )

        #_______________________________________________________________________________
        # Running the Time Simulation

        Multi_Agent = Power_System_Dynamics(env, hook, num_episodes = num_eps, return_Agents = true)
        Source = Multi_Agent.agents["classic"]["policy"].policy.Source

        #_______________________________________________________________________________
        # Plotting

        #= for eps in 1:num_eps

                plot_hook_results(hook = hook, 
                                        episode = eps,
                                        states_to_plot  = ["source1_v_C_filt_a"], 
                                        actions_to_plot = [],  
                                        p_to_plot       = [1 2 3], 
                                        q_to_plot       = [1 2 3], 
                                        vrms_to_plot    = [1 2 3], 
                                        irms_to_plot    = [1 2 3],
                                        freq_to_plot    = [1 2 3])
        end =#

        #_______________________________________________________________________________
        # Tests 

        D = [6.710040724346115e-5 0.0003909248293754465; 202.64236728467552 6148.754619013457]
        I_kp = [0.003212489119409699; 0.0021416594129397997; 0.003212489119409699]
        I_ki = [0.34970957351408755; 0.23313971567605846; 0.34970957351408755]
        V_kp = [0.44465925382594856; 0.2964395025506326]
        V_ki = [8.78477386799335; 5.85651591199507]

        #CSV.write("Classical_Control_Unit_test.csv", hook.df)
        old_data = convert.(Float64, Matrix(CSV.read(joinpath(pwd(), "test\\Classical_Control_Unit_test.csv"), DataFrame))[1:end, 1:17])

        new_data = convert.(Float64, Matrix(hook.df)[1:end, 1:17])

        total_steps = Int(env.maxsteps)

        @test Source.D[2:end, :] ≈ D atol = 0.001
        @test Source.I_kp ≈ I_kp atol = 0.001
        @test Source.I_ki ≈ I_ki atol = 0.001
        @test Source.V_kp[2:end] ≈ V_kp atol = 0.001
        @test Source.V_ki[2:end] ≈ V_ki atol = 0.001
        @test new_data ≈ old_data atol = 0.001
        @test new_data[1:total_steps, 2:end] ≈ new_data[total_steps + 1:2*total_steps, 2:end] atol = 0.001
        @test new_data[1:total_steps, 2:end] ≈ new_data[2*total_steps + 1:end, 2:end] atol = 0.001
end =#

#= @testset "Ornstein_Uhlenbeck" begin

        #_______________________________________________________________________________
        # Network Parameters 

        #-------------------------------------------------------------------------------
        # Time simulation

        Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
        t_end    = 10     # total run time, seconds
        num_eps  = 1       # number of episodes to run

        #-------------------------------------------------------------------------------
        # Connectivity Matrix

        CM = [ 0. 0. 0. 1.
                0. 0. 0. 2.
                0. 0. 0. 3.
                -1. -2. -3. 0.]

        #-------------------------------------------------------------------------------
        # Cable Impedances

        cable_list = []

        cable = Dict()
        cable["R"]       = 0.208   # Ω, line resistance
        cable["L"]       = 0.00025 # H, line inductance
        cable["C"]       = 0.4e-3  # F, line capacitance
        cable["i_limit"] = 10e12   # A, line current limit

        push!(cable_list, cable, cable, cable)

        #-------------------------------------------------------------------------------
        # Sources

        source_list = []

        source = Dict()

        source["mode"]         = "Swing"
        source["control_type"] = "classic"

        source["fltr"]         = "LCL"  # Filter type

        source["pwr"]          = 100e3  # Rated Apparent Power, VA

        source["v_pu_set"]     = 1.00   # Voltage Set Point, p.u.
        source["v_δ_set"]      = -1.0      # Voltage Angle, degrees
        source["vdc"]          = 800    # V

        push!(source_list, source)

        source = Dict()

        source["mode"]         = "Voltage"
        source["control_type"] = "classic"

        source["fltr"]         = "LC"   # Filter type

        source["pwr"]          = 150e3  # Rated Apparent Power, VA

        source["v_pu_set"]     = 1.00   # Voltage Set Point, p.u.
        source["v_δ_set"]      = 3.5      # Voltage Angle, degrees
        source["vdc"]          = 800    # V

        push!(source_list, source)

        source = Dict()

        source["mode"]         = "PQ"
        source["control_type"] = "classic"

        source["fltr"]         = "L"    # Filter type

        source["pwr"]          = 100e3  # Rated Apparent Power, VA
        source["p_set"]        = 50e3   # Real Power Set Point, Watt
        source["q_set"]        = -10e3   # Imaginary Power Set Point, VAi

        source["v_pu_set"]     = 1.00   # Voltage Set Point, p.u.
        source["vdc"]          = 800    # V

        source["std_asy"]      = 5e3   # Asymptotic Standard Deviation
        source["σ"]            = 25e3    # Brownian motion scale i.e. ∝ diffusion, volatility parameter
        source["Δt"]           = 0.02   # Time Step, seconds
        source["X₀"]           = 30e3   # Initial Process Values, Watt
        source["k"]            = 1      # Interpolation degree
        source["γ"]            = 50e3   # Asymptotoic Mean

        push!(source_list, source)

        #-------------------------------------------------------------------------------
        # Loads

        R_load, L_load, _, _ = Parallel_Load_Impedance(100e3, 0.95, 230)

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
        grid["ramp_end"] = 0.03
        grid["process_start"] = 0.04
        grid["f_grid"] = 50 

        #-------------------------------------------------------------------------------
        # Amalgamation

        parameters = Dict()

        parameters["source"] = source_list
        parameters["cable"] = cable_list
        parameters["load"] = load_list
        parameters["grid"] = grid

        #_______________________________________________________________________________
        # Defining the environment

        env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 0)

        #_______________________________________________________________________________
        # Setting up data hooks

        hook = DataHook(collect_sources  = [1 2 3],
                        collect_vrms_ids = [1 2 3], 
                        collect_irms_ids = [1 2 3], 
                        collect_pq_ids   = [1 2 3],
                        collect_freq     = [1 2 3],
                        collect_θ        = [1 2 3],
                        )

        #_______________________________________________________________________________
        # Running the Time Simulation

        Multi_Agent = Power_System_Dynamics(env, hook, num_episodes = num_eps, return_Agents = true)
        Source = Multi_Agent.agents["classic"]["policy"].policy.Source

        #_______________________________________________________________________________
        # Plotting

        #= for eps in 1:num_eps

                plot_hook_results(hook = hook, 
                                        episode = eps,
                                        states_to_plot  = [], 
                                        actions_to_plot = [],  
                                        p_to_plot       = [3], 
                                        q_to_plot       = [3], 
                                        vrms_to_plot    = [], 
                                        irms_to_plot    = [],
                                        freq_to_plot    = [3],
                                        θ_to_plot       = [1 2 3])
        end =#

        #_______________________________________________________________________________
        # Tests 

        step = Int(parameters["source"][3]["Δt"]/Timestep)
        start = Int(parameters["grid"]["process_start"]/Timestep) + 1 + step

        new_data = convert.(Float64, Matrix(hook.df)[start:step:end, 7])
        new_angles = convert.(Float64, Matrix(hook.df)[:, 18:19])

        angles_eval = ones(size(new_angles, 1), 2)
        angles_eval[:,1] = parameters["source"][1]["v_δ_set"].*angles_eval[:,1]
        angles_eval[:,2] = parameters["source"][2]["v_δ_set"].*angles_eval[:,2]

        stats = fit(Normal{Float32}, new_data)

        @test 1 ≈ stats.μ/parameters["source"][3]["γ"] atol = 0.01
        @test 1 ≈ stats.σ/parameters["source"][3]["std_asy"] atol = 0.1
        @test new_angles ≈ angles_eval atol = 0.001

        return nothing
end =#

@testset "Saturation" begin

        return nothing
end
#_______________________________________________________________________________
# Network Configuration 

#-------------------------------------------------------------------------------
# Time simulation

Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_end    = 0.1     # total run time, seconds

#-------------------------------------------------------------------------------
# Connectivity Matrix

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

R_load, L_load, _, _ = Parallel_Load_Impedance(100e3, 0.99, 230)

length = 1
parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => 200e3, 
                                        "mode" => "Semi-Synchronverter", 
                                        "v_pu_set" => 1.05),

                        Dict{Any, Any}("pwr" => 200e3, 
                                        "mode" => "PQ", 
                                        "p_set" => -54.68e3, # making this slightly less/more, means that the voltage control loop recovers
                                        "q_set" => 100e3),
                        ],
        "cable"   => Any[
                        Dict{Any, Any}("R" => length*0.208, 
                                        "L" => length*0.00025, 
                                        "C" => length*0.4e-3, 
                                        "i_limit" => 10e4,),
                        ],
        "grid" => Dict{Any, Any}("ramp_end" => 0.04, "process_start"=> 0.06)
    )
#_______________________________________________________________________________
# Defining the environment

env = SimEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 2)

#_______________________________________________________________________________
# Setting up data hooks

hook = DataHook(collect_vrms_ids = [1 2], 
                collect_irms_ids = [1 2], 
                collect_pq_ids   = [1 2],
                collect_freq     = [1 2],
                collect_sources  = [1 2],
                collect_θ        = [1 2],
                collect_debug    = [1 2 5 6])

#_______________________________________________________________________________
# Running the Time Simulation

Multi_Agent = Power_System_Dynamics(env, hook; return_Agents = true)
Source = Multi_Agent.agents["classic"]["policy"].policy.Source

#_______________________________________________________________________________
# Plotting

plot_hook_results(hook = hook, 
                    states_to_plot  = [], 
                    actions_to_plot = [],  
                    p_to_plot       = [1 2], 
                    q_to_plot       = [1 2], 
                    vrms_to_plot    = [1 2], 
                    irms_to_plot    = [1 2],
                    freq_to_plot    = [],
                    θ_to_plot       = [])