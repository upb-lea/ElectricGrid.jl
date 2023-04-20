using ElectricGrid
using Test
using Logging
using CSV
using DataFrames
using Distributions

@testset "ClassicalControllers_Dynamics" begin

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

    #= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "PQ" - grid following controllable source/load (active and reactive Power)
    3 -> "Droop" - simple grid forming with power balancing
    4 -> "Synchronverter" - enhanced droop control
    =#

    source_list = []

    source = Dict()

    source["mode"]         = "PQ"
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

    source["mode"]         = "Droop"
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

    source["mode"]         = "Synchronverter"
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

    R_load, L_load, _, _ = ParallelLoadImpedance(100e3, 0.6, 230)

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
    grid["Δfmax"] = 0.5 # The drop (increase) in frequency that causes a 100% increase (decrease) in power
    grid["ΔEmax"] = 5 # The drop (increase) in rms voltage that causes a 100% increase (decrease) in reactive power (from nominal)

    #-------------------------------------------------------------------------------
    # Amalgamation

    parameters = Dict()

    parameters["source"] = source_list
    parameters["cable"] = cable_list
    parameters["load"] = load_list
    parameters["grid"] = grid

    #_______________________________________________________________________________
    # Defining the environment

    env = ElectricGridEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 0)

    #_______________________________________________________________________________
    # Setting up data hooks

    hook = DataHook(collect_sources     = [1 2 3],
                    v_mag_inv            = [1 2 3],
                    v_mag_poc            = [1 2 3],
                    i_mag_inv            = [1 2 3],
                    i_mag_poc            = [1 2 3],
                    power_pq_inv         = [1 2 3],
                    power_pq_poc         = [1 2 3],
                    freq                 = [1 2 3],
                    )

    #_______________________________________________________________________________
    # initialising the agents

    Multi_Agent = SetupAgents(env)
    Source = Multi_Agent.agents["classic"]["policy"].policy.Source

    #_______________________________________________________________________________
    # running the time simulation

    hook = Simulate(Multi_Agent, env, num_episodes = num_eps, hook = hook)

    #_______________________________________________________________________________
    # Plotting

    #= for eps in 1:num_eps

            RenderHookResults(hook = hook,
                                    episode = eps,
                                    states_to_plot  = ["source1_v_C_filt_a"],
                                    actions_to_plot = [],
                                    power_p_inv     = [1 2 3],
                                    power_q_inv     = [1 2 3],
                                    power_p_poc     = [1 2 3],
                                    power_q_poc     = [1 2 3],
                                    v_mag_inv       = [1 2 3],
                                    v_mag_poc       = [1 2 3],
                                    i_mag_inv       = [1 2 3],
                                    i_mag_poc       = [1 2 3],
                                    freq            = [1 2 3])
    end =#

    #_______________________________________________________________________________
    # Tests

    D = [303.9635509270133 9223.131928520184; 202.64236728467552 6148.754619013457]
    I_kp = [0.003212489119409699; 0.0021416594129397997; 0.003212489119409699]
    I_ki = [0.34970957351408755; 0.23313971567605846; 0.34970957351408755]
    V_kp = [0.44465925382594856; 0.2964395025506326]
    V_ki = [8.78477386799335; 5.85651591199507]

    #CSV.write("Classical_Control_Unit_test.csv", hook.df)
    old_data = convert.(Float64, Matrix(CSV.read("./test/Classical_Control_Unit_test.csv", DataFrame))[1:end, 1:17])

    new_data = convert.(Float64, Matrix(hook.df)[1:end, 1:17])

    total_steps = Int(env.maxsteps)

    @test Source.D[2:end, :] ≈ D atol = 0.001
    @test Source.I_kp ≈ I_kp atol = 0.001
    @test Source.I_ki ≈ I_ki atol = 0.001
    @test Source.V_kp[2:end] ≈ V_kp atol = 0.001
    @test Source.V_ki[2:end] ≈ V_ki atol = 0.001
    @test ((1 .+ new_data)./(1 .+ old_data)) ≈ ones(size(new_data, 1),size(new_data, 2)) atol = 0.01

    @test new_data[1:total_steps, 2:end] ≈ new_data[total_steps + 1:2*total_steps, 2:end] atol = 0.001
    @test new_data[1:total_steps, 2:end] ≈ new_data[2*total_steps + 1:end, 2:end] atol = 0.001
end

@testset "Saturation" begin

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

    R_load, L_load, _, _ = ParallelLoadImpedance(100e3, 0.99, 230)

    length = 1
    parameters = Dict{Any, Any}(
            "source" => Any[
                            Dict{Any, Any}("pwr" => 200e3,
                                            "mode" => "Semi-Synchronverter",
                                            "v_pu_set" => 1.0),

                            Dict{Any, Any}("pwr" => 200e3,
                                            "mode" => "PQ",
                                            "fltr" => "LCL",
                                            "p_set" => -40.9e3, # making this slightly less/more, means that the voltage control loop recovers
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

    env = ElectricGridEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 0)

    #_______________________________________________________________________________
    # Setting up data hooks

    hook = DataHook(collect_sources  = [1 2],
                    v_mag_inv            = [1 2],
                    v_mag_poc            = [1 2],
                    i_mag_inv            = [1 2],
                    i_mag_poc            = [1 2],
                    power_pq_inv         = [1 2],
                    power_pq_poc         = [1 2],
                    freq             = [1 2],
                    angles           = [1 2],
                    i_sat            = [1 2],
                    v_sat            = [1],
                    i_err_t          = [1 2],
                    v_err_t          = [1])

    #_______________________________________________________________________________
    # initialising the agents

    Multi_Agent = SetupAgents(env)
    Source = Multi_Agent.agents["classic"]["policy"].policy.Source

    #_______________________________________________________________________________
    # running the time simulation

    hook = Simulate(Multi_Agent, env, hook = hook)

    #_______________________________________________________________________________
    # Plotting

    #= RenderHookResults(hook = hook,
                    states_to_plot  = [],
                    actions_to_plot = [],
                    power_p_inv     = [],
                    power_q_inv     = [],
                    power_p_poc     = [],
                    power_q_poc     = [],
                    v_mag_inv           = [1 2],
                    v_mag_poc           = [1 2],
                    i_mag_inv           = [1 2],
                    i_mag_poc           = [1 2],
                    i_sat           = [1 2],
                    v_sat           = [1],
                    i_err_t         = [1 2],
                    v_err_t         = [1],
                    freq            = [],
                    angles          = []) =#

    #CSV.write("Control_Saturation_Unit_test.csv", hook.df[end-100:end,:])
    old_data = convert.(Float64, Matrix(CSV.read("./test/Control_Saturation_Unit_test.csv", DataFrame))[1:end, 1:17])

    new_data = convert.(Float64, Matrix(hook.df)[end-100:end, 1:17])

    @test ((1 .+ new_data)./(1 .+ old_data)) ≈ ones(size(new_data, 1),size(new_data, 2)) atol = 0.01

    return nothing
end

@testset "OrnsteinUhlenbeck_Filters_Angles" begin

    #_______________________________________________________________________________
    # Network Parameters

    #-------------------------------------------------------------------------------
    # Time simulation

    Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
    t_end    = 20     # total run time, seconds
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
    cable["R"]       = 0.01   # Ω, line resistance
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
    source["Δt"]           = 0.01   # Time Step, seconds
    source["X₀"]           = 50e3   # Initial Process Values, Watt
    source["k"]            = 1      # Interpolation degree
    source["γ"]            = 50e3   # Asymptotoic Mean

    push!(source_list, source)

    #-------------------------------------------------------------------------------
    # Loads

    R_load, L_load, _, _ = ParallelLoadImpedance(100e3, 0.95, 230)

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

    env = ElectricGridEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 0)

    #_______________________________________________________________________________
    # Setting up data hooks

    hook = DataHook(collect_sources  = [1 2 3],
                        v_mag_inv            = [1 2 3],
                        i_mag_inv            = [1 2 3],
                        power_pq_inv     = [1 2 3],
                        power_pq_poc     = [1 2 3],
                        freq             = [1 2 3],
                        angles           = [1 2 3],
                        )

    #_______________________________________________________________________________
    # initialising the agents

    Multi_Agent = SetupAgents(env)
    Source = Multi_Agent.agents["classic"]["policy"].policy.Source

    #_______________________________________________________________________________
    # running the time simulation

    hook = Simulate(Multi_Agent, env, num_episodes = num_eps, hook = hook)

    #_______________________________________________________________________________
    # Plotting

    #= for eps in 1:num_eps

            RenderHookResults(hook = hook,
                                    episode = eps,
                                    states_to_plot  = [],
                                    actions_to_plot = [],
                                    power_p_inv         = [3],
                                    power_q_inv         = [3],
                                    v_mag_inv           = [],
                                    i_mag_inv           = [],
                                    freq            = [3],
                                    angles          = [1 2 3])
    end =#

    #_______________________________________________________________________________
    # Tests

    s1_L1 = 0.0006830651506262048
    s1_R1 = 0.13661303012524095
    s1_L2 = 9.341360873776272e-5
    s1_R2 = 0.018682721747552544
    s1_C = 7.706152043421691e-5
    s1_R_C = 0.3442162015385205

    s2_L1 = 0.0004553767670841366
    s2_R1 = 0.09107535341682732
    s2_C = 0.00011559228065132534
    s2_R_C = 0.22947746769234703

    s3_L1 = 0.0006830651506262048
    s3_R1 = 0.13661303012524095

    step = Int(parameters["source"][3]["Δt"]/Timestep)
    start = Int(parameters["grid"]["process_start"]/Timestep) + 1 + step

    new_data = convert.(Float64, (hook.df[!,"source3_p_inv"][start:step:end]))
    new_angles = convert.(Float64, Matrix(hook.df[!,["source1_θ", "source2_θ"]]))

    angles_eval = ones(size(new_angles, 1), 2)
    angles_eval[:,1] = parameters["source"][1]["v_δ_set"].*angles_eval[:,1]
    angles_eval[:,2] = parameters["source"][2]["v_δ_set"].*angles_eval[:,2]

    stats = fit(Normal{Float32}, new_data)

    @test 1 ≈ stats.μ/parameters["source"][3]["γ"] atol = 0.05
    @test 1 ≈ stats.σ/parameters["source"][3]["std_asy"] atol = 0.1
    @test new_angles ≈ angles_eval atol = 0.001

    @test s1_L1 ≈ env.nc.parameters["source"][1]["L1"] atol = 0.00001
    @test s1_R1 ≈ env.nc.parameters["source"][1]["R1"] atol = 0.00001
    @test s1_L2 ≈ env.nc.parameters["source"][1]["L2"] atol = 0.00001
    @test s1_R2 ≈ env.nc.parameters["source"][1]["R2"] atol = 0.00001
    @test s1_C ≈ env.nc.parameters["source"][1]["C"] atol = 0.00001
    @test s1_R_C ≈ env.nc.parameters["source"][1]["R_C"] atol = 0.00001

    @test s2_L1 ≈ env.nc.parameters["source"][2]["L1"] atol = 0.00001
    @test s2_R1 ≈ env.nc.parameters["source"][2]["R1"] atol = 0.00001
    @test s2_C ≈ env.nc.parameters["source"][2]["C"] atol = 0.00001
    @test s2_R_C ≈ env.nc.parameters["source"][2]["R_C"] atol = 0.00001

    @test s3_L1 ≈ env.nc.parameters["source"][3]["L1"] atol = 0.00001
    @test s3_R1 ≈ env.nc.parameters["source"][3]["R1"] atol = 0.00001

    #= println("stats.μ = ", stats.μ)
    println("γ = ", parameters["source"][3]["γ"])
    println("stats.σ = ", stats.σ)
    println("std_asy = ", parameters["source"][3]["std_asy"])
    println()
    println(env.nc.parameters["source"][1]["L1"])
    println(env.nc.parameters["source"][1]["R1"])
    println(env.nc.parameters["source"][1]["L2"])
    println(env.nc.parameters["source"][1]["R2"])
    println(env.nc.parameters["source"][1]["C"])
    println(env.nc.parameters["source"][1]["R_C"])
    println()
    println(env.nc.parameters["source"][2]["L1"])
    println(env.nc.parameters["source"][2]["R1"])
    println(env.nc.parameters["source"][2]["C"])
    println(env.nc.parameters["source"][2]["R_C"])
    println()
    println(env.nc.parameters["source"][3]["L1"])
    println(env.nc.parameters["source"][3]["R1"]) =#

    return nothing
end

@testset "Droop_frequency_v_mag" begin

    #_______________________________________________________________________________
    # Network Configuration

    #-------------------------------------------------------------------------------
    # Time simulation

    Timestep = 100e-6  # time step, seconds ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
    t_end    = 2.0     # total run time, seconds
    num_eps  = 1       # number of episodes to run

    #-------------------------------------------------------------------------------
    # Connectivity Matrix

    CM = [ 0. 0. 0. 1.
            0. 0. 0. 2.
            0. 0. 0. 3.
            -1. -2. -3. 0.]

    #= CM = [ 0. 0. 1.
            0. 0. 2.
            -1. -2. 0.] =#

    #= CM = [ 0. 1.
            -1. 0.] =#

    #-------------------------------------------------------------------------------
    # Parameters

    #= Modes:
    1 -> "Swing" - voltage source without dynamics (i.e. an Infinite Bus)
    2 -> "PQ" - grid following controllable source/load (active and reactive Power)
    3 -> "Droop" - simple grid forming with power balancing
    4 -> "VSG" - enhanced droop control
    =#

    R_load, L_load, _, _ = ParallelLoadImpedance(100e3, 0.99, 100)

    parameters = Dict{Any, Any}(
            "source" => Any[
                            Dict{Any, Any}("pwr" => 200e3,
                                            "mode" => "VSG",
                                            "v_pu_set" => 1.0,
                                            "vdc" => 800,
                                            "τv"       => 0.02,  # Time constant of the voltage loop, seconds
                                            "τf"       => 0.02,  # Time constant of the frequency loop, seconds
                                            "i_rip" => 0.02),
                            Dict{Any, Any}("pwr" => 150e3,
                                            "mode" => "Droop",
                                            "v_pu_set" => 1.0,
                                            "vdc" => 800,
                                            "τv"       => 0.02,  # Time constant of the voltage loop, seconds
                                            "τf"       => 0.002,  # Time constant of the frequency loop, seconds
                                            "i_rip" => 0.02),
                            Dict{Any, Any}("pwr" => 100e3,
                                            "mode" => "PQ",
                                            "p_set" => 40e3,
                                            "q_set" => 20e3,
                                            "i_rip" => 0.1),
                            ],
            "load"   => Any[
                            Dict{Any, Any}("impedance" => "RL",
                                            "R" => R_load,
                                            "L" => L_load),
                            ],
            "cable"   => Any[
                            Dict{Any, Any}("R" => 0.1,
                                            "L" => 0.25e-4,
                                            "C" => 0.1e-5,
                                            "i_limit" => 10e4),
                            Dict{Any, Any}("R" => 0.1,
                                            "L" => 0.25e-4,
                                            "C" => 0.1e-5,
                                            "i_limit" => 10e4,),
                            Dict{Any, Any}("R" => 0.1,
                                            "L" => 0.25e-4,
                                            "C" => 0.1e-5,
                                            "i_limit" => 10e4,),
                            ],
            "grid" => Dict{Any, Any}("ramp_end" => 0.04,
                                    "f_grid" => 60,
                                    "process_start" => 1.0,
                                    "v_rms" => 100,
                                    "Δfmax" => 0.5, # The drop (increase) in frequency that causes a 100% increase (decrease) in active power (from nominal)
                                    "ΔEmax" => 10, # The drop (increase) in rms voltage that causes a 100% increase (decrease) in reactive power (from nominal)
                                    )
    )
    #_______________________________________________________________________________
    # Defining the environment

    env = ElectricGridEnv(ts = Timestep, CM = CM, parameters = parameters, t_end = t_end, verbosity = 0, action_delay = 1)

    #_______________________________________________________________________________
    # initialising the agents

    Multi_Agent = SetupAgents(env)
    Source = Multi_Agent.agents["classic"]["policy"].policy.Source

    #_______________________________________________________________________________
    # running the time simulation

    hook = Simulate(Multi_Agent, env, num_episodes = num_eps)

    total_steps = Int(env.maxsteps)
    #_______________________________________________________________________________
    # Plotting

    #= RenderHookResults(hook = hook,
                    states_to_plot  = [],
                    actions_to_plot = [],
                    power_p_inv     = [1 2 3],
                    power_q_inv     = [1 2 3],
                    power_p_poc     = [1 2 3],
                    power_q_poc     = [1 2 3],
                    v_mag_inv       = [1 2 3],
                    v_mag_poc       = [1 2 3],
                    i_mag_inv       = [],
                    i_mag_poc       = [],
                    freq            = [1 2 3],
                    angles          = [1 2 3],
                    i_sat           = [],
                    v_sat           = [],
                    i_err_t         = [],
                    v_err_t         = []) =#

    #_______________________________________________________________________________
    # tests
    total_steps = Int(env.maxsteps)
    step_start = convert(Int, round(env.nc.parameters["grid"]["process_start"]/Timestep)) - 5

    # VSG
    freq_start_1 = hook.df[!,"source1_freq"][step_start]
    freq_end_1 = hook.df[!,"source1_freq"][total_steps]
    freq_Δ_1 = 100*(freq_start_1 - freq_end_1)/(env.nc.parameters["grid"]["Δfmax"]*env.nc.parameters["grid"]["f_grid"])

    vrms_start_1 = hook.df[!,"source1_v_mag_poc"][step_start]
    vrms_end_1 = hook.df[!,"source1_v_mag_poc"][total_steps]
    vrms_Δ_1 = 100*(vrms_start_1 - vrms_end_1)/(env.nc.parameters["grid"]["ΔEmax"]*env.nc.parameters["grid"]["v_rms"]*env.nc.parameters["source"][1]["v_pu_set"])

    p_start_1 = hook.df[!,"source1_p_poc"][step_start]
    p_end_1 = hook.df[!,"source1_p_poc"][total_steps]
    p_Δ_1 = - p_start_1 + p_end_1

    q_start_1 = hook.df[!,"source1_q_poc"][step_start]
    q_end_1 = hook.df[!,"source1_q_poc"][total_steps]
    q_Δ_1 = - q_start_1 + q_end_1

    dP_1 = p_Δ_1/Source.S[1]

    dQ_1 = q_Δ_1/Source.S[1]

    # Droop
    freq_start_2 = hook.df[!,"source2_freq"][step_start]
    freq_end_2 = hook.df[!,"source2_freq"][total_steps]
    freq_Δ_2 = 100*(freq_start_2 - freq_end_2)/(env.nc.parameters["grid"]["Δfmax"]*env.nc.parameters["grid"]["f_grid"])

    vrms_start_2 = hook.df[!,"source2_v_mag_poc"][step_start]
    vrms_end_2 = hook.df[!,"source2_v_mag_poc"][total_steps]
    vrms_Δ_2 = 100*(vrms_start_2 - vrms_end_2)/(env.nc.parameters["grid"]["ΔEmax"]*env.nc.parameters["grid"]["v_rms"]*env.nc.parameters["source"][2]["v_pu_set"])

    p_start_2 = hook.df[!,"source2_p_poc"][step_start]
    p_end_2 = hook.df[!,"source2_p_poc"][total_steps]
    p_Δ_2 = - p_start_2 + p_end_2

    q_start_2 = hook.df[!,"source2_q_poc"][step_start]
    q_end_2 = hook.df[!,"source2_q_poc"][total_steps]
    q_Δ_2 = - q_start_2 + q_end_2

    dP_2 = p_Δ_2/Source.S[2]

    dQ_2 = q_Δ_2/Source.S[2]

    #= @show freq_Δ_1
    @show dP_1

    @show vrms_Δ_1
    @show dQ_1

    @show freq_Δ_2
    @show dP_2

    @show vrms_Δ_2
    @show dQ_2

    @show hook.df[!,"source3_p_inv"][total_steps]
    @show hook.df[!,"source3_q_inv"][total_steps]
    =#
    @test freq_Δ_1/dP_1 ≈ 1 atol = 0.001
    @test vrms_Δ_1/dQ_1 ≈ 1 atol = 0.01

    @test freq_Δ_2/dP_2 ≈ 1 atol = 0.001
    @test vrms_Δ_2/dQ_2 ≈ 1 atol = 0.001

    @test freq_end_1/env.nc.parameters["grid"]["f_grid"] > 0.99
    @test freq_end_1/env.nc.parameters["grid"]["f_grid"] < 1.01

    @test freq_end_2/env.nc.parameters["grid"]["f_grid"] > 0.99
    @test freq_end_2/env.nc.parameters["grid"]["f_grid"] < 1.01

    @test vrms_end_1/env.nc.parameters["grid"]["v_rms"] > 0.99
    @test vrms_end_1/env.nc.parameters["grid"]["v_rms"] < 1.01

    @test vrms_end_2/env.nc.parameters["grid"]["v_rms"] > 0.99
    @test vrms_end_2/env.nc.parameters["grid"]["v_rms"] < 1.01

    @test hook.df[!,"source3_p_inv"][total_steps]/Source.pq0_set[3, 1] ≈ 1 atol = 0.001
    @test hook.df[!,"source3_q_inv"][total_steps]/Source.pq0_set[3, 2] ≈ 1 atol = 0.001

    return nothing
end
