using Test
using Dare
using MAT
using ReinforcementLearning




@testset "env_1source_1load" begin

    t_final = 0.0003 #time in seconds, total simulation run time
    ts = 1e-6
    t = 0:ts:t_final # time

    fs = 1/ts # Hz, Sampling frequency of controller ~ 15 kHz < fs < 50kHz

    CM = [0. 1.
    -1. 0.] 

    parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("fltr"=>"LC", "pwr"=>10000e3, "control_type" =>"classic", "source_type"=>"ideal", "mode" => 8,  "R1"=>1.1e-3, "L1"=>1e-3, "C"=>1e-3, "R_C"=>7e-3, "vdc"=>800, "v_limit"=>10000, "i_limit"=>10000)
                        ],
        "load"   => Any[
                        Dict{Any, Any}("impedance"=>"RLC", "R"=>100, "L"=>1e-2, "C"=>1e-2, "pf"=>0.8, "v_limit"=>10000, "i_limit"=>10000) 
                        ],
        "cable" => Any[
                        Dict{Any, Any}("len"=>1, "R"=>1e-3, "L"=>1e-4, "C"=>1e-4, "i_limit"=>10000),
                        ],
        "grid"   => Dict{Any, Any}("fs"=>50.0, "phase"=>3, "v_rms"=>230, "f_grid" => 50, "ramp_end"=>0.0)
    )


    env = SimEnv(ts = ts, use_gpu = false, CM = [0 1;1 0], num_sources = 1, num_loads = 1, verbosity = 0,parameters = parameters, maxsteps = length(t), action_delay = 1)

    plt_state_ids = ["source1_i_L1_a", "source1_v_C_filt_a",  "source1_v_C_cables_a", "cable1_i_L_a", "load1_v_C_total_a", "load1_i_L_a"]               
    plt_action_ids = ["source1_u_a"]
    hook = DataHook(collect_state_ids = plt_state_ids, collect_action_ids = plt_action_ids)

    Power_System_Dynamics(env, hook)

    idx_end = 300
    test_state_ids = ["next_state_source1_i_L1_a", "next_state_source1_v_C_filt_a", "next_state_source1_v_C_cables_a", "next_state_cable1_i_L_a", "next_state_load1_v_C_total_a", "next_state_load1_i_L_a"]  
    X_dare = Matrix(hook.df[!, test_state_ids][1:idx_end,:])
    X_malab = matread("./test/env_test_state_1source_1load1e6.mat")
        
    #=println("MATLAB:")
    display(X_malab["X_matlab"][1:idx_end,:])
    println()
    println("DARE:")
    display(X_dare)
    println() =#

    @test X_dare≈X_malab["X_matlab"][1:idx_end,:] atol=0.001   # 1e-6
end


@testset "env_2source_1load" begin

    CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.]


    parameters = Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}("fltr"=>"LC", "pwr"=>10000e3, "control_type" =>"classic", "source_type"=>"ideal", "mode" => 8,  "R1"=>1.1e-3, "L1"=>1e-3, "C"=>1e-3, "R_C"=>7e-3, "vdc"=>800, "v_limit"=>10000, "i_limit"=>10000),
                    Dict{Any, Any}("fltr"=>"LCL", "pwr"=>10000e3, "control_type" =>"classic", "source_type"=>"ideal", "mode" => 8,  "R1"=>1.1e-3, "L1"=>1e-3, "C"=>1e-3, "R_C"=>7e-3, "L2"=>1e-3, "R2"=>1.1e-3,  "vdc"=>800, "v_limit"=>10000, "i_limit"=>10000)
                    ],
    "load"   => Any[
                    Dict{Any, Any}("impedance"=>"RLC", "R"=>100, "L"=>1e-2, "C"=>1e-2, "pf"=>0.8, "v_limit"=>10000, "i_limit"=>10000) 
                    ],
    "cable" => Any[
                    Dict{Any, Any}("len"=>1, "R"=>1e-3, "L"=>1e-4, "C"=>1e-4, "i_limit"=>10000),
                    Dict{Any, Any}("len"=>1, "R"=>1e-3, "L"=>1e-4, "C"=>1e-4, "i_limit"=>10000),
                    ],
    "grid"   => Dict{Any, Any}("fs"=>50.0, "phase"=>3, "v_rms"=>230, "f_grid" => 50, "ramp_end"=>0.0)      
    )

    env = SimEnv(ts = 1e-6, use_gpu = false, CM = CM, num_sources = 2, num_loads = 1, parameters = parameters, maxsteps = 300, action_delay = 1, verbosity = 0)

    hook = DataHook(collect_sources = [1,2], collect_loads = [1], collect_cables = [1,2])

    Power_System_Dynamics(env, hook)

    idx_end = 300
    test_state_ids = ["next_state_source1_i_L1_a", "next_state_source1_v_C_filt_a", "next_state_source1_v_C_cables_a", "next_state_cable1_i_L_a", "next_state_load1_v_C_total_a", "next_state_load1_i_L_a", "next_state_source2_i_L1_a", "next_state_source2_v_C_filt_a", "next_state_source2_i_L2_a","next_state_source2_v_C_cables_a", "next_state_cable2_i_L_a",]  
    X_dare = Matrix(hook.df[!,test_state_ids][1:idx_end,:])
    X_malab = matread("./test/env_test_state_2source_1load1e6.mat")

    @test X_dare≈X_malab["X_matlab"][1:idx_end,:] atol=.001   # 1e-6

end

@testset "env_2source" begin
    CM = [ 0. 1.
        -1. 0.]

    parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => 200e3, "mode" => 8, "fltr" => "L", "L1" => 1e-3, "R1" => 1.1e-3, "i_limit"=>10e6),
                        Dict{Any, Any}("pwr" => 200e3, "mode" => 8, "fltr" => "L", "L1" => 1e-3, "R1" => 1.1e-3, "i_limit"=>10e6),
                        ],
        "cable"   => Any[
                        Dict{Any, Any}("R" => 1e-3, "L" => 1e-4, "C" => 1e-4, "i_limit" => 10e4,),
                        ],
        "grid" => Dict{Any, Any}("ramp_end" => 0.0)
    )
    env = SimEnv(ts = 1e-6, CM = CM, parameters = parameters, verbosity = 0, maxsteps = 600, action_delay = 1)

    hook = DataHook(collect_sources  = [1 2],
                    collect_cables = [1])
    
    Power_System_Dynamics(env, hook)

    test_state_ids = ["next_state_source1_i_L1_a", "next_state_source2_i_L1_a",   "next_state_source1_v_C_cables_a", "next_state_cable1_i_L_a", "next_state_source2_v_C_cables_a"]
    idx_end = 300
    X_dare = Matrix(hook.df[!, test_state_ids][1:idx_end,:])
    X_malab = matread("./test/env_test_state_2source1e6.mat")

    @test X_dare≈X_malab["X_matlab"][1:idx_end,:] atol=0.01   # 1e-6
end;

@testset "env_3source_2load" begin
    CM = [0.  0.  0. 1. 2.
        0.  0.  0. 3. 4.
        0.  0.  0. 5. 6.
       -1. -3. -5. 0. 0.
       -2. -4. -6. 0. 0.]
    # Filter
    R_1 = 1.1e-3;
    L_1 = 1e-4;
    R_c = 7e-3;
    C_1 = 1e-5;


    C_b = 1e-8/2;
    L_b = 1e-6;
    R_b = 5;

    R_l = 1e6;
    C_l = 1e-4;
    L_l = 1e-5;
       

    parameters = Dict{Any, Any}(
                                "source" => Any[
                                                Dict{Any, Any}("fltr"=>"LC", "pwr"=>10000e3, "control_type" =>"classic", "source_type"=>"ideal", "mode" => 8,  "R1"=>R_1, "L1"=>L_1, "C"=>C_1, "R_C"=>R_c, "vdc"=>800, "v_limit"=>10000, "i_limit"=>10e8),
                                                Dict{Any, Any}("fltr"=>"LCL", "pwr"=>10000e3, "control_type" =>"classic", "source_type"=>"ideal", "mode" => 8,  "R1"=>R_1, "L1"=>L_1, "C"=>C_1, "R_C"=>R_c, "L2"=>L_1, "R2"=>R_1,  "vdc"=>800, "v_limit"=>10000, "i_limit"=>10e8),
                                                Dict{Any, Any}("pwr" => 10000e3, "mode" => 8, "fltr" => "L", "L1" => L_1, "R1" => R_1, "i_limit"=>10e8),
                                                ],
                                "load"   => Any[
                                                Dict{Any, Any}("impedance"=>"RLC", "R"=>R_l, "L"=>L_l, "C"=>C_l, "pf"=>0.8, "v_limit"=>10e8, "i_limit"=>10e8), 
                                                Dict{Any, Any}("impedance"=>"R", "R"=>R_l, "v_limit"=>10e4, "i_limit"=>10e4) 
                                                ],
                                "cable" => Any[
                                                Dict{Any, Any}("len"=>1, "R"=>R_b, "L"=>L_b, "C"=>C_b, "i_limit"=>10e8),
                                                Dict{Any, Any}("len"=>1, "R"=>R_b, "L"=>L_b, "C"=>C_b, "i_limit"=>10e8),
                                                Dict{Any, Any}("len"=>1, "R"=>R_b, "L"=>L_b, "C"=>C_b, "i_limit"=>10e8),
                                                Dict{Any, Any}("len"=>1, "R"=>R_b, "L"=>L_b, "C"=>C_b, "i_limit"=>10e8),
                                                Dict{Any, Any}("len"=>1, "R"=>R_b, "L"=>L_b, "C"=>C_b, "i_limit"=>10e8),
                                                Dict{Any, Any}("len"=>1, "R"=>R_b, "L"=>L_b, "C"=>C_b, "i_limit"=>10e8),
                                                ],
                                "grid"   => Dict{Any, Any}("fs"=>50.0, "phase"=>3, "v_rms"=>230, "f_grid" => 50, "ramp_end"=>0.0)      
                                )
    #__________________________________________________________________
    # Defining the environment

    env = SimEnv(ts = 1e-6, CM = CM, parameters = parameters, maxsteps=300,  verbosity = 0)

    #_______________________________________________________________________________
    # Setting up data hooks

    hook = DataHook(collect_sources  = [1 2 3],
                    collect_loads = [1], 
                    collect_cables = [1, 5])

    #_______________________________________________________________________________
    # Running the Time Simulation

    Power_System_Dynamics(env, hook)

    #_______________________________________________________________________________
    # Plotting
    test_state_ids = ["next_state_source1_i_L1_a", "next_state_source1_v_C_filt_a", "next_state_source2_v_C_cables_a", "next_state_source3_i_L1_a", "next_state_load1_i_L_a", "next_state_source2_i_L1_a", "next_state_cable1_i_L_a", "load1_v_C_total_a",   "next_state_cable5_i_L_a"]
    
    idx_end = 300
    X_dare = Matrix(hook.df[!,test_state_ids][1:idx_end,:])
    X_malab = matread("./test/env_test_state_3source_2load1e6.mat")

    println("MATLAB:")
    display(X_malab["X_matlab"][1:idx_end,:])
    println()
    println("DARE:")
    display(X_dare)
    println()

    @test X_dare≈X_malab["X_matlab"][1:idx_end,:] atol=14 # such high, since the steady state current is very high, so talking about a toleranz below 1 % of max values
end


