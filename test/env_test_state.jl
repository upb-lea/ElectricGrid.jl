using Test
using Dare
using MAT
using ReinforcementLearning

vars = matread("./test/env_test_state_1source_1load1e6.mat")


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


    env = SimEnv(ts = ts, use_gpu = false, CM = [0 1;1 0], num_sources = 1, num_loads = 1, parameters = parameters, maxsteps = length(t), action_delay = 1)

    state_ids = get_state_ids(env.nc)
    action_ids = get_action_ids(env.nc)

    #_______________________________________________________________________________
    # Setting up the Classical Sources

    Animo = NamedPolicy("classic", Classical_Policy(env))

    state_ids_classic = Animo.policy.state_ids
    action_ids_classic = Animo.policy.action_ids

    Multi_Agents = Dict()
    Multi_Agent_list = []

    polc = Dict()

    polc["policy"] = Animo
    polc["state_ids"] = state_ids_classic
    polc["action_ids"] = action_ids_classic

    Multi_Agents[nameof(Animo)] = polc

    ma = MultiAgentGridController(Multi_Agents, action_ids)

    agentname = "agent"

    #_______________________________________________________________________________
    #%% Setting up data hooks


    plt_state_ids = ["source1_v_C_filt_a", "source1_i_L1_a", "source1_v_C_cables_a", "cable1_i_L_a", "load1_v_C_total_a", "load1_i_L_a"]               
    plt_action_ids = ["source1_u_a"]#"source1_u_a", "u_v1_b", "u_v1_c"]
    hook = DataHook(collect_state_ids = plt_state_ids, collect_action_ids = plt_action_ids,  collect_sources = [],
    collect_cables = [], collect_vrms_ids = [1], collect_irms_ids = [], collect_pq_ids = [], collect_vdq_ids = [], collect_idq_ids = [],
    save_best_NNA = false, collect_reference = false, plot_rewards = false)

    #_______________________________________________________________________________
    # Starting time simulation
    num_eps = 1
    RLBase.run(ma, env, StopAfterEpisode(num_eps), hook); # this can be replaced by the method called: Power_System_Dynamics(env, hook)


    for eps in 1:num_eps

        plot_hook_results(; hook = hook, states_to_plot = plt_state_ids, actions_to_plot = plt_action_ids, episode = eps, 
        pq_to_plot = [], vrms_to_plot = [], irms_to_plot = [], vdq_to_plot = [], idq_to_plot = [])
    end
    display(hook.df)
    println(env.state_ids)
    idx_end = 300
    X_dare = [hook.df[!,"source1_i_L1_a"][2:idx_end,:] hook.df[!, "source1_v_C_filt_a"][2:idx_end,:] hook.df[!,"source1_v_C_cables_a"][2:idx_end,:] hook.df[!,"cable1_i_L_a"][2:idx_end,:] -hook.df[!,"load1_v_C_total_a"][2:idx_end,:] -hook.df[!,"load1_i_L_a"][2:idx_end,:]]

    # TODO: why do we have 2 zero lines in X_dare?
    # TODO: Why is load1_i_L_a in X_dare off?

    println("MATLAB:")
    display(vars["X_matlab"][1:idx_end-1,:])
    println()
    println("DARE:")
    display(X_dare)
    println()   

    #@test X_dare≈vars["X_matlab"][1:idx_end-1,:] atol=12  # 1e-4
    #@test X_dare≈vars["X_matlab"][1:idx_end-1,:] atol=0.1   # 1e-5
    @test X_dare≈vars["X_matlab"][1:idx_end-1,:] atol=0.001   # 1e-6
end