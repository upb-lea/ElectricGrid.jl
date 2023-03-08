using JEG;


print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n\n")

#_______________________________________________________________________________
# Parameters - Time simulation
t_final = 0.03 #time in seconds, total simulation run time
ts = 1e-4
t = 0:ts:t_final # time

fs = 1/ts # Hz, Sampling frequency of controller ~ 15 kHz < fs < 50kHz

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


env = ElectricGridEnv(ts = ts, use_gpu = false, CM = CM, num_sources = 2, num_loads = 1, parameters = parameters, maxsteps = length(t), action_delay = 1)

#_______________________________________________________________________________
#%% Setting up data hooks


plt_state_ids = ["source1_v_C_filt_a", "source1_i_L1_a", "source1_v_C_cables_a", "cable1_i_L_a", "load1_v_C_total_a", "load1_i_L_a"]               
plt_action_ids = ["source1_u_a", "source2_u_a",]
hook = DataHook(collect_sources = [1,2])

#_______________________________________________________________________________
# Starting time simulation
num_eps = 1
#RLBase.run(ma, env, StopAfterEpisode(num_eps), hook);

Power_System_Dynamics(env, hook)

RenderHookResults(; hook = hook, states_to_plot = ["source1_v_C_filt_a","source2_v_C_filt_a",], actions_to_plot = plt_action_ids)

idx_end = 10
#=
function collect_matrix(plt_state_ids, idx_end)
    X_JEG2 = []
    for s_id in plt_state_ids
        if s_id == plt_state_ids[1]
            X_JEG2 = [hook.df[!,s_id][2:idx_end,:]']
        else
            hcat(X_JEG2, hook.df[!,s_id][2:idx_end,:]')
        end
    end
    return X_JEG2
end
=#


X_JEG2 = []
for s_id in plt_state_ids
    if s_id == plt_state_ids[1]
        println("Hallo")
        push!(X_JEG2, [hook.df[!,s_id][2:idx_end,:]'])
        println(X_JEG2)
    else
        hcat(X_JEG2, hook.df[!,s_id][2:idx_end,:]')
    end
end


#X_d = collect_matrix(plt_state_ids, idx_end)

println()
println("JEG:")
display(X_JEG2)
println()  


print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
