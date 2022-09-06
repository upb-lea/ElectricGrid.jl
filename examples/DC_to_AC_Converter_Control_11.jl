# activate .
# instantiate
# add OrdinaryDiffEq@5.71.2 SciMLBase@1.47.0 DiffEqBase@6.84.0 ReinforcementLearning@0.10.1

using DrWatson
@quickactivate ("dare")

using ControlSystems
using Plots
using ReinforcementLearning
using LinearAlgebra
using FFTW
using IntervalSets
using CUDA

include(srcdir("Classical_Control.jl"))
include(srcdir("Power_System_Theory.jl"))

include(srcdir("data_hook.jl"))
include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))

include(srcdir("Classical_Control_Plots.jl"))

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

#_______________________________________________________________________________
# Parameters - Time simulation
Timestep = 100 #time step in μs ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_final = 0.4 #time in seconds, total simulation run time

ts = Timestep*1e-6
t = 0:ts:t_final # time

fs = 1/ts # Hz, Sampling frequency of controller ~ 15 kHz -> 50kHz

#_______________________________________________________________________________
# Setting up the Sources

# srated = env.nc.parameters["source"][num_source_to_control]["pwr"]

num_sources = 2

Source = Source_Controller(t_final, fs, num_sources, delay = 1)

Mode_Keys = collect(keys(Source.Modes))

#=
    1 -> "Voltage Control Mode"
    2 -> "PQ Control Mode"
    3 -> "Droop Control Mode"
    4 -> "Synchronverter Mode"
    5 -> "Swing Mode"
=#

Source_Initialiser(Source, Mode_Keys[4], num_source = 1, Srated = 150e3)
Source_Initialiser(Source, Mode_Keys[4], num_source = 2, Srated = 150e3)

#_______________________________________________________________________________
# State space representation

CM = [ 0. 0. 1.
        0. 0. 2
        -1. -2. 0.]

#= CM = [0. 1.
   -1. 0.] =#

parameters = Dict()
source_list = []
source = Dict()

source["fltr"] = "LC"
source["R1"] = Source.Rf[1]
source["R_C"] = 0.0006
source["L1"] = Source.Lf[1]
source["C"] = Source.Cf[1]

push!(source_list, source, source)
#push!(source_list, source)

load_list = []
load = Dict()

R1, L, X, Z = Load_Impedance_2(100e3, 0.6, 230)
#R2, C, X, Z = Load_Impedance(50e3, -0.9999, 230)
load["impedance"] = "RL"
load["R"] = R1
load["L"] = L
#load["C"] = 0.0000000001;
push!(load_list, load)

cable_list = []

# Network Cable Impedances
l = 0.1 # length in km
cable = Dict()
cable["R"] = 0.208*l # Ω, line resistance 0.722#
cable["L"] = 0.00025*l # H, line inductance 0.264e-3#
cable["C"] = 0.4e-6*l # 0.4e-6#
push!(cable_list, cable, cable)
#push!(cable_list, cable)

parameters["source"] = source_list
parameters["cable"] = cable_list
parameters["load"] = load_list;
parameters["grid"] = Dict("fs" => fs, "phase" => 3, "v_rms" => 230);

# Define environment
env = SimEnv(reward_function = reward,  v_dc = 1, ts = ts, use_gpu = false, 
CM = CM, num_sources = num_sources, num_loads = 1, 
parameters = parameters, maxsteps = Source.N_cntr - 1)

# find the indices in the state vector that correspond to the inverters
Collect_IDs(env, Source)

Animo = Classical_Policy(action_space = action_space(env), Source = Source)

#%% Starting time simulation
println("\nHere we go.\n")
reset!(env)
# reset source as well.
# current limitation for PQ mode

plt_state_ids = ["i_f1_a", "i_f1_b", "i_f1_c"] 
#plt_state_ids = ["u_1_a", "u_1_b", "u_1_c", "u_2_a", "u_2_b", "u_2_c"] 
#plt_action_ids = ["u_v1_a", "u_v1_b", "u_v1_c"]
#plt_action_ids = ["u_v1_a", "u_v1_b", "u_v1_c", "u_v2_a", "u_v2_b", "u_v2_c"]
hook = DataHook(collect_state_ids = plt_state_ids#= , collect_action_ids = plt_action_ids =#)

run(Animo, env, StopAfterEpisode(1), hook)

#= @time begin

    println("Progress : 0.0 %")

    for i in 1:Source.N_cntr-1

        # Progress Bar
        if i > 1 && floor((10*t[i]/t_final)) != floor((10*t[i - 1]/t_final))
            flush(stdout)
            println("Progress : ", 10*floor((10*t[i]/t_final)), " %")
        end

        # Control System _______________________________________________________

        if t[i] > t_final/2
            nm_src = 1 # changing the power set points of the 2nd source
            Animo.Source.pq0_set[nm_src, 1] = -0e3 # W, Real Power
            #Animo.Source.pq0_set[nm_src, 2] = 30e3 # VAi, Imaginary Power
        end

        action = Animo(env)

        #-----------------------------------------------------------------------

        # System Dynamics ______________________________________________________

        env(action)
    end

    println("Progress : 100.0 %\n")
end =#

#%% Plots

Plot_I_dq0(0, 5000, Animo.Source, num_source = 1)

Plot_V_dq0(0, 5000, Animo.Source, num_source = 1)

Inst_Vout_Vref(5, 20, Animo.Source, env, num_source = 1)
Inst_Vout_Vref(5, 20, Animo.Source, env, num_source = 2)

Inst_Iout_Iref(10, 20, Animo.Source, env, num_source = 1)
Inst_Iout_Iref(10, 20, Animo.Source, env, num_source = 2)

Plot_PLL(0, 500, Animo.Source, env, num_source = 1, ph = 1)

Plot_Irms(0, 5000, Animo.Source, num_source = 1)

Plot_Vrms(0, 5000, Animo.Source, num_source = 1)
Plot_Vrms(0, 5000, Animo.Source, num_source = 2)

Plot_Real_Imag_Active_Reactive(0, 5000, Animo.Source, num_source = 1)
Plot_Real_Imag_Active_Reactive(0, 5000, Animo.Source, num_source = 2)

#Plot_fft(0, 1, Env, Source, num_node = 2, num_source = 2)

# Save plots
#_______________________________________________________________________________

#=
savefig(p_I_dq0, "p_I_dq0.png")
savefig(p_V_dq0, "p_V_dq0.png")
savefig(p_v_cntr, "p_v_cntr.png")
savefig(p_i_cntr, "p_i_cntr.png")
savefig(p_pll, "p_pll.png")
savefig(p_p, "p_p.png")
savefig(p_v_rms_ang, "p_v_rms_ang.png")
savefig(p_i_rms_ang, "p_i_rms_ang.png")
savefig(p_p_real_imag_act_react, "p_p_real_imag_act_react.png")
savefig(p_fft, "p_fft.png")
=#

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
