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
include(srcdir("Classical_Control_Plots.jl"))
include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))
include(srcdir("sin_policy.jl"))

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

#_______________________________________________________________________________
# Parameters - Time simulation
Timestep = 75 #time step in μs ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_final = 0.4 #time in seconds, total simulation run time

#_______________________________________________________________________________
# Environment Calcs

Ts = Timestep*1e-6
t = 0:Ts:t_final # time

f_cntr = 1/Ts # Hz, Sampling frequency of controller ~ 15 kHz -> 50kHz

#_______________________________________________________________________________
# Setting up the Sources

num_sources = 2

Source = Source_Controller(t_final, f_cntr, num_sources, delay = 1)

#=
    Typical values for the frequency droop are a 100% increase in power for a
    frequency decrease between 3% and 5% (from nominal values)
=#

Source.Δfmax = 0.03*50/100 # The drop in frequency, Hz, which will cause a 100% increase in active power
Source.ΔEmax = 0.05*230/100 # The drop in rms voltage, which will cause a 100% decrease in reactive power
τ = 1.0

Source.τv = τ # time constant of the voltage loop
Source.τf = τ

Mode_Keys = collect(keys(Source.Modes))

#=
    1 -> "Voltage Control Mode"
    2 -> "PQ Control Mode"
    3 -> "Droop Control Mode"
    4 -> "Synchronverter Mode"
    5 -> "Swing Mode"
=#

Source_Initialiser(Source, Mode_Keys[5], num_source = 1, Srated = 150e3)
Source_Initialiser(Source, Mode_Keys[5], num_source = 2, Srated = 150e3)

#_______________________________________________________________________________
# Circuit Elements Calcs

#_______________________________________________________________________________
# State space representation

CM = [ 0. 0. 1.
        0. 0. 2
        -1. -2. 0.]

parameters = Dict()
source_list = []
source = Dict()

source["fltr"] = "LC"
source["R1"] = 0.4 #Source.Rf[1]
source["R_C"] = 0.0006
source["L1"] = 2.3e-3#0.001*Source.Lf[1]
source["C"] = 1e-6 #0.001*Source.Cf[1]

push!(source_list, source, source)

load_list = []
load = Dict()

R1, L, X, Z = Load_Impedance(50e3, 0.9999, 230)
R2, C, X, Z = Load_Impedance(50e3, -0.9999, 230)
load["impedance"] = "R"
load["R"] = 1000000;
#load["L"] = L;
#load["C"] = 0.1;
push!(load_list, load);

cable_list = []

# Network Cable Impedances
l = 0.01 # length in km
cable = Dict()
cable["R"] = 0.722# 0.222*l # Ω, line resistance 0.722#
cable["L"] = 0.264e-3# 0.0024*l # H, line inductance 0.264e-3#
cable["C"] = 0.4e-6# 7.3077e-9*l # 0.4e-6#
push!(cable_list, cable, cable)

parameters["source"] = source_list
parameters["cable"] = cable_list
parameters["load"] = load_list;
parameters["grid"] = Dict("fs" => f_cntr, "phase" => 3, "v_rms" => 230);

# Define environment
env = SimEnv(reward_function = reward,  v_dc = 1000, ts = Ts, use_gpu = false, 
CM = CM, num_sources = num_sources, num_loads = 1, parameters = parameters, maxsteps = Source.N_cntr - 1)

A, B, C, D = get_sys(env.nc)

ns = length(A[1,:]) # get num of states
ni = length(B[1,:]) # get num of inputs

# find the indices in the state vector that correspond to the inverters
Collect_IDs(env, Source)
# Source.V_poc_loc = [3 6; 13 16; 23 26] # ID's at which nodes the sources are located
# Source.I_poc_loc = [1 4; 11 14; 21 24]
# Source.I_inv_loc = [1 4; 11 14; 21 24]

Animo = Classical_Policy(action_space = action_space(env), Source = Source)

#%% Starting time simulation
println("\nHere we go.\n")
reset!(env)
# reset source as well.
# current limitation for PQ mode

input_action_a = zeros(Source.N_cntr-1)
input_action_b = zeros(Source.N_cntr-1)
input_action_c = zeros(Source.N_cntr-1)

env_action_a = zeros(Source.N_cntr-1)
env_action_b = zeros(Source.N_cntr-1)
env_action_c = zeros(Source.N_cntr-1)

vout_a = zeros(Source.N_cntr-1)
vout_b = zeros(Source.N_cntr-1)
vout_c = zeros(Source.N_cntr-1)

@time begin

    println("Progress : 0.0 %")

    for i in 1:Source.N_cntr-1

        # Progress Bar
        if i > 1 && floor((10*t[i]/t_final)) != floor((10*t[i - 1]/t_final))
            flush(stdout)
            println("Progress : ", 10*floor((10*t[i]/t_final)), " %")
        end

        # Control System _______________________________________________________

        if t[i] > t_final/2
            #num_source = 2 # changing the power set points of the 2nd source
            #policy.Source.pq0_set[num_source, 1] = -50e3 # W, Real Power
            #policy.Source.pq0_set[num_source, 2] = 30e3 # VAi, Imaginary Power
        end

        #-----------------------------------------------------------------------

        # System Dynamics ______________________________________________________

        #action = Animo(env)

        u = [1*sin.(50*2*pi*t[i]) for j = 1:3]
        action = vcat(u,u)

        env(action)

        s = 1
        input_action_a[i] = action[s + Source.num_sources*(1 - 1)]
        input_action_b[i] = action[s + Source.num_sources*(2 - 1)]
        input_action_c[i] = action[s + Source.num_sources*(3 - 1)]

        env_action_a[i] = env.action[1 + 3*(s - 1)]
        env_action_b[i] = env.action[2 + 3*(s - 1)]
        env_action_c[i] = env.action[3 + 3*(s - 1)]

        vout_a[i] = env.state[Source.V_poc_loc[1, s]]
        vout_b[i] = env.state[Source.V_poc_loc[2, s]]
        vout_c[i] = env.state[Source.V_poc_loc[3, s]]

    end

    println("Progress : 100.0 %\n")
end

#%% Plots

num_source = 1
T_plot_start = 0
T_plot_end = 10

if T_plot_end > t_final*Source.fsys
    T_plot_end = t_final*Source.fsys
end

Nps = Source.f_cntr
N_plot_start = convert(Int64, round((T_plot_start/Source.fsys  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/Source.fsys  - 1/Nps)*Nps))
range = N_plot_start:N_plot_end

v_out = plot(t[range], vout_a[range], label = "a",
            xlabel = "time", ylabel = "V poc", title = "Env.State")
v_out = plot!(t[range], vout_b[range], label = "b")
v_out = plot!(t[range], vout_c[range], label = "c")
display(v_out)
 
u = plot(t[range], input_action_a[range], label = "a",
            xlabel = "time", ylabel = "V inv", title = "Policy Action")
u = plot!(t[range], input_action_b[range], label = "b")
u = plot!(t[range], input_action_c[range], label = "c")
display(u)

u = plot(t[range], env_action_a[range], label = "a",
            xlabel = "time", ylabel = "V inv", title = "Env.Action")
u = plot!(t[range], env_action_b[range], label = "b")
u = plot!(t[range], env_action_c[range], label = "c")
display(u)

#= v_out = plot(t[range], Source.V_filt_poc[num_source, 1, range], label = "a",
            xlabel = "time", ylabel = "V poc", title = "Inverter Output")
v_out = plot!(t[range], Source.V_filt_poc[num_source, 2, range], label = "b")
v_out = plot!(t[range], Source.V_filt_poc[num_source, 3, range], label = "c")
display(v_out)
 
u = plot(t[range], Animo.Source.Vd_abc_new[num_source, 1, range], label = "a",
            xlabel = "time", ylabel = "V inv", title = "Action")
u = plot!(t[range], Animo.Source.Vd_abc_new[num_source, 2, range], label = "b")
u = plot!(t[range], Animo.Source.Vd_abc_new[num_source, 3, range], label = "c")
display(u) =#

#Plot_I_dq0(0, 5000, Animo.Source, num_source = 2)

#Plot_V_dq0(0, 5000, Animo.Source, num_source = 2)

#Inst_Vout_Vref(0, 20, Animo.Source, env, num_source = 1)
#Inst_Vout_Vref(0, 20, Animo.Source, env, num_source = 2)

#Inst_Iout_Iref(0, 20, Animo.Source, env, num_source = 1)
#Inst_Iout_Iref(0, 20, Animo.Source, env, num_source = 2)

#Plot_PLL(0, 500, Animo.Source, env, num_source = 2, ph = 2)

#Plot_Irms(0, 5000, Animo.Source, num_source = 2)

#Plot_Vrms(0, 5000, Animo.Source, num_source = 1)
#Plot_Vrms(0, 5000, Animo.Source, num_source = 2)

#Plot_Real_Imag_Active_Reactive(0, 5000, Animo.Source, num_source = 1)
#Plot_Real_Imag_Active_Reactive(0, 5000, Animo.Source, num_source = 2)

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
