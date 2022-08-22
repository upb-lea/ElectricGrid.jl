# activate .
# instantiate

using DrWatson
@quickactivate "dare"

using Plots
using LinearAlgebra
using FFTW
using ControlSystems
using ReinforcementLearning
using IntervalSets
using CUDA

include(srcdir("Classical_Control.jl"))
include(srcdir("Power_System_Theory.jl"))
include(srcdir("Classical_Control_Plots.jl"))
include(srcdir("nodeconstructor.jl"))
include(srcdir("env.jl"))

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

#_______________________________________________________________________________
# Parameters - Time simulation
Timestep = 75 #time step in μs ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_final = 1.5 #time in seconds, total simulation run time

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

Source_Initialiser(Source, Mode_Keys[3], num_source = 1, Prated = 150e3)
Source_Initialiser(Source, Mode_Keys[2], num_source = 2, Prated = 100e3)

#_______________________________________________________________________________
# Circuit Elements Calcs

num_loads = 1

# Load 1 Impedance
Vo_rms = 230 # rms output voltage

SL1 = 75e3 # VA, 3-ph Apparent Power
pf1 = 0.6 # power factor
SL2 = 50e3 # VA, 3-ph Apparent Power
pf2 = -0.9999 # power factor

# Network Cable Impedances
l = 0.5 # length in km
Lt1 = 0.0024*l # H, line inductance
Lt2 = 0.0024*l
Rt1 = 0.222*l # Ω, line resistance
Rt2 = 0.222*l

#_______________________________________________________________________________
# State space representation
A, B, C, D, B2, D2 = Two_Sources_One_Load(Source, Vo_rms, SL1, pf1, SL2, pf2, Lt1, Lt2, Rt1, Rt2)

Env = Environment(t_final, Ts, A, B, C, D, num_sources, num_loads)
Source.V_poc_loc = [3 4; 11 12; 19 20] # ID's at which nodes the sources are located
Source.I_poc_loc = [5 6; 13 14; 21 22]
Source.I_inv_loc = [1 2; 9 10; 17 18]

CM = [ 0. 0. 1.
        0. 0. 2
        -1. -2. 0.]

parameters = Dict()
source_list = []
source = Dict()

source["pwr"] = 45000.0
source["v_rip"] = 0.01556109320329396
source["vdc"] = 750
source["i_rip"] = 0.10108821490394984
source["fltr"] = "LC"
source["R1"] = 0.4022094955070556
source["R_C"] = 0.0006447094780419011
source["L1"] = 0.001005523738767639
#source["R2"] = 0.4022094955070556   # needed for LCL
#source["L2"] = 0.001005523738767639
source["C"] = 2.302533850149647e-5;

push!(source_list, source, source);

load_list = []
load = Dict()

load["impedance"] = "RLC"
load["R"] = 30236.0;
load["L"] = 57.042;
load["C"] = 39.18;
push!(load_list, load);

cable_list = []

cable = Dict()
cable["R"] = 6.84059
cable["L"] = 0.00250127
cable["C"] = 3.7898e-6;
push!(cable_list, cable, cable);

parameters["source"] = source_list
parameters["cable"] = cable_list
parameters["load"] = load_list;
parameters["grid"] = Dict("fs" => 10000.0, "phase" => 3, "v_rms" => 230);

#######################################################################################
# Define grid using random initialization
power_grid = NodeConstructor(num_sources = 2, num_loads = 1, S2S_p=1, S2L_p=1, CM = CM, parameters = parameters);
A, B, C, D = get_sys(power_grid)
ns = length(A[1,:]) # get num of states
ni = length(B[1,:]) # get num of inputs
x0 = [0.0 for i = 1:ns]
env = SimEnv(A = A, B = B, C = C, D = D, x0 = x0, state_ids = get_state_ids(power_grid), rewardfunction = reward)
#Collect_IDs(env, Source, 2, 0, 0)

policy = Classical_Policy(action_space = action_space(env), Source = Source)

#%% Starting time simulation
println("\nHere we go.\n")
reset!(env)
# reset source as well.
# current limitation for PQ mode

@time begin

    println("Progress : 0.0 %")

    for i in 1:Env.N-1

        # Progress Bar
        if i > 1 && floor((10*t[i]/t_final)) != floor((10*t[i - 1]/t_final))
            flush(stdout)
            println("Progress : ", 10*floor((10*t[i]/t_final)), " %")
        end

        # Control System _______________________________________________________

        if t[i] > t_final/2
            num_source = 2 # changing the power set points of the 2nd source
            #Source.pq0_set[num_source, 1] = -50e3 # W, Real Power
            Source.pq0_set[num_source, 2] = 30e3 # VAi, Imaginary Power
        end

        #-----------------------------------------------------------------------

        Action = Classical_Control(Source, Env)

        # System Dynamics ______________________________________________________

        #env(Action)
        Evolution(Env, Action)

    end

    println("Progress : 100.0 %\n")
end

#%% Plots

Plot_I_dq0(0, 5000, Source, num_source = 2)

Plot_V_dq0(0, 5000, Source, num_source = 2)

Inst_Vout_Vref(0, 20, Source, Env, num_source = 2)

Inst_Iout_Iref(0, 20, Source, Env, num_source = 2)

Plot_PLL(0, 500, Source, Env, num_source = 2, ph = 2)

Plot_Irms(0, 5000, Source, num_source = 2)

Plot_Vrms(0, 5000, Source, num_source = 1)
Plot_Vrms(0, 5000, Source, num_source = 2)

Plot_Real_Imag_Active_Reactive(0, 5000, Source, num_source = 1)
Plot_Real_Imag_Active_Reactive(0, 5000, Source, num_source = 2)

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

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")
