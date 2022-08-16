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
Timestep = 100 #time step in μs ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_final = 2 #time in seconds, total simulation run time

#_______________________________________________________________________________
# Environment Calcs
Nps = (1/(Timestep*1e-6)) # time intervals
μps = 1/Nps # time step
t = 0:μps:t_final # time

f_cntr = 1/μps # Hz, Sampling frequency of controller ~ 15 kHz -> 50kHz

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
Source_Initialiser(Source, "PQ Control Mode", num_source = 2, Prated = 100e3, Qrated = 20e3)

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

#=
    CM = [ 0  0  1
        0  0  2
        -1 -2  0]

    power_grid = NodeConstructor(num_sources = 2, num_loads = 1, S2S_p = 1, S2L_p = 1)
    A, B, C, D = get_sys(power_grid)
=#

Env = Environment(t_final, μps, A, B, C, D, num_sources, num_loads)
Env.V_poc_loc = [3 4; 11 12; 19 20] # ID's at which nodes the are sources located
Env.I_poc_loc = [5 6; 13 14; 21 22]
Env.I_inv_loc = [1 2; 9 10; 17 18]

ns = length(A[1,:]) # get num of states
x0 = [0.0 for i = 1:ns] # initial conditions
env = SimEnv(A = A, B = B2, C = C, D = D2, x0 = x0)

reset!(env)

#%% Starting time simulation
println("\nHere we go.\n")

@time begin

    println("Progress : 0.0 %")

    for i in 1:Env.N-1

        # Progress Bar
        if i > 1 && floor((10*t[i]/t_final)) != floor((10*t[i - 1]/t_final))
            flush(stdout)
            println("Progress : ", 10*floor((10*t[i]/t_final)), " %")
        end

        # Environment __________________________________________________________

        Measurements(Env)

        # Control System _______________________________________________________

        if t[i] > t_final/2
            Source.pq0_set[2,:] = [-50e3; 10e3; 0]
        end

        #-----------------------------------------------------------------------

        Action = Classical_Policy(Source, Env)

        #=
            num_srce = 1

            Source_Interface(Env, Source, active = 0, fc = 2000)

            #Swing_Mode(Source, num_srce, ramp = 1)
            #Voltage_Control_Mode(Source, num_srce, ramp = 1)
            #PQ_Control_Mode(Source, num_srce, pq0)
            Droop_Control_Mode(Source, num_srce, ramp = 1, t_end = 1.0)
            #Synchronverter_Mode(Source, num_srce)

            Phase_Locked_Loop_3ph(Source, num_srce)
            #Phase_Locked_Loop_1ph(Source, num_srce, i, Kp = 0.5, Ki = 5, ph = 1)
            #Phase_Locked_Loop_1ph(Source, num_srce, i, Kp = 0.5, Ki = 5, ph = 2)
            #Phase_Locked_Loop_1ph(Source, num_srce, i, Kp = 0.5, Ki = 5, ph = 3)

            #-----------------------------------------------------------------------
            num_srce = 2

            Phase_Locked_Loop_3ph(Source, num_srce)
            #Voltage_Control_Mode(Source, num_srce, ramp = 1)
            #PQ_Control_Mode(Source, num_srce, pq0)
            #Droop_Control_Mode(Source, num_srce, ramp = 1, t_end = 1.0)
            Synchronverter_Mode(Source, num_srce, pq0_ref = pq0, ramp = 1, t_end = 1.0)
            #Self_Synchronverter_Mode(Source, num_srce, pq0_ref = pq0)

            Action = Env_Interface(Env, Source, num_srce)
        =#

        # System Dynamics ______________________________________________________

        Evolution(Env, Action)

    end

    println("Progress : 100.0 %\n")
end

#%% Plots

Plot_I_dq0(0, 5000, Source, Env, num_source = 2)

Plot_V_dq0(0, 5000, Source, Env, num_source = 2)

Inst_Vout_Vref(0, 50, Source, Env, num_node = 2, num_source = 2)

Inst_Iout_Iref(0, 50, Source, Env, num_source = 2, num_node = 2)

Plot_PLL(0, 500, Source, Env, num_source = 2, ph = 2)

Plot_Irms(0, 5000, Env, num_node = 2)

Plot_Vrms(0, 5000, Env, Source, num_node = 1, num_source = 1)
Plot_Vrms(0, 5000, Env, Source, num_node = 2, num_source = 2)

Plot_Real_Imag_Active_Reactive(0, 5000, Env, Source, num_node = 1, num_source = 1)
Plot_Real_Imag_Active_Reactive(0, 5000, Env, Source, num_node = 2, num_source = 2)

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
