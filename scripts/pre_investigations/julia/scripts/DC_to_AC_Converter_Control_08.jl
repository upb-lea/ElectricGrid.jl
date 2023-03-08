using Plots
using LinearAlgebra
using FFTW

include("Classical_Control.jl")
include("Power_System_Theory.jl")
include("Classical_Control_Plots.jl")

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

#_______________________________________________________________________________
# Parameters - Time simulation
Timestep = 100 #time step in μs
t_final = 1 #time in seconds, total simulation run time

num_sources = 2
num_loads = 1

# Load 1 Impedance
Vo_rms = 230 # rms output voltage

SL1 = 100e3 # VA, 3-ph Apparent Power
pf1 = 0.6 # power factor
SL2 = 1e3 # VA, 3-ph Apparent Power
pf2 = -0.1 # power factor

#_______________________________________________________________________________
# Environment Calcs
Nps = (1/(Timestep*1e-6)) # time intervals
μps = 1/Nps # time step
t = 0:μps:t_final # time

f_cntr = 1/μps # Hz, Sampling frequency of controller ~ 15 kHz -> 50kHz

Source = Source_Controller(t_final, f_cntr, num_sources)

#_______________________________________________________________________________
# Circuit Elements Calcs

RL1, LL, XL, ZL1 = Load_Impedance(SL1, pf1, Vo_rms)
RL2, CL, XC, ZL2 = Load_Impedance(SL2, pf2, Vo_rms)

# Reference currents
I_ref = Vo_rms/ZL1
I_ref_p = abs(I_ref)*sqrt(2) # peak reference current
I_ref_rms = abs(I_ref) # peak reference current
I_ref_ang = angle(I_ref)

# Inv 1 Source Filter Calcs
Lf1, Cf1, fc1 = Filter_Design(Source, 1, 0.15, 0.01537)
Lf2, Cf2, fc2 = Filter_Design(Source, 2, 0.15, 0.01537)

# Network Cable Impedances
l = 1 # length in km
Lt1 = 0.0024*l
Lt2 = 0.0024*l
Rt1 = 0.222*l
Rt2 = 0.222*l

#_______________________________________________________________________________
# State space representation
#A, B, C, D = Simple_State_Space(0, Lf1, Cf1, LL, RL)
A, B, C, D = Two_Sources_One_Load(Lf1, Lf2, Cf1, Cf2, LL, CL, RL1, RL2, Lt1, Lt2, Rt1, Rt2)

Env = Environment(t_final, μps, A, B, C, D, num_sources, num_loads)
Env.V_poc_loc = [3 4; 11 12; 19 20] # at which nodes are the sources located
Env.I_poc_loc = [5 6; 13 14; 21 22]
Env.I_inv_loc = [1 2; 9 10; 17 18]

Env.Δfmax = 0.03
Env.ΔEmax = 5
Source_Initialiser(Source, Env; num_source = 1, Srated = 10e3)
Source_Initialiser(Source, Env; num_source = 2, Srated = 90e3)

#_______________________________________________________________________________
# Control Reference signals
V_ref_p = sqrt(2)*Vo_rms
θtest = (2*π*50*t).%(2*π)
Source.V_ref[1, 1, :] = V_ref_p*sin.(θtest)
Source.V_ref[1, 2, :] = V_ref_p*sin.(θtest .- 120*π/180)
Source.V_ref[1, 3, :] = V_ref_p*sin.(θtest .+ 120*π/180)

Source.V_ref[2, 1, :] = V_ref_p*sin.(θtest)
Source.V_ref[2, 2, :] = V_ref_p*sin.(θtest .- 120*π/180)
Source.V_ref[2, 3, :] = V_ref_p*sin.(θtest .+ 120*π/180)

Source.I_ref[1, 1, :] = I_ref_p*sin.(θtest .+ I_ref_ang)
Source.I_ref[1, 2, :] = I_ref_p*sin.(θtest .- 120*π/180 .+ I_ref_ang)
Source.I_ref[1, 3, :] = I_ref_p*sin.(θtest .+ 120*π/180 .+ I_ref_ang)

Source.I_ref[2, 1, :] = I_ref_p*sin.(θtest .+ I_ref_ang)
Source.I_ref[2, 2, :] = I_ref_p*sin.(θtest .- 120*π/180 .+ I_ref_ang)
Source.I_ref[2, 3, :] = I_ref_p*sin.(θtest .+ 120*π/180 .+ I_ref_ang)

#%% Starting time simulation
println("\nHere we go.\n")

println("Progress : 0.0 %")

for i in 1:Env.N-1

    if i > 1 && floor((10*t[i]/t_final)) != floor((10*t[i - 1]/t_final))
        flush(stdout)
        println("Progress : ", 10*floor((10*t[i]/t_final)), " %")
    end

    # Environment ___________________________________________________________

    Measurements(Env, i)

    # Control System ___________________________________________________________

    #---------------------------------------------------------------------------
    # Source.θpll[num_srce,1,i] - for grid following,
    # Env.θs[i] - for grid forming
    # Source.θ_droop[num_srce, 1, i] - for droop

    num_srce = 1
    V_poc = Env.x[Env.V_poc_loc[: , num_srce], i] # a, b, c components
    I_poc = Env.x[Env.I_poc_loc[: , num_srce], i]
    I_inv = Env.x[Env.I_inv_loc[: , num_srce], i]
    P_Q = Env.p_q_inst[num_srce, :, i]

    Phase_Locked_Loop_3ph(Source, num_srce, i, V_poc)
    Droop_Control(Source, num_srce, i, P_Q, Env)

    #θt = Env.θs[i]
    θt = Source.θ_droop[num_srce, 1, i]
    #θt = Source.θpll[num_srce, 1 , i]

    Voltage_Controller(Kp = 0.01, Ki = 10, Source, num_srce, i, V_poc, θt)
    Current_Controller(Kp = 0.5, Ki = 25, Source, num_srce, i, I_inv, θt)

    #---------------------------------------------------------------------------
    num_srce = 2
    V_poc = Env.x[Env.V_poc_loc[: , num_srce], i]
    I_poc = Env.x[Env.I_poc_loc[: , num_srce], i]
    I_inv = Env.x[Env.I_inv_loc[: , num_srce], i]
    P_Q = Env.p_q_inst[num_srce, :, i]

    Phase_Locked_Loop_3ph(Source, num_srce, i, V_poc, Kp = 1, Ki = 100)
    Droop_Control(Source, num_srce, i, P_Q, Env)

    #θt = Env.θs[i]
    θt = Source.θ_droop[num_srce, 1, i]
    #θt = Source.θpll[num_srce, 1 , i]

    Voltage_Controller(Kp = 0.01, Ki = 10, Source, num_srce, i, V_poc, θt)
    #PQ_Control(Source, num_srce, i, I_poc, V_poc, θt)
    Current_Controller(Kp = 0.5, Ki = 10, Source, num_srce, i, I_inv, θt)

    # System Dynamics __________________________________________________________

    Evolution(Env, Source, i)

end

println("Progress : 100.0 %")

#%% Plots

Plot_I_dq0(0, 20, Source, Env, num_source = 2)

Plot_V_dq0(0, 20, Source, Env, num_source = 2)

Inst_Vout_Vref(30, 50, Source, Env, num_node = 2, num_source = 2)

Inst_Iout_Iref(0, 50, Source, Env, num_source = 2, num_node = 2)

Plot_PLL(0, 20, Source, Env, num_source = 2)

Plot_Droop(0, 20, Source, Env, num_source = 2)

Plot_Irms(0, 10, Env, num_node = 1)

Plot_Vrms(10, 50, Env, Source, num_node = 2, num_source = 2)

#Plot_Real_Imag_Active_Reactive(0, 250, Env, Source, num_node = 1, num_source = 1)
Plot_Real_Imag_Active_Reactive(0, 50, Env, Source, num_node = 2, num_source = 2)

Plot_fft(0, 1, Env, Source, num_node = 2, num_source = 2)

# Calculation Checks
Z1 = conj((3)*Env.V_ph[1, 1,2,Env.N-1]^2/(Env.P[1, 4, Env.N-1] + 1im*Env.Q[1, 4, Env.N-1]))
S1 = conj((3)*Env.V_ph[1, 1,2,Env.N-1]^2/(RL2 + 1im*XL))
I1 = Env.V_ph[1, 1, 2, Env.N-1]/(abs(RL1 + 1im*XL))

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
