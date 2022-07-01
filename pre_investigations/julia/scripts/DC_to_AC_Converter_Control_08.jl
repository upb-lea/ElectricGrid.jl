using Plots
using LinearAlgebra
using FFTW

include("DC_to_AC_Converter_Controls_01.jl")
include("Classical_Control_Plots.jl")

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

fsys = 50 # the frequency of the reference/control waveform
Vo_rms = 230 # rms output voltage
Vo_p = sqrt(2)*Vo_rms

#_______________________________________________________________________________
# Parameters - Time simulation
Timestep = 1e6/10e3 #time step in μs
t_final = 1 #time in seconds, total simulation run time

# Impedance Calculations
SL = 50e3 #VA, 3-ph Apparent Power
pf = 0.6 #power factor

#_______________________________________________________________________________
# Environment Calcs
Nps = (1/(Timestep*1e-6)) # time intervals
μps = 1/Nps # time step
t = 0:μps:t_final # time

f_cntr = 1/μps # Hz, Sampling frequency of controller ~ 15 kHz -> 50kHz

N = length(t) # number of samples

Source_Control = Source_Controller(t_final, f_cntr, 1)

#_______________________________________________________________________________
# Circuit Elements Calcs
PL = pf*SL #W, Active Power - average instance power
QL = sqrt(SL^2 - PL^2) #VAr, Reactive Power

RL = real(3*Vo_rms^2/(PL + 1im*QL)) # Load/Line Resistor
LL = -1*imag(3*Vo_rms^2/(PL + 1im*QL))/(2*π*fsys) # Load/Line Reactor
XL = (2*π*fsys)*LL
ZL = RL + 1im*XL #Load Impedance

# Reference currents
I_ref = Vo_rms/ZL
I_ref_p = abs(I_ref)*sqrt(2) # peak reference current
I_ref_rms = abs(I_ref) # peak reference current
I_ref_ang = angle(I_ref)

# Inv 1 Source Filter Calcs
Lf, Cf, fc = Filter_Design(Source_Control, 1, 0.15, 0.01537)

#_______________________________________________________________________________
# State space representation
a = [[0 -1/Lf 0];
    [1/Cf 0 -1/Cf];
    [0 1/LL -RL/LL]]
a2 = [[0 -1/Lf 0];
    [1/Cf 0 -1/Cf];
    [0 1/LL -RL*2/LL]]
A = [a zeros(3,3) zeros(3,3);
    zeros(3,3) a zeros(3,3);
    zeros(3,3) zeros(3,3) a]

b = [1/Lf; 0; -1/LL]
B = [b; b; b]

c = [[0 -1 0];
    [1 0 -1];
    [0 1 -RL]]
C = [c zeros(3,3) zeros(3,3);
    zeros(3,3) c zeros(3,3);
    zeros(3,3) zeros(3,3) c]
d = [1.0; 0.0; -1.0]
D = [d; d; d]

num_sources = 1
num_loads = 1
Env = Environment(t_final, μps, A, B, C, D, num_sources, num_loads)

#_______________________________________________________________________________
# Control Reference signals

Source_Control.V_ref[1, 1, :] = Vo_p*sin.(2*π*Env.fs.*t)
Source_Control.V_ref[1, 2, :] = Vo_p*sin.(2*pi*Env.fs.*t .- 120*π/180)
Source_Control.V_ref[1, 3, :] = Vo_p*sin.(2*pi*Env.fs.*t .+ 120*π/180)

Source_Control.I_ref[1, 1, :] = I_ref_p*sin.(2*π*Env.fs.*t .+ I_ref_ang)
Source_Control.I_ref[1, 2, :] = I_ref_p*sin.(2*pi*Env.fs.*t .- 120*π/180 .+ I_ref_ang)
Source_Control.I_ref[1, 3, :] = I_ref_p*sin.(2*pi*Env.fs.*t .+ 120*π/180 .+ I_ref_ang)


#%% Starting time simulation
println("\nHere we go.\n")

println("Progress : 0.0 %")

for i in 1:N-1

    if i > 1 && floor((10*t[i]/t_final)) != floor((10*t[i - 1]/t_final))
        flush(stdout)
        println("Progress : ", 10*floor((10*t[i]/t_final)), " %")
    end

    Measurements(Env, i)

    # System Dynamics __________________________________________________________

    Discrete_time(Env, i)

    # Control System ___________________________________________________________

    V_poc = [Env.x[2, i]; Env.x[5, i]; Env.x[8, i]] # a, b, c components
    I_inv = [Env.x[1, i]; Env.x[4, i]; Env.x[7, i]]

    Phase_Locked_Loop_3ph(Source_Control, 1, i, V_poc)

    #Source_Control.θpll[1,1,i] - for grid following

    Voltage_Controller(Source_Control, 1, i, V_poc, Env.θs[i])
    Current_Controller(Source_Control, 1, i, I_inv, Env.θs[i])

    # Inverter Voltages - Control Actions
    dly = 1
    Env.u[1, i + dly] = Source_Control.Vd_abc_new[1, 1, i]
    Env.u[4, i + dly] = Source_Control.Vd_abc_new[1, 2, i]
    Env.u[7, i + dly] = Source_Control.Vd_abc_new[1, 3, i]

    # 2nd Source
    on = 0
    Env.u[3, i + dly] = on*Vo_p*sin(Env.θs[i])
    Env.u[6, i + dly] = on*Vo_p*sin(Env.θs[i] - 120*π/180)
    Env.u[9, i + dly] = on*Vo_p*sin(Env.θs[i] + 120*π/180)

end

println("Progress : 100.0 %")

#%% Plots

#Plot_I_dq0(1, 0, 10, Source_Control)

#Plot_V_dq0(1, 0, 10, Source_Control)

#Inst_Vout_Vref(1, 1, 0, 10, Source_Control, Env)

#Inst_Iout_Vref(1, 1, 0, 20, Source_Control, Env)

Plot_PLL(1, 0, 20, Source_Control, Env)

#Plot_Irms(1, 0, 10, Env)

Plot_Vrms(1, 1, 0, 10, Env, Source_Control)

Plot_Real_Imag_Active_Reactive(1, 1, 0, 10, Env, Source_Control)

#Plot_fft(1, 1, 0, 1, Env, Source_Control)

# Calculation Checks
Z1 = conj((3)*Env.V_ph[1, 1,2,N-1]^2/(Env.P[1, 4, N-1] + 1im*Env.Q[1, 4, N-1]))
S1 = conj((3)*Env.V_ph[1, 1,2,N-1]^2/(RL + 1im*XL))
I1 = Env.V_ph[1, 1, 2, N-1]/(abs(RL + 1im*XL))

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
