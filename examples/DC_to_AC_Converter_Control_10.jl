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

source["pwr"] = 150e3
source["v_rip"] = 0.15
source["vdc"] = 800
source["i_rip"] = 0.01537
source["fltr"] = "LC"
source["R1"] = Source.Rf[1]
source["R_C"] = 0.0006
source["L1"] = Source.Lf[1]
source["C"] = Source.Cf[1];

push!(source_list, source, source);

load_list = []
load = Dict()

R1, L_C, X, Z = Load_Impedance(50e3, 0.6, 230)
R2, L_C, X, Z = Load_Impedance(50e3, -0.99, 230)
load["impedance"] = "RLC"
load["R"] = R1 + R2;
load["L"] = L_C;
load["C"] = L_C;
push!(load_list, load);

cable_list = []

# Network Cable Impedances
l = 0.5 # length in km
cable = Dict()
cable["R"] = 0.222*l # Ω, line resistance
cable["L"] = 0.0024*l # H, line inductance
cable["C"] = 3.7898e-6;
push!(cable_list, cable, cable);

parameters["source"] = source_list
parameters["cable"] = cable_list
parameters["load"] = load_list;
parameters["grid"] = Dict("fs" => f_cntr, "phase" => 3, "v_rms" => 230);

# Define grid using random initialization
power_grid = NodeConstructor(num_sources = 2, num_loads = 1, S2S_p = 1, S2L_p = 1, CM = CM, parameters = parameters);
A, B, C, D = get_sys(power_grid)

ns = length(A[1,:]) # get num of states
ni = length(B[1,:]) # get num of inputs
x0 = [0.0 for i = 1:ns]
env = SimEnv(A = A, B = B, C = C, D = D, x0 = x0, state_ids = get_state_ids(power_grid), rewardfunction = reward)

# find the indices in the state vector that correspond to the inverters
Collect_IDs(env, Source, power_grid.num_fltr_LCL, power_grid.num_fltr_LC, power_grid.num_fltr_L)

c_policy = Classical_Policy(action_space = action_space(env), Source = Source)
s_policy = sin_policy(action_space = action_space(env), ts = Ts)

# Helper logging definitions till history works
state_a = zeros(Source.N_cntr)
state_b = zeros(Source.N_cntr)
state_c = zeros(Source.N_cntr)
c_u_a = zeros(Source.N_cntr)
c_u_b = zeros(Source.N_cntr)
c_u_c = zeros(Source.N_cntr)
s_u_a = zeros(Source.N_cntr)
s_u_b = zeros(Source.N_cntr)
s_u_c = zeros(Source.N_cntr)

#%% Starting time simulation
println("\nHere we go.\n")
reset!(env)
# reset source as well.
# current limitation for PQ mode

@time begin

    println("Progress : 0.0 %")

    for i in 1:Source.N_cntr-1

        # Progress Bar
        if i > 1 && floor((10*t[i]/t_final)) != floor((10*t[i - 1]/t_final))
            flush(stdout)
            #println("Progress : ", 10*floor((10*t[i]/t_final)), " %")
        end

        # Control System _______________________________________________________

        if t[i] > t_final/2
            num_source = 2 # changing the power set points of the 2nd source
            #policy.Source.pq0_set[num_source, 1] = -50e3 # W, Real Power
            #policy.Source.pq0_set[num_source, 2] = 30e3 # VAi, Imaginary Power
        end

        #-----------------------------------------------------------------------

        # System Dynamics ______________________________________________________

        c_action = c_policy(env)
        s_action = s_policy(env)
        env(s_action)
        num_source = 1
        id_a = Source.V_poc_loc[1, num_source]
        id_b = Source.V_poc_loc[2, num_source]
        id_c = Source.V_poc_loc[3, num_source]
        state_a[i] = env.state[id_a]
        state_b[i] = env.state[id_b]
        state_c[i] = env.state[id_c]
        c_u_a[i] = c_action[1 + 3*(num_source - 1)]
        c_u_b[i] = c_action[2 + 3*(num_source - 1)]
        c_u_c[i] = c_action[3 + 3*(num_source - 1)]
        s_u_a[i] = s_action[1 + 3*(num_source - 1)]
        s_u_b[i] = s_action[2 + 3*(num_source - 1)]
        s_u_c[i] = s_action[3 + 3*(num_source - 1)]

    end

    println("Progress : 100.0 %\n")
end


v_out = plot(t, state_b, xlabel="time", ylabel="state")
v_out = plot!(t, s_u_b, xlabel="time")
#v_out = plot(t, Source.V_filt_poc[num_source, 1, :])
#v_out = plot!(t, Source.V_filt_poc[num_source, 2, :])
#v_out = plot!(t, Source.V_filt_poc[num_source, 3, :])
display(v_out)
 
u = plot(t, c_u_a, xlabel="time", ylabel="action", label = "clas")
#u = plot(t, c_u_b, xlabel="time")
#u = plot(t, c_u_c, xlabel="time")
#u = plot!(t, s_u_a, label = "sine")
#u = plot!(t, s_u_b, label = "sine")
#u = plot!(t, s_u_c, label = "sine")
num_source = 1
ph = 1
#u = plot!(t, c_policy.Source.Vd_abc_new[num_source, ph, 1:end-1], label = "inside")
#display(u)

#%% Plots


#Plot_I_dq0(0, 5000, c_policy.Source, num_source = 2)

#Plot_V_dq0(0, 5000, c_policy.Source, num_source = 2)

#Inst_Vout_Vref(10, 20, c_policy.Source, env, num_source = 2)
#=
Inst_Iout_Iref(0, 20, c_policy.Source, env, num_source = 2)

Plot_PLL(0, 500, c_policy.Source, env, num_source = 2, ph = 2)

Plot_Irms(0, 5000, c_policy.Source, num_source = 2)

Plot_Vrms(0, 5000, c_policy.Source, num_source = 1)
Plot_Vrms(0, 5000, polic_policycy.Source, num_source = 2)

Plot_Real_Imag_Active_Reactive(0, 5000, c_policy.Source, num_source = 1)
Plot_Real_Imag_Active_Reactive(0, 5000, c_policy.Source, num_source = 2)
=#

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
