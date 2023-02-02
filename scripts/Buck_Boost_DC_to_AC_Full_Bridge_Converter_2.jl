using SpecialFunctions
using LinearAlgebra
using Plots
using StatsBase
using LinearAlgebra
using FFTW
using Expokit #for exponentials of sparse matrices
#using GraphMakie
#using GLMakie, Makie
#using LightGraphs

function Simpsons_Rule(f, start, finish, n, A, B, u)
    # numerical integration technique
    # n should be a positive even integer
    # f is another function that returns the value of an equation

    S = zeros(size(A, 1))
    Δ = (finish - start)/n

    for i = 0:n

        if i == 0 || i == n

            S = S + f(i*Δ + start, A, B, u)
            #println("f(",i*Δ + start,") = ", f(i*Δ + start))
        elseif mod(i,2) == 1

            S = S + 4*f(i*Δ + start, A, B, u)
            #println("f(",i*Δ + start,") = ", 4*f(i*Δ + start))
        elseif mod(i,2) == 0

            S = S + 2*f(i*Δ + start, A, B, u)
            #println("f(",i*Δ + start,") = ", 2*f(i*Δ + start))
        end
    end

    S = Δ*S/3

    return S
end

function Particular_Solution(A, B, u, start, finish)

    local xp

    xp = Simpsons_Rule(xp_f, start, finish, 4, A, B, u)
    xp = exp(finish*A)*xp

    return xp
end

function xp_f(t, A, B, u)
    f = expmv(-t, A, B*u)
    return f
end

function Continuous_time(A, B, u, t, μ, x0)

    xp = Array{Float64, 1}(undef, size(A,1))

    if u == 0 || isempty(B)
        xp = fill!(xp, 0.0)
    else
        xp = Particular_Solution(A, B, u, t, t + μ)
    end

    x = exp(μ*A)*x0 + xp

    return x
end

function Discrete_time(Ad, Bd, u, k, x0)

    xp = Array{Float64, 1}(undef, size(Ad,1))
    xp = fill!(xp, 0.0)

    if isempty(Bd) || u == 0
        xp = xp
    else
        for j in 0:(k-1)
            xp = xp + (Ad^((k-1) - j))*u*Bd
        end
    end

    x = (Ad^k)*x0 + xp

    # Brute force - not discrete time matrices
    #dx_dt = A*x0+ u*B
    #x = x0 + μ.*dx_dt

    return x
end

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

#%% Parameters - Buck Boost stage
#______________________________________________________________________________
Vd = 50.0 #input DC voltage
VCB = 400.0 #the reference DC voltage of the buck-boost stage

RB = 0.025 # parasitic resistance of inductor
LB = 625e-6 # Buck Boost Inductor

ΔVCB_VCB = 0.01
# The ratio of the ripple voltage to the capacitor voltage that we use to design CB

fB = 10000.0 # the frequency of the saw tooth waveform

#%% Parameters - DC to AC Full Bridge stage
#______________________________________________________________________________
Vo = 230 # rms output voltage
fr = 50 # the frequency of the reference/control waveform

RL = 3.1738 # the load resistor

ΔILf_ILf = 0.05
# specifying the maximum current ripple current through the filter inductor,
# allows us to design the filter inductor, Lf

ΔVCf_VCf = 0.025
# specifying the maximum current ripple current through the filter capacitor,
# allows us to design the filter capacitor, Cf

fs = 5000 # the frequency of the triangle waveform

tdΔ = 5e-6
#=
In practise the converter switching pairs do not switch directly on and off
immediately after each other. To avoid both pairs of switches being on at the
same time the pulse-width modulator 'waits' for a time, tdΔ, (the dead time)
from when the one pair has switched off until the following switch is switched
on. For IGBT converters the dead-time is typically 5e-6 and for MOSFET converters
it is typically 1e-6.
Dead-time causes the following two problems:
1. Deforming of the output voltage waveform
2. It limits the maximum value of the modulation index
=#

#%% Parameters - Time simulation
#______________________________________________________________________________

Timestep = 1 #time step in μs (should be less that dead time)
t_final = 0.3 #time is seconds

println("\n_____________________________________________________________________")
println("Simulation Parameters\n")
println("Total simulation time: ", round(t_final, digits = 2), " s")
println("Simulation time step: ", Timestep, " μs")

#%% Buck Boost theoretical design calculations
#______________________________________________________________________________
TB = 1/fB # period of the saw tooth waveform
vsaw_p = 1 # peak of saw tooth signal
mB = vsaw_p/TB # the gradient of the saw waveform

IBo_Id = Vd/VCB # The approximate ratio of the output current to the input voltage
RL_dc = RL
IBo = VCB/RL # The output current of the buck boost stage,
# IBo = ILf_rms

DB = (VCB)/(VCB + Vd) # Duty cycle,
# If the parasitic resistance, RB, of the inductor, LB, is neglected,
# the output voltage, VCB, tends to infinity as D-> 1
cntr_B = DB*vsaw_p # control signal for buck boost

ΔVCB = ΔVCB_VCB * VCB # Ripple in output voltage, the voltage over the DC capacitor
ΔILB = DB*TB*Vd/LB # Ripple in inductor current
CB = (DB*TB*VCB)/(RL_dc*ΔVCB)
τ = RL_dc*CB # the time constant of the filter capacitor and load resistor.

IB_dis = (VCB*(1 - DB)^2 * TB)/(2*LB)
# this value of the output current from the buck boost stage defines
# the boundary between continuous and discontinuous current
RB_dis = VCB/IB_dis

#=
While the transient response of the buck boost converter may be discontinuous,
if the load decreases to such an extent that the effective DC resistance
is greater than RB_dis, then the discontinous effects will be permanent
=#

# Other Steady state / average theoretical values
# The average inductor current is IL = Id + Io, Vd*Id = Vo*Io

Id = (DB/(1 - DB))*IBo
ILB = Id + IBo # = IBo(1 - DB)
VLB = Vd*DB - VCB*(1 - DB)

#= Effect of the parasitic resistance, RB
When the duty-cycle value approaches 1 the inductor current become very large.
The average inductor current is equal to the average supply current plus the
average load current.
It causes large losses in the inductor.
The result is then that the value of the output voltage is lower than the
theoretic prediction.
=#

println("\n_____________________________________________________________________")
println("Buck Boost - Design Parameters\n")

println("Input DC voltage, Vd: ", Vd, " V")
println("Output voltage: ", VCB, " V")
println("Duty cycle, D: ", round(DB, digits = 2))
println("Voltage ripple, Δvc: ",
round(ΔVCB, digits = 2), " V")
println("Inductor current ripple, ΔIL: ",
round(ΔILB, digits = 2), " A")
println("Average inductor current*, ILB: ",
round(ILB, digits = 2), " A") #ma
println("Ripple voltage to Output voltage, ΔVc/Vc: ",
ΔVCB_VCB*100, " %")
println("Discontinuous current boundary, R: ",
round(RB_dis, digits = 4), " Ω")
println("Filter capacitor: ",
round(CB*1e6, digits = 2), " μF")
println("Time constant of filter capacitor, τ = RC: ",
round(τ, digits = 4), " s") #typically the transient period is taken to be 4*τ

#%% Full Bridge Converter theoretical calculations
#______________________________________________________________________________

Ts = 1/fs # period of the triangle waveform
# During one switching period, each switch switches on once and off once. Thus
# fs is also called the switching frequency
fs_p = 1 # the amplitude of the triangle wave
ms = 4*fs_p/Ts # the gradient of the triangle waveform

mf = fs/fr # the frequency modulation index

vo_p = Vo*sqrt(2) # the amplitude of the output voltage
io_p = vo_p/RL # the amplitude of the load current

ma = vo_p/VCB # fr_p/fs_p the modulation index
#=
We can implement over-modulation, increasing the modulation index, to above 1
The amplitude of the sinusoidal reference is thus larger than that of the
triangle waveform.
The advantage is that the amplitude of the fundamental component is now larger
than VcB
The disadvantage is that the quality of the output decreases.
In the case of overmodulation, the switching frequeny becomes lower because the
converter does not switch for segments of a cycle.
When ma is very large the converter switches at fr (e.g. 50Hz). This is known
as square-wave switching. In this case the amplitude of the fundamental is
given by Vcf = (4*VcB)/π
=#

fr_p = ma*fs_p # the amplitude of the reference signal

ΔILf_max = ΔILf_ILf*io_p
# we can assume that vo is constant over a switching period and positive
# we find that the maximum current ripple is when the duty cycle of a pair of
# switches is at 50%
Lf = VCB/(2*fs*ΔILf_max)

ΔVCf_max = ΔVCf_VCf*Vo
Cf = ΔILf_max/(8*ΔVCf_max*fs)

fc = 1/(2*π*sqrt(Lf*Cf))
# the cut-off frequency of the filter circuit

println("\n_____________________________________________________________________")
println("AC-DC Converter - Design Parameters\n")

println("RMS output voltage, Vo: ", Vo, " V")
println("Frequency of the reference/control waveform, fr: ", fr, " Hz")
println("Switching frequency of the triangle waveform, fs: ", fs, " Hz")
println("Output load resistor: ",
RL, " Ω")
println("Frequency modulation index, mf: ", round(mf, digits = 2))
println("Modulation index, ma: ", round(ma, digits = 2))
println("Filter inductor: ",
round(Lf*1e3, digits = 2), " mH")
println("Filter capacitor: ",
round(Cf*1e6, digits = 2), " μF")
println("Filter cut-off frequency: ",
round(fc, digits = 2), " Hz")
println("Maximum Voltage ripple, Δvo: ",
round(ΔVCf_max, digits = 2), " V")
println("Maximum Filter Inductor current ripple, ΔIf: ",
round(ΔILf_max, digits = 2), " A")
println("Maximum Ripple voltage to Output voltage, ΔVo/Vo: ",
ΔVCf_VCf*100, " %")
println("Maximum Ripple current to Inductor current, ΔILF/ILF: ",
ΔILf_ILf*100, " %")

#%% initialisation
#______________________________________________________________________________

# simulation
Nps = (1/(Timestep*1e-6)) - 1 # time intervals, samples in a second + correction
μps = 1/Nps #corrected time step
t_final = (t_final - μps)
t = 0:μps:t_final # time

N = length(t) # number of samples

# buck boost control
Saw = Array{Float64, 1}(undef, N)
ControlB = Array{Int64, 1}(undef, N)

Saw = fill!(Saw, 0.0)
ControlB = fill!(ControlB, 0)

# DC_AC control
Rf = Array{Float64, 1}(undef, N)
Tri = Array{Float64, 1}(undef, N)
ControlA = Array{Int64, 1}(undef, N)
ControlA = fill!(ControlA, 0)
toggl = 1
dead_count = 0

# state space variables
iLB = Array{Float64, 1}(undef, N)
vCB = Array{Float64, 1}(undef, N)
iLf = Array{Float64, 1}(undef, N)
vCf = Array{Float64, 1}(undef, N)

iLB = fill!(iLB, 0)
vCB = fill!(vCB, 0)
iLf = fill!(iLf, 0)
vCf = fill!(vCf, 0)

# output variables
vLB = Array{Float64, 1}(undef, N)
iCB = Array{Float64, 1}(undef, N)
vLf = Array{Float64, 1}(undef, N)
iCf = Array{Float64, 1}(undef, N)

# initial conditions

iLB[1] = 0
vCB[1] = 0
iLf[1] = 0
vCf[1] = 0

println("\n_____________________________________________________________________")
println("Initial conditions\n")
println("Buck Boost Inductor Current, iLB[t = 0]: ", iLB[1], " A")
println("Buck Boost Capacitor Voltage, vCB[t = 0]: ", vCB[1], " V")
println("DC-AC Converter Inductor Current, iLf[t = 0]: ", iLf[1], " A")
println("DC-AC Converter Capacitor Voltage, vCf[t = 0]: ", vCf[1], " V")

# state space representations
#=
dx/dt = A*x + B
K = exp(-t0*A)*(x0 - xp)
x(t) = exp(t*A)*K + xp

Vector K is calculated from initial conditions
The particular solution can be determined through inspection by calculating
the values of the state variables in the steady-state. It is given by the
inductor current and capacitor voltage after the circuit has been in a specific
switching state for a long time.

Which can also be written as, if the system is a linear homogeneous ODE:
Ψ(t) = exp(t*A)
Ψ^-1(t0) = exp(-t0*A)
x(t, t0, x0) = Ψ(t)*Ψ^-1(t0)*x0 = Φ(t,t0)*x0 ∀t,
where Φ(t,t0) is sometimes called the state transition matrix.
The state transition matrix is uniquely determined by A and independent of
the particular choice of Ψ(t).

Furthermore, Φ(t0,t0) = I, the identity matrix
dΦ(t,t0)/dt = A*Φ(t,t0)
Φ(t,τ) = Φ(t,σ)*Φ(σ,β)*Φ(β,τ)
Φ(t,t0)^-1 = Φ(t0,t), i.e. going backwards in time will give the same initial state

The particular solution can also be found through:
xp = ∫Φ(t,τ)*B(τ) dτ, integrating from t0 to t
the particular solution is due to the forcing term B(t)

For time-invariant A, a fundamental matrix is the matrix exponential
Ψ(t) = exp(A*t)
The state transition matrix is
Φ(t,t0) = exp(A(t - t0))

In power electronic power supplies we constantly switch between two different
linear circuits
The one circuit is represented by state matrices A1, B1, and C1; and the other
circuit by state matrices A2, B2, and C2.
The only thing "connecting" the two circuits is the fact that end values of the
state variables of the one circuit becomes the starting values of the following
circuit. The end values of the state variables for the one switching state
becomes the starting values for the following state.
=#

# initialise
# Td is closed
A01 = [[(-RB/LB) 0 0 0]; # TA+ and TB- closed
        [0 0 1/CB 0];
        [0 -1/Lf 0 -1/Lf];
        [0 0 1/Cf -1/(RL*Cf)]]
A02 = [[(-RB/LB) 0 0 0]; # TA- and TB+ closed
        [0 0 -1/CB 0];
        [0 1/Lf 0 -1/Lf];
        [0 0 1/Cf -1/(RL*Cf)]]

# Td is open
A11 = [[(-RB/LB) -1/LB 0 0]; # TA+ and TB- closed
        [1/CB 0 1/CB 0];
        [0 -1/Lf 0 -1/Lf];
        [0 0 1/Cf -1/(RL*Cf)]]
A12 = [[(-RB/LB) -1/LB 0 0]; # TA- and TB+ closed
        [1/CB 0 -1/CB 0];
        [0 1/Lf 0 -1/Lf];
        [0 0 1/Cf -1/(RL*Cf)]]

#first order liner non-homogeneous equation
B01 = [1/LB; 0; 0; 0]
B02 = [1/LB; 0; 0; 0]

#first order liner homogeneous equation - no particular solution
B11 = [0; 0; 0; 0]
B12 = [0; 0; 0; 0]

C01 = [[-RB 0 0 0];
        [0 0 1 0];
        [0 -1 0 -1];
        [0 0 1 -1/RL]]
C02 = [[-RB 0 0 0];
        [0 0 -1 0];
        [0 1 0 -1];
        [0 0 1 -1/RL]]

C11 = [[-RB -1 0 0];
        [1 0 1 0];
        [0 -1 0 -1];
        [0 0 1 -1/RL]]
C12 = [[-RB -1 0 0];
        [1 0 -1 0];
        [0 1 0 -1];
        [0 0 1 -1/RL]]

D01 = [1; 0; 0; 0]
D02 = [1; 0; 0; 0]

D11 = [0; 0; 0; 0]
D12 = [0; 0; 0; 0]

#%% Discrete time matrices
#______________________________________________________________________________

Ad01 = exp(A01*μps)
Ad02 = exp(A02*μps)
Bd01 = inv(A01)*(Ad01 - Matrix(I, size(A01,1), size(A01,1)))*B01
Bd02 = inv(A02)*(Ad02 - Matrix(I, size(A02,1), size(A02,1)))*B02

Ad11 = exp(A11*μps)
Ad12 = exp(A12*μps)
Bd11 = inv(A11)*(Ad11 - Matrix(I, size(A11,1), size(A11,1)))*B11
Bd12 = inv(A12)*(Ad12 - Matrix(I, size(A12,1), size(A12,1)))*B12

#%% Starting time simulation
#______________________________________________________________________________

for i in 1:N

    # Buck Boost Control
    #___________________________________________________________________________
    if i < N
        if ms/Nps + Saw[i] >= vsaw_p && Saw[i] != vsaw_p
            Saw[i+1] = vsaw_p
        elseif Saw[i] >= vsaw_p
            Saw[i+1] = mB/Nps
        else
            Saw[i+1] = mB/Nps + Saw[i]
        end
    end

    if cntr_B >= Saw[i]
        ControlB[i] = 1
    end

    # DC to AC Control
    #___________________________________________________________________________
    global toggl, dead_count
    Rf[i] = fr_p*sin(2*pi*fr*t[i])

    if Tri[i] >= fs_p && toggl == 1
        toggl = -1
    elseif Tri[i] <= -fs_p && toggl == -1
        toggl = 1
    end

    if i < N
        Tri[i+1] = toggl*ms/Nps + Tri[i]
    end

    if Rf[i] >= Tri[i] && i > 1

        if ControlA[i-1] == -1
            dead_count = t[i]
        end

        if t[i] - dead_count < tdΔ
            ControlA[i] = 0
        else
            ControlA[i] = 1
        end
    elseif Rf[i] < Tri[i] && i > 1

        if ControlA[i-1] == 1
            dead_count = t[i]
        end

        if t[i] - dead_count < tdΔ
            ControlA[i] = -0
        else
            ControlA[i] = -1
        end
    end

    # by setting this time and ControlA value we set when the AC conversion starts
    if t[i] < 0
        ControlA[i] = -1 # 0. 1, or -1
    end
    # System dynamics
    #___________________________________________________________________________
    if i < N

        if ControlB[i] == 1 && i < N # if Td is closed

            if ControlA[i] == 1 # 1 # TA+ and TB- are closed, TA- and TB+ are open

                x0 = [iLB[i]; vCB[i]; iLf[i]; vCf[i]]

                #x = Continuous_time(A01, B01, Vd, t[i], μps, x0)
                x = Discrete_time(Ad01, Bd01, Vd, 1, x0)

                # Update state variable
                iLB[i + 1] = x[1];
                vCB[i + 1] = x[2];
                iLf[i + 1] = x[3];
                vCf[i + 1] = x[4];

                #  Outputs
                y = C01*x + D01*Vd

                # Update output variable
                vLB[i + 1] = y[1];
                iCB[i + 1] = y[2];
                vLf[i + 1] = y[3];
                iCf[i + 1] = y[4];

            elseif ControlA[i] == -1 # -1 # TA- and TB+ are closed, TA+ and TB- are open

                x0 = [iLB[i]; vCB[i]; iLf[i]; vCf[i]]

                # Exact Method
                #x = Continuous_time(A02, B02, Vd, t[i], μps, x0)
                x = Discrete_time(Ad02, Bd02, Vd, 1, x0)

                #  Outputs
                y = C02*x + D02*Vd

                # Update state variables
                iLB[i + 1] = x[1];
                vCB[i + 1] = x[2];
                iLf[i + 1] = x[3];
                vCf[i + 1] = x[4];

                # Update output variables
                vLB[i + 1] = y[1];
                iCB[i + 1] = y[2];
                vLf[i + 1] = y[3];
                iCf[i + 1] = y[4];

            elseif ControlA[i] == 0 # 0 # both pairs of switches are open

                x0 = [iLB[i]; vCB[i]; iLf[i]; vCf[i]]

                if x0[3] >= 0

                    #x = Continuous_time(A01, B01, Vd, t[i], μps, x0)
                    x = Discrete_time(Ad01, Bd01, Vd, 1, x0)

                    #  Outputs
                    y = C01*x + D01*Vd
                else

                    #x = Continuous_time(A02, B02, Vd, t[i], μps, x0)
                    x = Discrete_time(Ad02, Bd02, Vd, 1, x0)

                    #  Outputs
                    y = C02*x + D02*Vd
                end

                # Update state variables
                iLB[i + 1] = x[1];
                vCB[i + 1] = x[2];
                iLf[i + 1] = x[3];
                vCf[i + 1] = x[4];

                # Update output variable
                vLB[i + 1] = y[1];
                iCB[i + 1] = y[2];
                vLf[i + 1] = y[3];
                iCf[i + 1] = y[4];
            end

        elseif ControlB[i] == 0 && i < N # if Td is open

            if ControlA[i] == 1 # 1 # TA+ and TB- are closed, TA- and TB+ are open

                x0 = [iLB[i]; vCB[i]; iLf[i]; vCf[i]]

                # Exact Method
                #x = Continuous_time(A11, [], 0.0, t[i], μps, x0)
                x = Discrete_time(Ad11, [], 0.0, 1, x0)

                # Update state variable
                iLB[i + 1] = x[1];
                vCB[i + 1] = x[2];
                iLf[i + 1] = x[3];
                vCf[i + 1] = x[4];

                #  Outputs
                y = C11*x + D11*Vd

                # Update output variable
                vLB[i + 1] = y[1];
                iCB[i + 1] = y[2];
                vLf[i + 1] = y[3];
                iCf[i + 1] = y[4];

            elseif ControlA[i] == -1 # -1 # TA- and TB+ are closed, TA+ and TB- are open

                x0 = [iLB[i]; vCB[i]; iLf[i]; vCf[i]]

                # Exact Method
                #x = Continuous_time(A12, [], 0, t[i], μps, x0)
                x = Discrete_time(Ad12, [], 0.0, 1, x0)

                # Update state variable
                iLB[i + 1] = x[1];
                vCB[i + 1] = x[2];
                iLf[i + 1] = x[3];
                vCf[i + 1] = x[4];

                #  Outputs
                y = C12*x + D12*Vd

                # Update output variable
                vLB[i + 1] = y[1];
                iCB[i + 1] = y[2];
                vLf[i + 1] = y[3];
                iCf[i + 1] = y[4];

            elseif ControlA[i] == 0 # 0 # both pairs of switches are open

                x0 = [iLB[i]; vCB[i]; iLf[i]; vCf[i]]

                if x0[3] >= 0

                    # Exact Method
                    #x = Continuous_time(A11, [], 0, t[i], μps, x0)
                    x = Discrete_time(Ad11, [], 0.0, 1, x0)

                    #  Outputs
                    y = C11*x + D11*Vd
                else

                    # Exact Method
                    x = Continuous_time(A12, [], 0, t[i], μps, x0)
                    x = Discrete_time(Ad12, [], 0.0, 1, x0)

                    #  Outputs
                    y = C12*x + D12*Vd
                end

                # Update state variables
                iLB[i + 1] = x[1];
                vCB[i + 1] = x[2];
                iLf[i + 1] = x[3];
                vCf[i + 1] = x[4];

                # Update output variable
                vLB[i + 1] = y[1];
                iCB[i + 1] = y[2];
                vLf[i + 1] = y[3];
                iCf[i + 1] = y[4];
            end
        end

        if iLB[i + 1] < 0 # boundary between continuous and discontinuous current
            iLB[i + 1] = 0
        end
    end
end

# Harmonic Analysis
#______________________________________________________________________________
T_B = t_final*fB # periods for the buck boost stage, i.e. simulation time
T_r = t_final*fr # periods for the reference stage, i.e. simulation time

per = 5
T_plot_start = T_r-per
N_plot_start = convert(Int64, round((T_plot_start/fr)*Nps))
N_dif = N-N_plot_start

shifted_k = fftshift(fftfreq(length(N_plot_start:N), Nps))
VCf_shifted_fft = 2/N_dif*fftshift(fft(vCf[N_plot_start:N]))
ILf_shifted_fft = 2/N_dif*fftshift(fft(iLf[N_plot_start:N]))

k_start = length(VCf_shifted_fft)/2 + 1
k_start = convert(Int64, round(k_start))
k_end = k_start + length(VCf_shifted_fft)/400 # the /xxxx scales the x axis
k_end = convert(Int64, round(k_end))

freq = shifted_k[k_start:k_end]
VCf_harm_mag = abs.(VCf_shifted_fft[k_start:k_end])
VCf_harm_ang = angle.(VCf_shifted_fft[k_start:k_end])
ILf_harm_mag = abs.(ILf_shifted_fft[k_start:k_end])
ILf_harm_ang = angle.(ILf_shifted_fft[k_start:k_end])

vCf_m1 = VCf_harm_mag[per] #sometimes should be per + 1
vCf_a1 = VCf_harm_ang[per]
ILf_m1 = ILf_harm_mag[per]
ILf_a1 = ILf_harm_ang[per]

VCf_signal = vCf_m1.*cos.(2π * freq[per] .* t .+ (vCf_a1))
ILf_signal = ILf_m1.*cos.(2π * freq[per] .* t .+ (ILf_a1))

THD_vCf = 100*std(VCf_harm_mag[(per+2):end])*sqrt(length(VCf_harm_mag[(per+2):end]))/vCf_m1
# Design accuracy evaluation
#______________________________________________________________________________

T_plot_eval = T_B - 2
N_plot_eval = convert(Int64, round((T_plot_eval/fB  - 1/Nps)*Nps))

ΔvCB = maximum(vCB[N_plot_eval : N]) - minimum(vCB[N_plot_eval : N])
ΔiLB = maximum(iLB[N_plot_eval : N]) - minimum(iLB[N_plot_eval : N])
VCB_average = sum(vCB[N_plot_eval : N])/length(vCB[N_plot_eval : N])
iLB_average = sum(iLB[N_plot_eval : N])/length(iLB[N_plot_eval : N])

print("")
println("\n_____________________________________________________________________")
println("Buck Boost - Design Verification")

println("\nAverage Output DC voltage, VCB: ",
round(VCB_average, digits = 2), " V")
println("Output voltage ripple, Δvc: ",
round(ΔvCB, digits = 2), " V")
println("Inductor current ripple, ΔIL: ",
round(ΔiLB, digits = 2), " A")
println("Average inductor current, ILB: ",
round(iLB_average, digits = 2), " A")
println("Ripple voltage to Output voltage, ΔVc/Vc: ",
round(100*ΔvCB/VCB_average, digits = 2), " %")

T_plot_eval = T_r - 3
N_plot_eval = convert(Int64, round((T_plot_eval/fr  - 1/Nps)*Nps))

vo_average = sum(vCf[N_plot_eval : N])/length(vCf[N_plot_eval : N])
vo_peak = maximum(vCf[N_plot_eval : N])
vo_rms = std(vCf[N_plot_eval : N])

# first need to find where the ripple is the biggest
# The maximum value of the ripples occur when the duty cycle of switches are
# at 50%, d1 = 0.5. That is the ripple is largest when the reference signal
# has a zero crossing.

T_plot_eval_a = T_r - 0.5

T_s = T_plot_eval_a*fs/fr
T_plot_eval_bi = T_s - 0.55
T_plot_eval_bii = T_s + 0.55

N_plot_eval_bi = convert(Int64, round((T_plot_eval_bi/fs  - 1/Nps)*Nps))
N_plot_eval_bii = convert(Int64, round((T_plot_eval_bii/fs  - 1/Nps)*Nps))

d = (vCf[N_plot_eval_bi : N_plot_eval_bii] - VCf_signal[N_plot_eval_bi : N_plot_eval_bii])
e = (iLf[N_plot_eval_bi : N_plot_eval_bii] - ILf_signal[N_plot_eval_bi : N_plot_eval_bii])

Δvo = maximum(d) - minimum(d)
ΔiLf = maximum(e) - minimum(e)

println("\n_____________________________________________________________________")
println("DC to AC Conversion - Design Verification")

println("\nRMS Output Voltage, VCf: ",
round(vo_rms, digits = 2), " V")
println("Maximum Voltage ripple, Δvo: ",
round(Δvo, digits = 2), " V")
println("Maximum Filter Inductor current ripple, ΔIf: ",
round(ΔiLf, digits = 2), " A")
println("Maximum Ripple voltage to Output voltage*, ΔVo/Vo: ",
round(100*Δvo/vo_rms, digits = 2), " %")
println("Maximum Ripple current to Inductor current, ΔILF/ILF: ",
round(100*ΔiLf*RL/vo_peak, digits = 2), " %")
println("Voltage THD: ",
round(THD_vCf, digits = 2), " %")

# Plots
#______________________________________________________________________________

FFTs_p = plot(freq, VCf_harm_mag, seriestype = :scatter, line=:stem,
title = "Frequency Spectrum",
label = "Shifted FFT",
markershape = :circle, markersize = 2,
#xticks=shifted_k[1:10:100],
legend=true);

display(FFTs_p)

p2 = plot(t[N_plot_eval_bi : N_plot_eval_bii], e,
title = "Converter Ripple Current",
label = "iLf",
xlabel = "Time [s]",
ylabel = "Current [A]")

display(p2)

p3 = plot(t[N_plot_eval_bi : N_plot_eval_bii], d,
title = "Converter Ripple Voltage",
label = "VLf",
xlabel = "Time [s]",
ylabel = "Electrical Potential [V]")

display(p3)

T_plot1 = 3
N_plot1 = convert(Int64, round((T_plot1/fB  - 1/Nps)*Nps))
p1a = plot(t[1:N_plot1], Saw[1:N_plot1], label = "Saw",
xlabel = "Time [s]",
title = "Buck Boost Control")
p1a = plot!(t[1:N_plot1], cntr_B*ones(N_plot1), label = "vc")
p1a = plot!(t[1:N_plot1], ControlB[1:N_plot1], label = "Control")

T_plot2 = 3
N_plot2 = convert(Int64, round((T_plot2/fs  - 1/Nps)*Nps))
p1b = plot(t[1:N_plot2], Rf[1:N_plot2], label = "Reference",
xlabel = "Time [s]",
title = "DC-AC Converter Control")
p1b = plot!(t[1:N_plot2], Tri[1:N_plot2], label = "Triangle")
p1b = plot!(t[1:N_plot2], ControlA[1:N_plot2], label = "Control")

PWM = plot(p1a, p1b, layout = (2, 1), legend = true)

display(PWM)

T_plot_start = 0
T_plot_end = T_B
N_plot_start = convert(Int64, round((T_plot_start/fB  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fB  - 1/Nps)*Nps))

p2 = plot(t[N_plot_start : N_plot_end], iLB[N_plot_start : N_plot_end],
title = "Buck Boost Inductor Current",
label = false,
xlabel = "Time [s]",
ylabel = "Current [A]")

display(p2)

p3 = plot(t[N_plot_start : N_plot_end], vLB[N_plot_start : N_plot_end],
title = "Buck Boost Inductor Voltage",
xlabel = "Time [s]",
ylabel = "Electrical Potential [V]",
label = false)

#display(p3)

p4 = plot(t[N_plot_start : N_plot_end], vCB[N_plot_start : N_plot_end],
title = "Buck Boost Output Voltage",
label = false,
xlabel = "Time [s]",
ylabel = "Electrical Potential [V]")

display(p4)

T_plot_start = T_r-2
T_plot_end = T_r
N_plot_start = convert(Int64, round((T_plot_start/fr  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fr  - 1/Nps)*Nps))

p5 = plot(t[N_plot_start : N_plot_end], iLf[N_plot_start : N_plot_end],
title = "DC-AC Converter Inductor current",
label = "output",
xlabel = "Time [s]",
ylabel = "Current [A]")
p5 = plot!(t[N_plot_start:N], ILf_signal[N_plot_start:N],
label = "fundamental")

display(p5)

p6 = plot(t[N_plot_start : N_plot_end], vCf[N_plot_start : N_plot_end],
title = "DC-AC Converter Output Voltage",
label = "output",
xlabel = "Time [s]",
ylabel = "Electrical Potential [V]")
p6 = plot!(t[N_plot_start : N_plot_end], -vo_p.*Rf[N_plot_start : N_plot_end]./fr_p,
label = "reference")
p6 = plot!(t[N_plot_start:N], VCf_signal[N_plot_start:N],
label = "fundamental")

display(p6)

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
