using SpecialFunctions
using LinearAlgebra
using Plots
using StatsBase
using LinearAlgebra
using FFTW
using Expokit #for exponentials of sparse matrices

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
Vd = 200.0 #input DC voltage
VCB = 800.0 #the reference DC voltage of the buck-boost stage

RB = 0.025 #0.025 parasitic resistance of inductor
LB = 625e-6 # Buck Boost Inductor

ΔVCB_VCB = 0.01 #0.01
# The ratio of the ripple voltage to the capacitor voltage that we use to design CB

fB = 20000.0 # the frequency of the saw tooth waveform

#%% Parameters - DC to AC Full Bridge stage
#______________________________________________________________________________
Vo = 230 # rms output voltage
fr = 50 # the frequency of the reference/control waveform
fsys = fr
tr = 1/fr # period length

# Impedance Calculations
SL = 50e3 #VA, 3-ph Apparent Power
pf = 0.9999 #power factor
PL = pf*SL #W, Active Power - average instance power
QL = sqrt(SL^2 - PL^2) #VAr, Reactive Power

RL = real(3*Vo^2/(PL + 1im*QL)) # Load Resistor
Ll = -1*imag(3*Vo^2/(PL + 1im*QL))/(2*π*fsys) # Load Reactor
XL = (2*π*fsys)*Ll
ZL = RL + 1im*XL #Load Impedance

#RL = 5 # the load resistor
#Ll = 1e-6 # the load inductor
Rf = 0.05 #0.05 parasitc resistance of filter inductor

ΔILf_ILf = 0.15
# specifying the maximum current ripple current through the filter inductor,
# allows us to design the filter inductor, Lf

ΔVCf_VCf = 0.01537
# specifying the maximum current ripple current through the filter capacitor,
# allows us to design the filter capacitor, Cf

fs = 15000 # the frequency of the triangle waveform

tdΔ = -5e-6
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

Timestep = 0.01 #time step in μs (should be less than dead time)
t_final = 0.2 #time in seconds

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
IBo = 3*VCB/abs(ZL) # The output current of the buck boost stage,
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

#%% Full Bridge Converter theoretical calculations
#______________________________________________________________________________

Ts = 1/fs # period of the triangle waveform
# During one switching period, each switch switches on once and off once. Thus
# fs is also called the switching frequency
fs_p = 1 # the amplitude of the triangle wave
ms = 4*fs_p/Ts # the gradient of the triangle waveform

mf = fs/fr # the frequency modulation index

vo_p = Vo*sqrt(2) # the amplitude of the output voltage
io_p = vo_p/abs(ZL) # the amplitude of the load current

ma = 2*vo_p/VCB # fr_p/fs_p the modulation index
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

fr_p = 1.104*ma*fs_p # the amplitude of the reference signal

ΔILf_max = ΔILf_ILf*io_p
# we can assume that vo is constant over a switching period and positive
# we find that the maximum current ripple is when the duty cycle of a pair of
# switches is at 50%
Lf = VCB/(4*fs*ΔILf_max)

ΔVCf_max = ΔVCf_VCf*vo_p
Cf = ΔILf_max/(8*ΔVCf_max*fs)

fc = 1/(2*π*sqrt(Lf*Cf))
# the cut-off frequency of the filter circuit

println("\n_____________________________________________________________________")
println("AC-DC Converter - Design Parameters\n")

println("RMS output voltage, Vo: ", Vo, " V")
println("Frequency of the reference/control waveform, fr: ", fr, " Hz")
println("Switching frequency of the triangle waveform, fs: ", fs, " Hz")
println("Output load resistor: ", RL, " Ω")
println("Output load inductor: ", Ll, " H")
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

# Simulation
Nps = (1/(Timestep*1e-6)) - 1 # time intervals, samples in a second + correction
μps = 1/Nps #corrected time step
t_final = (t_final - μps)
t = 0:μps:t_final # time

N = length(t) # number of samples

# Buck Boost control
Saw = Array{Float64, 1}(undef, N)
ControlB = Array{Int64, 1}(undef, N)

Saw = fill!(Saw, 0.0)
ControlB = fill!(ControlB, 0)

# DC_AC control
Rf_a = Array{Float64, 1}(undef, N)
Rf_b = Array{Float64, 1}(undef, N)
Rf_c = Array{Float64, 1}(undef, N)
Tri_a = Array{Float64, 1}(undef, N)
Tri_b = Array{Float64, 1}(undef, N)
Tri_c = Array{Float64, 1}(undef, N)
ControlA_a = Array{Int64, 1}(undef, N)
ControlA_b = Array{Int64, 1}(undef, N)
ControlA_c = Array{Int64, 1}(undef, N)
ControlA_a = fill!(ControlA_a, 0)
ControlA_b = fill!(ControlA_b, 0)
ControlA_c = fill!(ControlA_c, 0)
toggl_a = 1
toggl_b = 1
toggl_c = 1
dead_count_a = 0
dead_count_b = 0
dead_count_c = 0

# state space variables

x = Array{Float64, 2}(undef, 12, N)
x = fill!(x, 0)

# x = [iLB; vCBp; vCBn; iLf_a; vCf_a; iLf_b;
# vCf_b; iLf_c; vCf_c; io_a; io_b; io_c]

dims = 12

# output variables
y = Array{Float64, 2}(undef, dims, N)
y = fill!(y, 0)

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

A = Array{Float64, 4}(undef, 2, 8, dims, dims)
B = Array{Float64, 2}(undef, 2, dims)
C = Array{Float64, 4}(undef, 2, 8, dims, dims)
D = Array{Float64, 2}(undef, 2, dims)

CB = CB/2 #removing the two to avoid carrying it into matrices
# initialise
# Td is closed
# TA+ closed, TB+ closed, TC+ closed
A[1,1,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 -1/CB 0 -1/CB 0 -1/CB 0 0 0 0];#vCBp
                [0 0 0 0 0 0 0 0 0 0 0 0];#vCBn - should all be 0, not 1, but can invert now
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA+ closed, TB+ closed, TC- closed
A[1,2,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 -1/CB 0 -1/CB 0 0 0 0 0 0];#vCBp
                [0 0 0 0 0 0 0 -1/CB 0 0 0 0];#vCBn
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA+ closed, TB- closed, TC+ closed
A[1,3,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 -1/CB 0 0 0 -1/CB 0 0 0 0];#vCBp
                [0 0 0 0 0 -1/CB 0 0 0 0 0 0];#vCBn
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA+ closed, TB- closed, TC- closed
A[1,4,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 -1/CB 0 0 0 0 0 0 0 0];#vCBp
                [0 0 0 0 0 -1/CB 0 -1/CB 0 0 0 0];#vCBn
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA- closed, TB+ closed, TC+ closed
A[1,5,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 0 0 -1/CB 0 -1/CB 0 0 0 0];#vCBp
                [0 0 0 -1/CB 0 0 0 0 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA- closed, TB+ closed, TC- closed
A[1,6,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 0 0 -1/CB 0 0 0 0 0 0];#vCBp
                [0 0 0 -1/CB 0 0 0 -1/CB 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA- closed, TB- closed, TC+ closed
A[1,7,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 0 0 0 0 -1/CB 0 0 0 0];#vCBp
                [0 0 0 -1/CB 0 -1/CB 0 0 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA- closed, TB- closed, TC- closed
A[1,8,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 0 0 0 0 0 0 0 0 0];#vCBp
                [0 0 0 -1/CB 0 -1/CB 0 -1/CB 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c

# Td is open
# TA+ closed, TB+ closed, TC+ closed
A[2,1,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 -1/CB 0 -1/CB 0 -1/CB 0 0 0 0];#vCBp
                [1/CB 0 0 0 0 0 0 0 0 0 0 0];#vCBn - should all be 0, not 1, but can invert now
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA+ closed, TB+ closed, TC- closed
A[2,2,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 -1/CB 0 -1/CB 0 0 0 0 0 0];#vCBp
                [1/CB 0 0 0 0 0 0 -1/CB 0 0 0 0];#vCBn
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA+ closed, TB- closed, TC+ closed
A[2,3,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 -1/CB 0 0 0 -1/CB 0 0 0 0];#vCBp
                [1/CB 0 0 0 0 -1/CB 0 0 0 0 0 0];#vCBn
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA+ closed, TB- closed, TC- closed
A[2,4,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 -1/CB 0 0 0 0 0 0 0 0];#vCBp
                [1/CB 0 0 0 0 -1/CB 0 -1/CB 0 0 0 0];#vCBn
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA- closed, TB+ closed, TC+ closed
A[2,5,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 0 0 -1/CB 0 -1/CB 0 0 0 0];#vCBp
                [1/CB 0 0 -1/CB 0 0 0 0 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA- closed, TB+ closed, TC- closed
A[2,6,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 0 0 -1/CB 0 0 0 0 0 0];#vCBp
                [1/CB 0 0 -1/CB 0 0 0 -1/CB 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA- closed, TB- closed, TC+ closed
A[2,7,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 0 0 0 0 -1/CB 0 0 0 0];#vCBp
                [1/CB 0 0 -1/CB 0 -1/CB 0 0 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c
# TA- closed, TB- closed, TC- closed
A[2,8,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 0 0 0 0 0 0 0 0 0];#vCBp
                [1/CB 0 0 -1/CB 0 -1/CB 0 -1/CB 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -RL/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -RL/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -RL/Ll]]#io_c

#first order liner non-homogeneous equation
B[1,:] = [1/LB; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]

#first order liner homogeneous equation - no particular solution
B[2,:] = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]

# Td is closed
# TA+ closed, TB+ closed, TC+ closed
C[1,1,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 -1/CB 0 -1/CB 0 -1/CB 0 0 0 0];#vCBp
                [1 1 1 1 1 1 1 1 1 1 1 1];#vCBn - should all be 0, not 1, but can invert now
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA+ closed, TB+ closed, TC- closed
C[1,2,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 -1/CB 0 -1/CB 0 0 0 0 0 0];#vCBp
                [0 0 0 0 0 0 0 -1/CB 0 0 0 0];#vCBn
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA+ closed, TB- closed, TC+ closed
C[1,3,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 -1/CB 0 0 0 -1/CB 0 0 0 0];#vCBp
                [0 0 0 0 0 -1/CB 0 0 0 0 0 0];#vCBn
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA+ closed, TB- closed, TC- closed
C[1,4,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 -1/CB 0 0 0 0 0 0 0 0];#vCBp
                [0 0 0 0 0 -1/CB 0 -1/CB 0 0 0 0];#vCBn
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA- closed, TB+ closed, TC+ closed
C[1,5,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 0 0 -1/CB 0 -1/CB 0 0 0 0];#vCBp
                [0 0 0 -1/CB 0 0 0 0 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA- closed, TB+ closed, TC- closed
C[1,6,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 0 0 -1/CB 0 0 0 0 0 0];#vCBp
                [0 0 0 -1/CB 0 0 0 -1/CB 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA- closed, TB- closed, TC+ closed
C[1,7,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [0 0 0 0 0 0 0 -1/CB 0 0 0 0];#vCBp
                [0 0 0 -1/CB 0 -1/CB 0 0 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA- closed, TB- closed, TC- closed
C[1,8,:,:] =   [[(-RB/LB) 0 0 0 0 0 0 0 0 0 0 0];#iLB
                [1 1 1 1 1 1 1 1 1 1 1 1];#vCBp
                [0 0 0 -1/CB 0 -1/CB 0 -1/CB 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c

# Td is open
# TA+ closed, TB+ closed, TC+ closed
C[2,1,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 -1/CB 0 -1/CB 0 -1/CB 0 0 0 0];#vCBp
                [1/CB 0 0 0 0 0 0 0 0 0 0 0];#vCBn - should all be 0, not 1, but can invert now
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA+ closed, TB+ closed, TC- closed
C[2,2,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 -1/CB 0 -1/CB 0 0 0 0 0 0];#vCBp
                [1/CB 0 0 0 0 0 0 -1/CB 0 0 0 0];#vCBn
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA+ closed, TB- closed, TC+ closed
C[2,3,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 -1/CB 0 0 0 -1/CB 0 0 0 0];#vCBp
                [1/CB 0 0 0 0 -1/CB 0 0 0 0 0 0];#vCBn
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA+ closed, TB- closed, TC- closed
C[2,4,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 -1/CB 0 0 0 0 0 0 0 0];#vCBp
                [1/CB 0 0 0 0 -1/CB 0 -1/CB 0 0 0 0];#vCBn
                [0 1/Lf 0 -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA- closed, TB+ closed, TC+ closed
C[2,5,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 0 0 -1/CB 0 -1/CB 0 0 0 0];#vCBp
                [1/CB 0 0 -1/CB 0 0 0 0 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA- closed, TB+ closed, TC- closed
C[2,6,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 0 0 -1/CB 0 0 0 0 0 0];#vCBp
                [1/CB 0 0 -1/CB 0 0 0 -1/CB 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 1/Lf 0 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA- closed, TB- closed, TC+ closed
C[2,7,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 0 0 0 0 -1/CB 0 0 0 0];#vCBp
                [1/CB 0 0 -1/CB 0 -1/CB 0 0 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 1/Lf 0 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c
# TA- closed, TB- closed, TC- closed
C[2,8,:,:] =   [[(-RB/LB) 1/LB -1/LB 0 0 0 0 0 0 0 0 0];#iLB
                [-1/CB 0 0 0 0 0 0 0 0 0 0 0];#vCBp
                [1/CB 0 0 -1/CB 0 -1/CB 0 -1/CB 0 0 0 0];#vCBn
                [0 0 1/Lf -Rf/Lf -1/Lf 0 0 0 0 0 0 0];#iLf_a
                [0 0 0 1/Cf 0 0 0 0 0 -1/Cf 0 0];#vCf_a
                [0 0 1/Lf 0 0 -Rf/Lf -1/Lf 0 0 0 0 0];#iLf_b
                [0 0 0 0 0 1/Cf 0 0 0 0 -1/Cf 0];#vCf_b
                [0 0 1/Lf 0 0 0 0 -Rf/Lf -1/Lf 0 0 0];#iLf_c
                [0 0 0 0 0 0 0 1/Cf 0 0 0 -1/Cf];#vCf_c
                [0 0 0 0 1/Ll 0 0 0 0 -1/Ll 0 0];#io_a
                [0 0 0 0 0 0 1/Ll 0 0 0 -1/Ll 0];#io_b
                [0 0 0 0 0 0 0 0 1/Ll 0 0 -1/Ll]]#io_c

D[1,:] = [1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
D[2,:] = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]


#%% Discrete time matrices
#______________________________________________________________________________

Ad = Array{Float64, 4}(undef, 2, 8, dims, dims)
Bd = Array{Float64, 3}(undef, 2, 8, dims)

for i in 1:2
    for j in 2:7
        Ad[i,j,:,:] = exp(A[i,j,:,:]*μps)
        Bd[i,j,:] = inv(A[i,j,:,:])*(Ad[i,j,:,:] - Matrix(I, dims, dims))*B[i,:]
    end
end


#%% Starting time simulation
#______________________________________________________________________________
println()
for i in 1:N

    if i > 1 && floor((10*t[i]/t_final)) != floor((10*t[i-1]/t_final))
        println("Progress : ", floor((100*t[i]/t_final)), " %")
    end

    global x, y

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
    global toggl_a, toggl_b, toggl_c
    global dead_count_a, dead_count_b, dead_count_c

    Rf_a[i] = fr_p*sin(2*pi*fr*t[i])
    Rf_b[i] = fr_p*sin(2*pi*fr*t[i] - 120*π/180)
    Rf_c[i] = fr_p*sin(2*pi*fr*t[i] + 120*π/180)

    # Phase a
    #___________________________________________________________________________
    if Tri_a[i] >= fs_p && toggl_a == 1
        toggl_a = -1
    elseif Tri_a[i] <= -fs_p && toggl_a == -1
        toggl_a = 1
    end

    if i < N
        Tri_a[i+1] = toggl_a*ms/Nps + Tri_a[i]
    end

    if Rf_a[i] >= Tri_a[i] && i > 1

        if ControlA_a[i-1] == -1
            dead_count_a = t[i]
        end

        if t[i] - dead_count_a < tdΔ
            ControlA_a[i] = 0
        else
            ControlA_a[i] = 1
        end
    elseif Rf_a[i] < Tri_a[i] && i > 1

        if ControlA_a[i-1] == 1
            dead_count_a = t[i]
        end

        if t[i] - dead_count_a < tdΔ
            ControlA_a[i] = -0
        else
            ControlA_a[i] = -1
        end
    end

    # Phase b
    #___________________________________________________________________________
    if Tri_b[i] >= fs_p && toggl_b == 1
        toggl_b = -1
    elseif Tri_b[i] <= -fs_p && toggl_b == -1
        toggl_b = 1
    end

    if i < N
        Tri_b[i+1] = toggl_b*ms/Nps + Tri_b[i]
    end

    if Rf_b[i] >= Tri_b[i] && i > 1

        if ControlA_b[i-1] == -1
            dead_count_b = t[i]
        end

        if t[i] - dead_count_b < tdΔ
            ControlA_b[i] = 0
        else
            ControlA_b[i] = 1
        end
    elseif Rf_b[i] < Tri_b[i] && i > 1

        if ControlA_b[i-1] == 1
            dead_count_b = t[i]
        end

        if t[i] - dead_count_b < tdΔ
            ControlA_b[i] = -0
        else
            ControlA_b[i] = -1
        end
    end

    # Phase c
    #___________________________________________________________________________
    if Tri_c[i] >= fs_p && toggl_c == 1
        toggl_c = -1
    elseif Tri_c[i] <= -fs_p && toggl_c == -1
        toggl_c = 1
    end

    if i < N
        Tri_c[i+1] = toggl_c*ms/Nps + Tri_c[i]
    end

    if Rf_c[i] >= Tri_c[i] && i > 1

        if ControlA_c[i-1] == -1
            dead_count_c = t[i]
        end

        if t[i] - dead_count_c < tdΔ
            ControlA_c[i] = 0
        else
            ControlA_c[i] = 1
        end
    elseif Rf_c[i] < Tri_c[i] && i > 1

        if ControlA_c[i-1] == 1
            dead_count_c = t[i]
        end

        if t[i] - dead_count_c < tdΔ
            ControlA_c[i] = -0
        else
            ControlA_c[i] = -1
        end
    end

    # by setting this time and ControlA value we set when the AC conversion starts
    if t[i] < 0
        ControlA_a[i] = -1 # 0. 1, or -1
        ControlA_b[i] = -1 # 0. 1, or -1
        ControlA_c[i] = -1 # 0. 1, or -1
    end
    # System dynamics
    #___________________________________________________________________________
    if i < N

        if ControlB[i] == 1 && i < N # if Td is closed

            if ControlA_a[i] == 1 && ControlA_b[i] == 1 && ControlA_c[i] == 1

                #x[:,i+1] = Continuous_time(A[1,1,:,:], B[1,:], Vd, t[i], μps, x[:,i])
                #x[:,i+1] = Discrete_time(Ad[1,1,:,:], Bd[1,1,:], Vd, 1, x[:,i])

                dx_dt = A[1,1,:,:]*x[:,i]+ Vd*B[1,:]
                x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[1,1,:,:]*x[:,i+1] + Vd .* D[1,:]

            elseif ControlA_a[i] == 1 && ControlA_b[i] == 1 && ControlA_c[i] == -1

                #x[:,i+1] = Continuous_time(A[1,2,:,:], B[1,:], Vd, t[i], μps, x[:,i])
                x[:,i+1] = Discrete_time(Ad[1,2,:,:], Bd[1,2,:], Vd, 1, x[:,i])

                #dx_dt = A[1,2,:,:]*x[:,i]+ Vd*B[1,:]
                #x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[1,2,:,:]*x[:,i+1] + Vd .* D[1,:]

            elseif ControlA_a[i] == 1 && ControlA_b[i] == -1 && ControlA_c[i] == 1

                #x[:,i+1] = Continuous_time(A[1,3,:,:], B[1,:], Vd, t[i], μps, x[:,i])
                x[:,i+1] = Discrete_time(Ad[1,3,:,:], Bd[1,3,:], Vd, 1, x[:,i])

                #dx_dt = A[1,3,:,:]*x[:,i]+ Vd*B[1,:]
                #x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[1,3,:,:]*x[:,i+1] + Vd .* D[1,:]

            elseif ControlA_a[i] == 1 && ControlA_b[i] == -1 && ControlA_c[i] == -1

                #x[:,i+1] = Continuous_time(A[1,4,:,:], B[1,:], Vd, t[i], μps, x[:,i])
                x[:,i+1] = Discrete_time(Ad[1,4,:,:], Bd[1,4,:], Vd, 1, x[:,i])

                #dx_dt = A[1,4,:,:]*x[:,i]+ Vd*B[1,:]
                #x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[1,4,:,:]*x[:,i+1] + Vd .* D[1,:]

            elseif ControlA_a[i] == -1 && ControlA_b[i] == 1 && ControlA_c[i] == 1

                #x[:,i+1] = Continuous_time(A[1,5,:,:], B[1,:], Vd, t[i], μps, x[:,i])
                x[:,i+1] = Discrete_time(Ad[1,5,:,:], Bd[1,5,:], Vd, 1, x[:,i])

                #dx_dt = A[1,5,:,:]*x[:,i]+ Vd*B[1,:]
                #x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[1,5,:,:]*x[:,i+1] + Vd .* D[1,:]

            elseif ControlA_a[i] == -1 && ControlA_b[i] == 1 && ControlA_c[i] == -1

                #x[:,i+1] = Continuous_time(A[1,6,:,:], B[1,:], Vd, t[i], μps, x[:,i])
                x[:,i+1] = Discrete_time(Ad[1,6,:,:], Bd[1,6,:], Vd, 1, x[:,i])

                #dx_dt = A[1,6,:,:]*x[:,i]+ Vd*B[1,:]
                #x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[1,6,:,:]*x[:,i+1] + Vd .* D[1,:]

            elseif ControlA_a[i] == -1 && ControlA_b[i] == -1 && ControlA_c[i] == 1

                #x[:,i+1] = Continuous_time(A[1,7,:,:], B[1,:], Vd, t[i], μps, x[:,i])
                x[:,i+1] = Discrete_time(Ad[1,7,:,:], Bd[1,7,:], Vd, 1, x[:,i])

                #dx_dt = A[1,7,:,:]*x[:,i]+ Vd*B[1,:]
                #x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[1,7,:,:]*x[:,i+1] + Vd .* D[1,:]

            elseif ControlA_a[i] == -1 && ControlA_b[i] == -1 && ControlA_c[i] == -1

                #x[:,i+1] = Continuous_time(A[1,8,:,:], B[1,:], Vd, t[i], μps, x[:,i])
                #x[:,i+1] = Discrete_time(Ad[1,8,:,:], Bd[1,8,:], Vd, 1, x[:,i])

                dx_dt = A[1,8,:,:]*x[:,i]+ Vd*B[1,:]
                x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[1,8,:,:]*x[:,i+1] + Vd .* D[1,:]

            elseif ControlA_a[i] == 0 && ControlA_b[i] == 0 && ControlA_c[i] == 0

                if 5 >= 0

                    #x = Continuous_time(A01, B01, Vd, t[i], μps, x0)
                    #x = Discrete_time(Ad01, Bd01, Vd, 1, x0)

                    #  Outputs
                    #y = C01*x + D01*Vd
                else

                    #x = Continuous_time(A02, B02, Vd, t[i], μps, x0)
                    #x = Discrete_time(Ad02, Bd02, Vd, 1, x0)

                    #  Outputs
                    #y = C02*x + D02*Vd
                end
            end

        elseif ControlB[i] == 0 && i < N # if Td is open

            if ControlA_a[i] == 1 && ControlA_b[i] == 1 && ControlA_c[i] == 1

                #x[:,i+1] = Continuous_time(A[2,1,:,:], B[2,:], Vd, t[i], μps, x[:,i])
                #x[:,i+1] = Discrete_time(Ad[2,1,:,:], Bd[2,1,:], Vd, 1, x[:,i])

                dx_dt = A[2,1,:,:]*x[:,i]+ Vd*B[2,:]
                x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[2,1,:,:]*x[:,i+1] + Vd .* D[2,:]

            elseif ControlA_a[i] == 1 && ControlA_b[i] == 1 && ControlA_c[i] == -1

                #x[:,i+1] = Continuous_time(A[2,2,:,:], B[2,:], Vd, t[i], μps, x[:,i])
                x[:,i+1] = Discrete_time(Ad[2,2,:,:], Bd[2,2,:], Vd, 1, x[:,i])

                #dx_dt = A[2,2,:,:]*x[:,i]+ Vd*B[2,:]
                #x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[2,2,:,:]*x[:,i+1] + Vd .* D[2,:]

            elseif ControlA_a[i] == 1 && ControlA_b[i] == -1 && ControlA_c[i] == 1

                #x[:,i+1] = Continuous_time(A[2,3,:,:], B[2,:], Vd, t[i], μps, x[:,i])
                x[:,i+1] = Discrete_time(Ad[2,3,:,:], Bd[2,3,:], Vd, 1, x[:,i])

                #dx_dt = A[2,3,:,:]*x[:,i]+ Vd*B[2,:]
                #x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[2,3,:,:]*x[:,i+1] + Vd .* D[2,:]

            elseif ControlA_a[i] == 1 && ControlA_b[i] == -1 && ControlA_c[i] == -1

                #x[:,i+1] = Continuous_time(A[2,4,:,:], B[2,:], Vd, t[i], μps, x[:,i])
                x[:,i+1] = Discrete_time(Ad[2,4,:,:], Bd[2,4,:], Vd, 1, x[:,i])

                #dx_dt = A[2,4,:,:]*x[:,i]+ Vd*B[2,:]
                #x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[2,4,:,:]*x[:,i+1] + Vd .* D[2,:]

            elseif ControlA_a[i] == -1 && ControlA_b[i] == 1 && ControlA_c[i] == 1

                #x[:,i+1] = Continuous_time(A[2,5,:,:], B[2,:], Vd, t[i], μps, x[:,i])
                x[:,i+1] = Discrete_time(Ad[2,5,:,:], Bd[2,5,:], Vd, 1, x[:,i])

                #dx_dt = A[2,5,:,:]*x[:,i]+ Vd*B[2,:]
                #x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[2,5,:,:]*x[:,i+1] + Vd .* D[2,:]

            elseif ControlA_a[i] == -1 && ControlA_b[i] == 1 && ControlA_c[i] == -1

                #x[:,i+1] = Continuous_time(A[2,6,:,:], B[2,:], Vd, t[i], μps, x[:,i])
                x[:,i+1] = Discrete_time(Ad[2,6,:,:], Bd[2,6,:], Vd, 1, x[:,i])

                #dx_dt = A[2,6,:,:]*x[:,i]+ Vd*B[2,:]
                #x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[2,6,:,:]*x[:,i+1] + Vd .* D[2,:]

            elseif ControlA_a[i] == -1 && ControlA_b[i] == -1 && ControlA_c[i] == 1

                #x[:,i+1] = Continuous_time(A[2,7,:,:], B[2,:], Vd, t[i], μps, x[:,i])
                x[:,i+1] = Discrete_time(Ad[2,7,:,:], Bd[2,7,:], Vd, 1, x[:,i])

                #dx_dt = A[2,7,:,:]*x[:,i]+ Vd*B[2,:]
                #x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[2,7,:,:]*x[:,i+1] + Vd .* D[2,:]

            elseif ControlA_a[i] == -1 && ControlA_b[i] == -1 && ControlA_c[i] == -1

                #x[:,i+1] = Continuous_time(A[2,8,:,:], B[2,:], Vd, t[i], μps, x[:,i])
                #x[:,i+1] = Discrete_time(Ad[2,8,:,:], Bd[2,8,:], Vd, 1, x[:,i])

                dx_dt = A[2,8,:,:]*x[:,i]+ Vd*B[2,:]
                x[:,i+1] = x[:,i] + μps.*dx_dt

                #  Outputs
                y[:,i+1] = C[2,8,:,:]*x[:,i+1] + Vd .* D[2,:]

            elseif ControlA_a[i] == 0 && ControlA_b[i] == 0 && ControlA_c[i] == 0

                if 5 >= 0

                    #x = Continuous_time(A01, B01, Vd, t[i], μps, x0)
                    #x = Discrete_time(Ad01, Bd01, Vd, 1, x0)

                    #  Outputs
                    #y = C01*x + D01*Vd
                else

                    #x = Continuous_time(A02, B02, Vd, t[i], μps, x0)
                    #x = Discrete_time(Ad02, Bd02, Vd, 1, x0)

                    #  Outputs
                    #y = C02*x + D02*Vd
                end
            end
        end

        if x[1, i+1] < 0 # boundary between continuous and discontinuous current
            x[1, i+1] = 0
        end
    end
end
println("Progress : 100.0 %")

#_______________________________________________________________________________
T_B = t_final*fB # periods for the buck boost stage, i.e. simulation time
T_r = t_final*fr # periods for the reference stage, i.e. simulation time

# Harmonic Analysis
#______________________________________________________________________________

per = 3
T_plot_start = T_r-per
N_plot_start = convert(Int64, round((T_plot_start/fr)*Nps))
N_dif = N - N_plot_start

shifted_k = fftshift(fftfreq(length(N_plot_start:N), Nps))
VCf_shifted_fft = 2/N_dif*fftshift(fft(x[5, N_plot_start:N]))
ILf_shifted_fft = 2/N_dif*fftshift(fft(x[4, N_plot_start:N]))

k_start = length(VCf_shifted_fft)/2 + 1
k_start = convert(Int64, round(k_start))
k_end = length(VCf_shifted_fft)
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

vCB = x[3, :] - x[2, :]
ΔvCB = maximum(vCB[N_plot_eval : N]) - minimum(vCB[N_plot_eval : N])
ΔiLB = maximum(x[1, N_plot_eval : N]) - minimum(x[1, N_plot_eval : N])
VCB_average = sum(vCB[N_plot_eval : N])/length(vCB[N_plot_eval : N])
iLB_average = sum(x[1, N_plot_eval : N])/length(x[1, N_plot_eval : N])

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

vo_average = sum(x[5, N_plot_eval : N])/length(x[5, N_plot_eval : N])
vo_peak = maximum(x[5, N_plot_eval : N])
vo_rms = std(x[5, N_plot_eval : N])

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

d = (x[5, N_plot_eval_bi : N_plot_eval_bii] - VCf_signal[N_plot_eval_bi : N_plot_eval_bii])
e = (x[4, N_plot_eval_bi : N_plot_eval_bii] - ILf_signal[N_plot_eval_bi : N_plot_eval_bii])

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
round(100*Δvo/vo_peak, digits = 2), " %")
println("Maximum Ripple current to Inductor current, ΔILF/ILF: ",
round(100*ΔiLf*RL/vo_peak, digits = 2), " %")
println("Voltage THD: ",
round(THD_vCf, digits = 2), " %")

# PWM
#_______________________________________________________________________________

FFTs_p = plot(freq, VCf_harm_mag, seriestype = :scatter, line=:stem, title = "Frequency Spectrum",
    label = "Shifted FFT",
    markershape = :circle, markersize = 2,
    xlim = (0, +2*fs),
    xticks = 0:500:2*fs,
    legend = true);

display(FFTs_p)

p2 = plot(t[N_plot_eval_bi : N_plot_eval_bii], e, title = "Converter Ripple Current",
    label = "iLf_a",
    xlabel = "Time [s]",
    ylabel = "Current [A]")

display(p2)

p3 = plot(t[N_plot_eval_bi : N_plot_eval_bii], d, title = "Converter Ripple Voltage",
    label = "VLf_a",
    xlabel = "Time [s]",
    ylabel = "Electrical Potential [V]")

display(p3)

T_plot1 = 3
N_plot1 = convert(Int64, round((T_plot1/fB  - 1/Nps)*Nps))
p_cntr_BB = plot(t[1:N_plot1], Saw[1:N_plot1], label = "Saw",
xlabel = "Time [s]",
title = "Buck Boost Control")
p_cntr_BB = plot!(t[1:N_plot1], cntr_B*ones(N_plot1), label = "vc")
p_cntr_BB = plot!(t[1:N_plot1], ControlB[1:N_plot1], label = "Control")

T_plot2 = 3
N_plot2 = convert(Int64, round((T_plot2/fs  - 1/Nps)*Nps))
# Phase a Control Signals
p_cntr_AC_a = plot(t[1:N_plot2], Rf_a[1:N_plot2], label = "Reference_a",
xlabel = "Time [s]",
title = "DC-AC Converter Control Phase A")
p_cntr_AC_a = plot!(t[1:N_plot2], Tri_a[1:N_plot2], label = "Triangle_a")
p_cntr_AC_a = plot!(t[1:N_plot2], ControlA_a[1:N_plot2], label = "Control_a")

# Phase b Control Signals
p_cntr_AC_b = plot(t[1:N_plot2], Rf_b[1:N_plot2], label = "Reference_b",
xlabel = "Time [s]",
title = "DC-AC Converter Control Phase B")
p_cntr_AC_b = plot!(t[1:N_plot2], Tri_b[1:N_plot2], label = "Triangle_b")
p_cntr_AC_b = plot!(t[1:N_plot2], ControlA_b[1:N_plot2], label = "Control_b")

# Phase c Control Signals
p_cntr_AC_c = plot(t[1:N_plot2], Rf_c[1:N_plot2], label = "Reference_c",
xlabel = "Time [s]",
title = "DC-AC Converter Control Phase C")
p_cntr_AC_c = plot!(t[1:N_plot2], Tri_c[1:N_plot2], label = "Triangle_c")
p_cntr_AC_c = plot!(t[1:N_plot2], ControlA_c[1:N_plot2], label = "Control_c")

PWM = plot(p_cntr_BB, p_cntr_AC_a, p_cntr_AC_b, p_cntr_AC_c,
layout = (4, 1),
legend = true,
size=(900,900))
display(PWM)

#%% Plots
# Buck Boost Voltage
#_______________________________________________________________________________
T_plot_start = 0
T_plot_end = T_r
N_plot_start = convert(Int64, round((T_plot_start/fr  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fr  - 1/Nps)*Nps))

vCB = x[3,N_plot_start : N_plot_end] - x[2,N_plot_start : N_plot_end]
p_V_BB = plot(t[N_plot_start : N_plot_end], vCB[N_plot_start : N_plot_end],
title = "Buck Boost Output Voltage",
label = false,
xlabel = "Time [s]",
ylabel = "Electrical Potential [V]",
grid = true,
foreground_color_grid = :black,
minorgrid = true,
thickness_scaling = 1.5,
legendfont = font(5),
size = (800, 600))
display(p_V_BB)

# DC AC Converter Voltage
#_______________________________________________________________________________
T_plot_start = 0
T_plot_end = T_r
N_plot_start = convert(Int64, round((T_plot_start/fr  + 1/Nps)*Nps))
N_plot_end = convert(Int64, round((T_plot_end/fr  - 1/Nps)*Nps))

p_V_AC_inst = plot(t[N_plot_start : N_plot_end], x[5,N_plot_start : N_plot_end],
title = "DC-AC Converter Instantaneous Voltage",
label = "Phase A",
xlabel = "Time [s]",
ylabel = "Electrical Potential [V]",
grid = true,
foreground_color_grid = :black,
minorgrid = true,
thickness_scaling = 1.5,
legendfont = font(5))
p_V_AC_inst = plot!(t[N_plot_start : N_plot_end], x[7,N_plot_start : N_plot_end],
label = "Phase B")
p_V_AC_inst = plot!(t[N_plot_start : N_plot_end], x[9,N_plot_start : N_plot_end],
label = "Phase C")

# Calculating Moving Window RMS values
#_______________________________________________________________________

T_rms = T_plot_start:0.1:Int(floor(T_plot_end))

V_rms = Array{Float64, 2}(undef, 3, length(T_rms))
V_rms = fill!(V_rms, 0)
I_rms = Array{Float64, 2}(undef, 3, length(T_rms))
I_rms = fill!(I_rms, 0)

for i in 1:length(T_rms)

    N_end = convert(Int64, round((T_rms[i]*tr)*Nps))
    T_eval = 1 # number of cycles to use for calculation
    N_start = N_end - convert(Int64, round((T_eval*tr)*Nps))

    if N_start < 1
        N_start = 1
    end

    V_rms[1,i] = std(x[5, N_start : N_end])
    V_rms[2,i] = std(x[7, N_start : N_end])
    V_rms[3,i] = std(x[9, N_start : N_end])
    I_rms[1,i] = std(x[10, N_start : N_end])
    I_rms[2,i] = std(x[11, N_start : N_end])
    I_rms[3,i] = std(x[12, N_start : N_end])
end

p_V_AC_rms = plot(tr*T_rms, V_rms[1,:],
title = "DC-AC Converter RMS Voltage",
label = "Phase A",
xlabel = "Time [s]",
ylabel = "Electrical Potential [Vrms]",
grid = true,
foreground_color_grid = :black,
minorgrid = true,
thickness_scaling = 1.5,
legendfont = font(5))
p_V_AC_rms = plot!(tr*T_rms, V_rms[2,:],
label = "Phase B")
p_V_AC_rms = plot!(tr*T_rms, V_rms[3,:],
label = "Phase C")

p_I_AC_rms = plot(tr*T_rms, I_rms[1,:],
title = "DC-AC Converter RMS Voltage",
label = "Phase A",
xlabel = "Time [s]",
ylabel = "Current [A rms]",
grid = true,
foreground_color_grid = :black,
minorgrid = true,
thickness_scaling = 1.5,
legendfont = font(5))
p_V_AC_rms = plot!(tr*T_rms, I_rms[2,:],
label = "Phase B")
p_V_AC_rms = plot!(tr*T_rms, I_rms[3,:],
label = "Phase C")

PWM = plot(p_V_AC_inst, p_V_AC_rms, p_I_AC_rms,
layout = (3, 1),
legend = true,
size = (900,900))
display(PWM)

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
