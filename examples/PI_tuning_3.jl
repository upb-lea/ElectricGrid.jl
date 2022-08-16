using Plots
using LinearAlgebra
using FFTW
using SymPy
using ControlSystems

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

# Parameters
Lf = 2.3e-3 #Source.Lf[1] #
Rf = 0.4
Cf = 10e-6 #Source.Cf[1] #
Ts = 100e-6 #1/Source.f_cntr #
Vdc = 1200 #Source.Vdc[1] #

ωc = 1/(sqrt(Lf*Cf)) # filter cut-off
println("\nFilter Cut-off, ωc = ", ωc, " rad/s")

#= Controller Definitions
    ωcp - gain crossover frequency, is the frequency where the phase margin is measured,
    which is a 0-dB Gain crossing frequency
    ωcg - phase crossover frequency, is the frequency where the gain margin is measure,
    which is a -π phase crossing frequency
    pm - phase margin, the phase margin is the difference between the phase of the response
    and -π when the loop gain is 1
    gm - gain margin, the amount of gain variance required to make the loop gain unity
    at the frequency ωcg where the phase angle is -π. In other words, the gain margin
    is 1/g if g is the gain at the -π phase frequency. Negative gain margins indicate
    that stability is lost by decreasing the gain, while positive gain margins indicate
    that stability is lost by decreasing the gain, while positive gain margins indicate
    that stability is lost by increasing the gain.
=#

#= Dead Time in Digital Control Loops
    If the control scheme is implemented on a microcontroller or microprocessor,
    then a certain time is required to process the control algorithm. Therefore,
    a measured value can affect the voltage reference only after this time period
    has passed. In an appropriate manner, all these processes are synchronised with
    the clock cycle given by the pulse width modulation or vector modulation. This
    way, the digital control loop introduces a dead time of one sampling step.
    Together with the ZoH which samples the signals at the PoC, a total dead time of
    1.5 sampling steps of the current control loop results.
=#

#--------------------------------------
# Current Controller

println("\n#--------------------------------------")
println("Current Controller\n")
# Short Circuit State Space
A = [-Rf/Lf]
B = [Vdc/(2*Lf)]
C = [1]
D = [0]
sys_sc = ss(A,B,C,D) # continuous

dly = 0.5*Ts
ZoH = tf([1], [Ts/2, 1]) # Transfer function of the sample and hold process + Pade
Pade = tf([-dly/2, 1], [dly/2, 1]) # Pure first-order delay approximation

SC = tf([Vdc/(2*Lf)], [1, Rf/Lf]) # = tf(sys_sc)
Gsc_ol = tf(sys_sc)*Pade#*ZoH

for i in 5:20

    global Gsc_ol, Gi_cl, Gpi, pm, ωp

    ωp = 2π/(i*Ts) # gain cross-over frequency
    pm = 60 # degrees, phase margin
    Gpi, kp, ki = loopshapingPI(Gsc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)

    Gi_cl = G_PS(Gsc_ol*Gpi, tf(1)) # closed loop transfer function

    if any(real(poles(Gi_cl)) .> 0) == false

        # all the poles are on the left side
        ωp = 2π/((i+1)*Ts) # gain cross-over frequency
        pm = 60 # degrees, phase margin
        Gpi_i, kp_i, ki_i = loopshapingPI(Gsc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)
        Gi_cl = G_PS(Gsc_ol*Gpi_i, tf(1))

        println("Gain Cross-over frequency, ωp = ", ωp, " rad/s")
        println("Gain Cross-over frequency, fp = ", ωp/2π, " Hz")
        println("Gain Cross-over frequency, ωp = 2π/(", i+1,"*Ts) rad/s")
        println("Phase Margin, pm = ", pm, " degrees")
        println("\nPI coefficient, kp = ", kp_i, " V/A")
        println("PI coefficient, ki = ", ki_i, " V/As")

        println("\nClosed Loop Poles: ")
        for j in 1:length(poles(Gi_cl))
            println("   ", j,". = ", poles(Gi_cl)[j])
        end

        break

    end
end

p1 = bodeplot(Gi_cl);
p2 = nyquistplot(Gsc_ol*Gpi, ylims = (-1.5,1.5), xlims = (-1.5,1.5), unit_circle = true);
p3 = pzmap(Gi_cl);

y = plot(p1, p2, p3, layout=(3, 1), size=(800,800))
display(y)
#=
#--------------------------------------
# Voltage Controller
println("\n#--------------------------------------")
println("Voltage Controller\n")
# Open Circuit State Space
A = [[-Rf/Lf -1/Lf];
    [1/Cf 0]]
B = [Vdc/(2*Lf); 0]
C = [[0 1];
    [0 0]]
D = [0; 0]
sys_oc = ss(A,B,C,D) # continuous

OC = tf([Vdc/(2*Lf*Cf)], [1, Rf/Lf, 1/(Lf*Cf)]) # = tf(sys_oc)[1,1]
ZoH = tf([1], [Ts/2, 1]) # Transfer function of the sample and hold process + Pade
Goc_ol = tf(sys_oc)[1,1]*ZoH*Gi_cl

for i in 5:60

    global Goc_ol, Gv_cl, Gpi_v, pm, ωp

    ωp = 2π/(i*Ts) # gain cross-over frequency
    #println("Gain Cross-over frequency, fp = ", ωp/2π, " Hz")
    #fp = 300/(i-4)
    pm = 60 # degrees, phase margin
    Gpi_v, kp_v, ki_v = loopshapingPI(Goc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)

    Gv_cl = G_PS(Goc_ol*Gpi_v, tf(1)) # closed loop transfer function

    if any(real(poles(Gv_cl)) .> 0) == false

        # all the poles are on the left side
        ωp = 2π/((i+3)*Ts) # gain cross-over frequency
        pm = 60 # degrees, phase margin
        Gpi_v, kp_v, ki_v = loopshapingPI(Goc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)
        Gv_cl = G_PS(Goc_ol*Gpi_v, tf(1))

        println("Gain Cross-over frequency, ωp = ", ωp, " rad/s")
        println("Gain Cross-over frequency, fp = ", ωp/2π, " Hz")
        println("Gain Cross-over frequency, ωp = 2π/(", i+1,"*Ts) rad/s")
        println("Phase Margin, pm = ", pm, " degrees")
        println("\nPI coefficient, kp = ", kp_v, " A/V")
        println("PI coefficient, ki = ", ki_v, " A/Vs")

        println("\nClosed Loop Poles: ")
        for j in 1:length(poles(Gv_cl))
            println("   ",j,". = ",poles(Gv_cl)[j])
        end

        break

    end
end
=#

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
