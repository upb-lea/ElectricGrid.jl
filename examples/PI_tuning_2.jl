using Plots
using LinearAlgebra
using FFTW
using SymPy
using ControlSystems

#= Theory 1: Controller Definitions
    ωcp - gain crossover frequency, is the frequency where the phase margin is measured,
    which is a 0-dB Gain crossing frequency
    ωcg - phase crossover frequency, is the frequency where the gain margin is measured,
    which is a -π phase crossing frequency
    pm - phase margin, the phase margin is the difference between the phase of the response
    and -π when the loop gain is 1
    gm - gain margin, the amount of gain variance required to make the loop gain unity
    at the frequency ωcg where the phase angle is -π. In other words, the gain margin
    is 1/g if g is the gain at the -π phase frequency. Negative gain margins indicate
    that stability is lost by decreasing the gain, while positive gain margins indicate
    that stability is lost by increasing the gain.
=#

#= Theory 2: Dead Time in Digital Control Loops
    If the control scheme is implemented on a microcontroller or microprocessor,
    then a certain time is required to process the control algorithm. Therefore,
    a measured value can affect the voltage reference only after this time period
    has passed. In an appropriate manner, all these processes are synchronised with
    the clock cycle given by the pulse width modulation or vector modulation. This
    way, the digital control loop introduces a dead time of one sampling step.
    Together with the ZoH which samples the signals at the PoC, a total dead time of
    1.5 sampling steps of the current control loop results.
=#

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

# Parameters


Timestep = 100 #time step in μs ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
f_cntr = 1/(Timestep*1e-6) # Hz, Sampling frequency of controller ~ 15 kHz -> 50kHz

#num_source = 1

Lf = Source.Lf[num_source] # 2.3e-3 # Henry, filter inductor 
Rf = Source.Rf[num_source]  # 400e-3 # Ohms, inductor parasitic resistance
Cf = Source.Cf[num_source] # 10e-6  # Farad, filter capacitor
Ts = Timestep*1e-6 # 0.1e-3 # Seconds, switching/sampling time step
Vdc = Source.Vdc[num_source] # 2*600 # Volts, DC bus 

fsys = 50 # Hz, System frequency

min_fp = 300 # Hz, minimum allowable gain cross-over frequency

max_i = convert(Int64, floor(1/(min_fp*Ts)))
ωc = 1/(sqrt(Lf*Cf)) # filter cut-off
println("\nFilter Cut-off, ωc = ", round(ωc, digits = 3), " [rad/s], fc = ", 
round(ωc/(1e3*2π), digits = 3), " [kHz]")
println("\nInverter Switching Frequency, fs = ", round(1/(1e3*Ts), digits = 3), " [kHz]")

#--------------------------------------
# Current Controller

println("\n#--------------------------------------")
println("Current Controller")
println("Short Circuit State Space\n")
# Short Circuit State Space
A = [-Rf/Lf]
B = [1/(Lf)]
C = [1]
D = [0]
sys_sc = ss(A, B, C, D) # continuous

test = tf([1], [Ts/2, 1])*tf([1, 0], [1, 0])

dly = 1*Ts
ZoH = tf([1], [Ts/2, 1]) # Transfer function of the sample and hold process
Pade = tf([-dly/2, 1], [dly/2, 1]) # Pure first-order delay approximation
PWM_gain = Vdc/2

#SC = tf([1/(1*Lf)], [1, Rf/Lf]) # = tf(sys_sc)
Gsc_ol = minreal(tf(sys_sc)*Pade*PWM_gain*ZoH) # Full transfer function of plant
#Gsc_ol = minreal(tf(sys_sc)*Pade*PWM_gain) # Full transfer function of plant
kp_i = 0.0
ki_i = 0.0
for i in 6:max_i

    global Gsc_ol, Gi_cl, Gpi_i, ωp, pm, kp_i, ki_i 

    ωp = 2π/(i*Ts) # gain cross-over frequency
    pm = 60 # degrees, phase margin
    Gpi_i, _, _ = loopshapingPI(Gsc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)

    Gi_cl = minreal(Gsc_ol*Gpi_i/(1 + Gsc_ol*Gpi_i)) # closed loop transfer function

    if any(real(poles(Gi_cl)) .> 0) == false

        # all the poles are on the left side
        ωp = 2π/((i)*Ts) # gain cross-over frequency
        pm = 60 # degrees, phase margin
        Gpi_i, kp_i, ki_i = loopshapingPI(Gsc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)
        Gi_cl = minreal(Gsc_ol*Gpi_i/(1 + Gsc_ol*Gpi_i))

        println("Gain Cross-over frequency, ωp = ", round(ωp), " [rad/s]")
        println("Gain Cross-over frequency, fp = ", round(ωp/2π), " [Hz]")
        println("Gain Cross-over frequency, ωp = 2π/(", i,"*Ts) [rad/s]")
        println("Design Phase Margin, pm = ", pm, " [degrees]")

        Validate = evalfr(Gsc_ol*Gpi_i, 1im*(ωp))
        println("\nMeasured gain @ crossover = ", round(10*log(10, abs(Validate[1]))), " [dB]")
        println("Measured phase margin = ", round(180 + angle(Validate[1])*180/π), " [degrees]")

        println("\nPI coefficient, kp = ", round(kp_i, digits = 3), " [A/V]")
        println("PI coefficient, ki = ", round(ki_i, digits = 3), " [A/Vs]")

        println("\nClosed Loop Poles: ")
        for j in 1:length(poles(Gi_cl))
            println("   ", j,". = ", round(poles(Gi_cl)[j], digits = 3))
        end

        #=
            println("\nClosed Loop Zeros: ")
            for j in 1:length(tzeros(Gi_cl))
                println("   ", j,". = ", round(tzeros(Gi_cl)[j], digits = 3))
            end
        =#
        break

    end

    if i == max_i
        println("\nHouston, we have a problem.")
    end
end

p1 = bodeplot(Gi_cl);
p2 = nyquistplot(Gsc_ol*Gpi_i, ylims = (-1.5, 1.5), xlims = (-1.5, 1.5), unit_circle = true);
p3 = pzmap(Gi_cl);

y = plot(p1, p2, p3, layout=(3, 1), size=(800, 800))
#display(y)

ζ = 0.2
ω = 1

B = [1]
A = [1, 2ζ*ω, ω^2]
P = tf(B, A)
#z = plot(step(P, 15))
#z = plot(step(Gi_cl, 2.5))
#display(z)

#--------------------------------------
# Voltage Controller
println("\n#--------------------------------------")
println("Voltage Controller")
println("Open Circuit State Space\n")
# Open Circuit State Space
A = [[-Rf/Lf -1/Lf];
    [1/Cf 0]]
B = [1/(Lf); 0]
C = [[0 1];
    [1 0]]
D = [0; 0]
sys_oc = ss(A, B, C, D) # continuous

dly = 1*Ts
ZoH = tf([1], [Ts/2, 1]) # Transfer function of the sample and hold process
Pade = tf([-dly/2, 1], [dly/2, 1]) # Pure first-order delay approximation
PWM_gain = Vdc/2

Gpi_i = tf([kp_i, ki_i], [1, 0])
Gi_oc_ol = tf(sys_oc)[2,1]*Gpi_i*ZoH*Gpi_i*PWM_gain*Pade
Gi_oc_cl = minreal(Gi_oc_ol/(1 + Gi_oc_ol)) # closed loop transfer function of the current controller and open circuit plant response

# It might be a good idea to check that the current controller is still stable
println("Current Controller Closed Loop Poles: ")
for j in 1:length(poles(Gi_oc_cl))
    println("   ", j,". = ", round(poles(Gi_oc_cl)[j], digits = 3))
end

#OC = tf([1/(Lf*Cf)], [1, Rf/Lf, 1/(Lf*Cf)]) # = tf(sys_oc)[1,1] - basic check
Goc_ol = minreal(tf(sys_oc)[1,1]*ZoH*Gi_oc_cl) # Full open loop transfer function of current controller and open circuit plant

Goc_ol = minreal(Gi_cl*tf([1], [Cf, 0]))

for i in max_i:-1:1

    global Goc_ol, Gv_cl, Gpi_v, ωp, pm

    ωp = 2π/(i*Ts) # gain cross-over frequency
    pm = 60 # degrees, phase margin
    Gpi_v, _, _ = loopshapingPI(Goc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)

    Gv_cl = minreal(Goc_ol*Gpi_v/(1 + Goc_ol*Gpi_v)) # closed loop transfer function

    if any(real(poles(Gv_cl)) .> 0) == false

        # all the poles are on the left side
        ωp = 2π/((i)*Ts) # gain cross-over frequency
        pm = 60 # degrees, phase margin
        Gpi_v, kp_v, ki_v = loopshapingPI(Goc_ol, ωp, rl = 1, phasemargin = pm, form = :parallel)
        Gv_cl = minreal(Goc_ol*Gpi_v/(1 + Goc_ol*Gpi_v))

        println("\nGain Cross-over frequency, ωp = ", round(ωp, digits = 3), " [rad/s]")
        println("Gain Cross-over frequency, fp = ", round(ωp/2π, digits = 3), " [Hz]")
        println("Gain Cross-over frequency, ωp = 2π/(", i,"*Ts) [rad/s]")
        println("Phase Margin, pm = ", pm, " [degrees]")

        Validate = evalfr(Goc_ol*Gpi_v, 1im*(ωp))
        println("\nMeasured gain @ crossover = ", round(10*log(10, abs(Validate[1]))), " [dB]")
        println("Measured phase margin = ", round(180 + angle(Validate[1])*180/π), " [degrees]")
        
        println("\nPI coefficient, kp = ", round(kp_v, digits = 3), " [A/V]")
        println("PI coefficient, ki = ", round(ki_v, digits = 3), " [A/Vs]")

        println("\nClosed Loop Poles: ")
        for j in 1:length(poles(Gv_cl))
            println("   ", j,". = ", round(poles(Gv_cl)[j], digits = 3))
        end

        break

    end

    if i == 1
        println("\nHouston, we have a problem.")
    end
end

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
