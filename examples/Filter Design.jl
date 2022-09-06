print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

#____________________________________________________________
# Inputs

Ir = 15/100 # i_rip in %
Vr = 1.537/100 # v_rip in %

Sr = 150e3 # Srated, apparent power in VA
fs = 10e3 # Hz, switching frequency
Vdc = 800 # V, dc voltage

Vrms = 230 # rms voltage

println("\nHalf-bridge inverter LC filter design:\n")
println("Inverter rated power, ", Sr/1e3, " kVA")
println("Inverter switching frequency, ", fs/1e3, " kHz")
println("DC voltage, ", Vdc, " V")
println("AC POC rms voltage, ", Vrms, " V\n")
println("Inductor ripple current, ", Ir*100, " %")
println("Capacitor ripple voltage, ", Vr*100, " %\n")

#____________________________________________________________
# Inductor Design
Vorms = Vrms*1.05
Vop = Vorms*sqrt(2)

Zl = 3*Vorms*Vorms/Sr

Iorms = Vorms/Zl
Iop = Iorms*sqrt(2)

ΔIlfmax = Ir*Iop

Lf = Vdc/(4*fs*ΔIlfmax)

println("Lf = ", round(Lf*1e6, digits = 3), " μH")

#____________________________________________________________
# Capacitor Design
Vorms = Vrms*0.95
Vop = Vorms*sqrt(2)

Zl = 3*Vorms*Vorms/Sr

Iorms = Vorms/Zl
Iop = Iorms*sqrt(2)
Ir_d = Vdc/(4*fs*Lf*Iop)
ΔIlfmax = Ir_d*Iop
ΔVcfmax = Vr*Vop

Cf = ΔIlfmax/(8*fs*ΔVcfmax)

println("Cf = ", round(Cf*1e6, digits = 3), " μF\n")

#____________________________________________________________
# Verifying Design

#= Theory
    The design should be valid for the entire range. That is, the if 
    the network voltage is at 0.95 p.u. then the ripple voltage should 
    be still be below the specified value. In other words, at 0.95 p.u. 
    the ripple voltage is maximised, which is why we designed the Cf at 
    0.95 p.u. For the ripple current we use the inverse methodology. The 
    ripple current is maximised at 1.05 p.u. 
=#

Vorms = Vrms*1.05
Vop = Vorms*sqrt(2)

Zl = 3*Vorms*Vorms/Sr

Iorms = Vorms/Zl
Iop = Iorms*sqrt(2)
ΔIlfmax = Ir*Iop

Ir_d = Vdc/(4*fs*Lf*Iop)

println("Ir_d = ", round(Ir_d*100, digits = 3), " % @ Vpoc = ", Vorms/Vrms, " p.u.")

if round(Ir_d, digits = 3) <= round(Ir, digits = 3)
    println("Lf has been correctly designed.\n")
else
    println("Lf has been incorrectly designed.\n")
end

Vorms = Vrms*0.95
Vop = Vorms*sqrt(2)

Zl = 3*Vorms*Vorms/Sr

Iorms = Vorms/Zl
Iop = Iorms*sqrt(2)
ΔIlfmax = Ir*Iop

Ir_d = Vdc/(4*fs*Lf*Iop)
Vr_d = Ir_d*Iop/(8*fs*Cf*Vop)

println("Vr_d = ", round(Vr_d*100, digits = 3), " % @ Vpoc = ", Vorms/Vrms, " p.u.")
if round(Vr_d, digits = 3) <= round(Vr, digits = 3)
    println("Cf has been correctly designed.\n")
else
    println("Cf has been incorrectly designed.\n")
end

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
