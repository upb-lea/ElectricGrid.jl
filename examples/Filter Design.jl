Ir = 15/100 # i_rip in %
Vr = 1.54/100 # v_rip in %

Sr = 50e3 # Srated, apparent power in VA
fs = 15e3 # Hz, switching frequency
Vdc = 800 # V, dc voltage

#____________________________________________________________
# Inductor Design
Vorms = 230*1.05
Vop = Vorms*sqrt(2)

Zl = 3*Vorms*Vorms/Sr

Iorms = Vorms/Zl
Iop = Iorms*sqrt(2)

ΔIlfmax = Ir*Iop

Lf = Vdc/(4*fs*ΔIlfmax)

println("Lf = ", Lf*1e6, " μH")

#____________________________________________________________
# Capacitor Design
Vorms = 230*0.95
Vop = Vorms*sqrt(2)

Zl = 3*Vorms*Vorms/Sr

Iorms = Vorms/Zl
Iop = Iorms*sqrt(2)
Ir = Vdc/(4*fs*Lf*Iop)
ΔIlfmax = Ir*Iop
ΔVcfmax = Vr*Vop

Cf = ΔIlfmax/(8*fs*ΔVcfmax)

println("Cf = ", Cf*1e6, " μF")

#____________________________________________________________
# Verifying Design
Vorms = 230*1.05
Vop = Vorms*sqrt(2)

Zl = 3*Vorms*Vorms/Sr

Iorms = Vorms/Zl
Iop = Iorms*sqrt(2)
ΔIlfmax = Ir*Iop

Ir = Vdc/(4*fs*Lf*Iop)
Vr = Ir*Iop/(8*fs*Cf*Vop)

println("Ir = ", Ir*100, " % @ ", Vorms/230, " p.u.")
println("Vr = ", Vr*100, " % @ ", Vorms/230, " p.u.")