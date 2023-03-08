using DrWatson
@quickactivate ("JEG")

using Plots

include(srcdir("power_system_theory.jl"))

Timestep = 100 #time step in μs ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_final = 0.4 #time in seconds, total simulation run time

ts = Timestep*1e-6
t = 0:ts:t_final # time
fsys = 50

N = length(t)
Vac = Array{Float64, 2}(undef, 3, N)
Iac = Array{Float64, 2}(undef, 3, N)

Vac[1, :] = 230*sqrt(2)*sin.(2π*fsys.*t)
Vac[2, :] = 230*sqrt(2)*sin.(2π*fsys.*t .- 120*π/180)
Vac[3, :] = 230*sqrt(2)*sin.(2π*fsys.*t .+ 120*π/180)

δ = 45*π/180
Iac[1, :] = 1*sqrt(2)*sin.(2π*fsys.*t .+ δ)
Iac[2, :] = 1*sqrt(2)*sin.(2π*fsys.*t .- 120*π/180 .+ δ)
Iac[3, :] = 1*sqrt(2)*sin.(2π*fsys.*t .+ 120*π/180 .+ δ)

pq0 = Array{Float64, 2}(undef, 3, N)
P = Array{Float64, 1}(undef, N)

for i in 1:N

    pq0[:, i] = pqTheory(Vac[:, i], Iac[:, i])

end

Vout = Plots.plot(t, Vac[1, :], xlabel = "time [s]", ylabel = "V")
Iout = Plots.plot(t, Iac[1, :], xlabel = "time [s]", ylabel = "A")

P = Vac[1, :].*Iac[1, :] .+ Vac[2, :].*Iac[2, :] .+ Vac[3, :].*Iac[3, :]
pout = Plots.plot(t, pq0[1, :], 
title = "Real Power", xlabel = "time [s]", ylabel = "W", label = "thinky")
pout = Plots.plot!(t, P, label = "simple")

qout = Plots.plot(t, pq0[2, :], 
title = "Imaginary Power", xlabel = "time [s]", ylabel = "VAi", label = "thinky")

display(pout)
display(qout)