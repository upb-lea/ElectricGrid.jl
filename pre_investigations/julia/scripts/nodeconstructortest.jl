using DrWatson
@quickactivate "MicroGridSimWithRL"

using ControlSystems
using JSON
using Plots
using LinearAlgebra
using PyCall
@pyinclude(srcdir("nodeconstructor.py"))
@pyinclude(srcdir("nodeconstructorcable.py"))

include(srcdir("custom_control.jl"))

ts = 1e-5
t_end = 0.003
cable = false

p = Dict()
p["R_source"] = 0.4
p["L_source"] = 2.3e-3
p["C_source"] = 10e-6
p["L_cable"] = 2.3e-3
p["R_cable"] = 0.4
p["R_load"] = 14
p["V_dc"] = 300

CM = [0 0 1
    0 0 2
    -1 -2 0]

if cable
    n = py"NodeConstructorCable"(2, 1, CM=CM)
else
    n = py"NodeConstructor"(2, 1, CM=CM, p)
end

A, B, C, D = n.get_sys()

ns = length(A[1,:])
na = length(B[1,:])
t = collect(0:ts:t_end)

Ad = exp(A*ts)
Bd = A \ (Ad - C) * B
sys_d = ss(Ad, Bd, C, D, ts)

x0 = [0.0 for i = 1:ns]
u = [230.0 for i = 1:length(t)]
uu = [u for i = 1:na ]
uuu = mapreduce(permutedims, vcat, uu)

result = custom_lsim(sys_d,uuu,t,x0=x0)
resulttoplot = result[2,:]

display(plot(resulttoplot))

ev = eigvals(A)

evR_ge_zero = 0
evR_gr = -999999.0
evI_gr = -999999.0
all_ev = Dict()
all_ev["re"] = []
all_ev["im"] = []

for e in ev
    if real(e) >= 0.0
        global evR_ge_zero += 1
    end
    if abs(imag(e)) >= evI_gr
        global evR_gr = real(e)
        global evI_gr = abs(imag(e))
    end
    append!(all_ev["re"], abs(real(e)))
    append!(all_ev["im"], abs(imag(e))/(2*pi))
end