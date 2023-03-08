using DrWatson
@quickactivate "MicroGridSimWithRL"

using PyCall
@pyinclude(srcdir("nodeconstructor.py"))

using JSON
using Plots
using BenchmarkTools
using Statistics
using DifferentialEquations
using CUDA

nodes = 4
t_end = 0.03

ts=1e-4

parameter=Dict()
parameter["R_source"] = 0.4
parameter["L_source"] = 2.3e-3
parameter["C_source"] = 10e-6
parameter["L_cable"] = 2.3e-3
parameter["R_cable"] = 0.4
parameter["R_load"] = 14
parameter["V_dc"] = 300

results = Dict()


CM_list = JSON.parsefile(srcdir("CM_matrices", "CM_nodes" * string(nodes) * ".json"))
CM = reduce(hcat, CM_list[1])'
CM = convert(Matrix{Int}, CM)
nc = py"NodeConstructor"(nodes, nodes, parameter, CM=CM)

A, B, C, D = nc.get_sys()

A = CuArray(Float32.(A))
B = CuArray(Float32.(B))
t_end = Float32(t_end)
t_start = Float32(0.0)

ns = length(A[1,:])
na = length(B[1,:])

function f!(du, u, p, t)
    du[:] = A * u + B * p
end

p = [230.0 for i = 1:na]
tspan = (t_start,t_end)
u0 = [0.0 for i = 1:ns]

p = CuArray(Float32.(p))
u0 = CuArray(Float32.(u0))

problem = ODEProblem(f!,u0,tspan,p)

tol = Float32(1e-6)

res = solve(problem, BS3(), reltol=tol, abstol=tol, saveat=ts)
#res = solve(problem, BS3(), adaptive=false, dt=ts, dtmin=ts, dtmax=ts)