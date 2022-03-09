using DrWatson
@quickactivate "MicroGridSimWithRL"

using PyCall
@pyinclude(srcdir("nodeconstructor.py"))

using ControlSystems
using JSON
using Plots
using LinearAlgebra

printit = false
discrete = false
num_cm = 1
num_mat_start = 110
num_mat_end = num_mat_start

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

for n=num_mat_start:num_mat_end
    if printit
        println(" ")
        println("--------- n=$n ---------")
        println(" ")
    end

    global results[n] = Dict()

    CM_list = JSON.parsefile(srcdir("CM_matrices", "CM_nodes" * string(n) * ".json"))

    for i=1:num_cm
        CM = reduce(hcat, CM_list[i])'
        CM = convert(Matrix{Int}, CM)

        nc = py"NodeConstructor"(n, n, parameter, CM=CM)

        global A, B, C, D = nc.get_sys()
        if discrete A = exp(A*ts) end

        #ns = length(A[1,:])
        #na = length(B[1,:])

        #notnull = count(i->(i!= 0), A)
        #null = count(i->(i== 0), A)
        #ev = eigvals(A)
    end
end
