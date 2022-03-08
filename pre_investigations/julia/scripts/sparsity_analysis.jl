using DrWatson
@quickactivate "MicroGridSimWithRL"

using PyCall
@pyinclude(srcdir("nodeconstructor.py"))

using ControlSystems
using JSON
using Plots
using LinearAlgebra

printit = true
discrete = true
num_cm = 30
num_mat_start = 8
num_mat_end = 8

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

        A, B, C, D = nc.get_sys()
        if discrete A = exp(A*ts) end

        ns = length(A[1,:])
        na = length(B[1,:])

        notnull = count(i->(i!= 0), A)
        null = count(i->(i== 0), A)
        ev = eigvals(A)

        if discrete
            evg1 = 0
            ev_gr = -999999.0
            for e in ev
                if abs(e) >= 1
                    evg1 += 1
                end
                if abs(e) >= ev_gr
                    ev_gr = abs(e)
                end
            end
        else
            evs_ge_zero = 0
            ev_gr = -999999.0
            
            for e in ev
                if real(e) >= 0.0
                    evs_ge_zero += 1
                end
                if real(e) >= ev_gr
                    ev_gr = real(e)
                end
            end
        end

        if printit
            println("A: $(ns)x$(ns)")
            println("Anzahl Nullen: $(null),   Anzahl Nicht-Null: $(notnull)")
            println("Ratio Nicht-Null/Null: $(round(notnull/null, digits=4))")
            println("Anteil Nullen: $(round((1-(notnull/null))*100, digits=2))%")
            if discrete
                println("Eigenwerte abs >= 1: $(evg1)")
                println("Größter abs EW: $(ev_gr)")
            else
                println("Eigenwerte >= 0: $(evs_ge_zero)")
                println("Größter EW: $(ev_gr)")
            end
            println(" ")
        end

        global results[n][i] = round((1-(notnull/null))*100, digits=2)
    end
end

arraytoplot = [results[n][1] for n=num_mat_start:num_mat_end]

plot(arraytoplot)