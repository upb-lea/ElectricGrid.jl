using DrWatson
@quickactivate "MicroGridSimWithRL"

using PyCall
@pyinclude(srcdir("nodeconstructor.py"))
@pyinclude(srcdir("nodeconstructorcable.py"))

using ControlSystems
using JSON
using Plots
using LinearAlgebra

printit = true
discrete = false
cable = true
num_cm = 1
num_mat_start = 2
num_mat_end = 30


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

global evI_list = []
global evR_list = []
global all_ev = Dict()

for n=num_mat_start:num_mat_end
    if printit
        println(" ")
        println("--------- n=$n ---------")
        println(" ")
    end

    global results[n] = Dict()

    if cable
        global num_cm = 1
    else
        CM_list = JSON.parsefile(srcdir("CM_matrices", "CM_nodes" * string(n) * ".json"))
    end

    for i=1:num_cm
        if cable
            nc = py"NodeConstructorCable"(n, n)
        else
            CM = reduce(hcat, CM_list[i])'
            CM = convert(Matrix{Int}, CM)

            nc = py"NodeConstructor"(n, n, parameter, CM=CM)
        end

        global A, B, C, D = nc.get_sys()
        if discrete A = exp(A*ts) end

        ns = length(A[1,:])
        na = length(B[1,:])

        notnull = count(i->(i!= 0), A)
        null = count(i->(i== 0), A)
        global ev = eigvals(A)

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
            evR_ge_zero = 0
            evR_gr = -999999.0
            evI_gr = -999999.0
            all_ev[n] = Dict()
            all_ev[n]["re"] = []
            all_ev[n]["im"] = []
            
            for e in ev
                if real(e) >= 0.0
                    evR_ge_zero += 1
                end
                if abs(imag(e)) >= evI_gr
                    evR_gr = real(e)
                    evI_gr = abs(imag(e))
                end
                append!(all_ev[n]["re"], abs(real(e)))
                append!(all_ev[n]["im"], abs(imag(e))/(2*pi))
            end
        end

        if printit
            println("A: $(ns)x$(ns)")
            println("Anzahl Nullen: $(null),   Anzahl Nicht-Null: $(notnull)")
            #println("Ratio Nicht-Null/Null: $(round(notnull/null, digits=4))")
            println("Anteil Nullen: $(round((null/(null + notnull))*100, digits=2))%")
            if discrete
                println("Eigenwerte abs >= 1: $(evg1)")
                println("Größter abs EW: $(ev_gr)")
            else
                println("Eigenwerte mit Realteil >= 0: $(evR_ge_zero)")
                println("Größter EW-Realteil: $(evR_gr)")
                println("Größter absoluter EW-Imaginärteil: $(evI_gr)")
            end
            println(" ")
        end

        if !discrete
            append!(evI_list, evI_gr/(2*pi))
            append!(evR_list, evR_gr)
        end
        global results[n][i] = round((null/(null + notnull))*100, digits=2)
    end
end

arraytoplot = [results[n][1] for n=num_mat_start:num_mat_end]

#plot(arraytoplot)
evI_list = Float64.(evI_list)
#histogram(evI_list)
plot(num_mat_start:num_mat_end, evI_list)
#ylims!((2000,6000))
#plot(num_mat_start:num_mat_end, evR_list)