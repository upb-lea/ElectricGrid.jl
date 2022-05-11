using DrWatson
@quickactivate "MicroGridSimWithRL"

using PyCall
@pyinclude(srcdir("nodeconstructor.py"))
@pyinclude(srcdir("nodeconstructorcable.py"))
@pyinclude(srcdir("nodeconstructorcableloads.py"))
include(srcdir("nodeconstructor.jl"))

using ControlSystems
using JSON
using Plots
using LinearAlgebra

printit = true
discrete = false
cable = true
cableloads = true
julia = true
cut_outliers = false
num_cm = 1
num_mat_start = 2  # hier auf 1 und dann num_LC = 0 ?!
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

global nc = nothing

for n=num_mat_start:num_mat_end
    if printit
        println(" ")
        println("--------- n=$n ---------")
        println(" ")
    end

    global results[n] = Dict()

    if cable
        CM_list = JSON.parsefile(srcdir("CM_matrices", "CM_nodes" * string(n) * ".json"))
        #global num_cm = 1
    else
        CM_list = JSON.parsefile(srcdir("CM_matrices", "CM_nodes" * string(n) * ".json"))
    end

    for i=1:num_cm
        if julia
            CM = reduce(hcat, CM_list[i])'
            CM = convert(Matrix{Int}, CM)

            global nc = NodeConstructor(num_source=n, num_load=n, CM=CM)
            global parameter = nc.parameter
        elseif cableloads
            global nc = py"NodeConstructorCableLoads"(n, n, CM=CM)
            global parameter = nc.parameter
        elseif cable
            CM = reduce(hcat, CM_list[i])'
            CM = convert(Matrix{Int}, CM)

            global nc = py"NodeConstructorCable"(n, n, CM=CM)
            global parameter = nc.parameter
        else
            CM = reduce(hcat, CM_list[i])'
            CM = convert(Matrix{Int}, CM)

            global nc = py"NodeConstructor"(n, n, parameter, CM=CM)
        end

        if julia
            global A, B, C, D = get_sys(nc)
        else
            global A, B, C, D = nc.get_sys()
        end

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
                    println(" EW gr 0 RE: $(e)")
                end
                if abs(imag(e)) >= evI_gr
                    if cut_outliers
                        if abs(imag(e)) <= 1000000
                            evR_gr = real(e)
                            evI_gr = abs(imag(e))
                        end
                    else
                        evR_gr = real(e)
                        evI_gr = abs(imag(e))
                    end
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
p1 = plot(num_mat_start:num_mat_end, evI_list,  ylabel="Im{EW_max_imag}")
#ylims!((2000,6000))
#plot(num_mat_start:num_mat_end, evR_list)

evR_list = Float64.(evR_list)
p2 = plot(num_mat_start:num_mat_end, evR_list, xlabel="Nodes", ylabel="Re{EW_max_imag}")

display(plot(p1,p2,layout=(2,1)))

print(all_ev)
println(" ")
println(" ")
println(all_ev[2]["im"])
#print(all_ev[2]["im"])

x_ax = []
y_ax = []
for n=num_mat_start:num_mat_end
    global x_ax = vcat(x_ax, repeat([n], length(all_ev[n]["im"])))
    global y_ax = vcat(y_ax, all_ev[n]["im"])

end


#x = length(all_ev[3]["im"])
#x_ax = repeat([3], x)

print(x_ax)
print(y_ax)

p3 = scatter(x_ax, y_ax, xlabel="Nodes", ylabel="Im{EWs(A)}")

display(plot(p3))