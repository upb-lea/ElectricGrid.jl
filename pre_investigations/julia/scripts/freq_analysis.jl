using DrWatson
@quickactivate "MicroGridSimWithRL"
include(srcdir("nodeconstructor.jl"))

"""
File for Frequ-Analysis using noddeconstructor.jl
    Different number of nodes (=2 -> 2 sources and 2 loads) can be investigated from
    num_mat_start to num_mat_end. Thereby the CM Matrices are loaded from files to have 
    reproduceable experiements. 
    In the End 3 Plots are printed:
     - The maximal accuring frequency per node and the corresponding real-part
     - all accuring frequencys per node (Scatter-plot)
     - A step response of the (currently) secound state (if LC-filter -> v_LC over C)
       therefore the last (!) model/state-space system is used at the moment
"""

using ControlSystems
using JSON
using Plots
using LinearAlgebra

discrete = false
julia = true
cut_outliers = false
num_cm = 1
num_mat_start = 1 
num_mat_end = 30

ts=1e-6

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

    CM_list = JSON.parsefile(srcdir("CM_matrices", "CM_nodes" * string(n) * ".json"))

    for i=1:num_cm
        CM = reduce(hcat, CM_list[i])'
        CM = convert(Matrix{Int}, CM)

        global nc = NodeConstructor(num_source=n, num_loads=n, CM=CM)
        global parameters = nc.parameters
        
        # for debug perpose print used params
        #println(parameters)
        #println(" ")

        global A, B, C, D = get_sys(nc)

        #println("")
        #println(A)
        #println("")

        if discrete A = exp(A*ts) end

        ns = length(A[1,:])
        na = length(B[1,:])

        notnull = count(i->(i!= 0), A)
        null = count(i->(i== 0), A)
        global ev = eigvals(A)
        global eigenvec = eigvecs(A)

        #println(" ")
        #println(ev)
        #println(" ")
        #println(eigenvec)
        #println(" ")
        
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
                println("Größter absoluter EW-Imaginärteil: $(evI_gr/2/pi)")
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
p1 = plot(num_mat_start:num_mat_end, evI_list,  ylabel="Im{EW_max_imag}/ (2* π)")
#ylims!((2000,6000))
#plot(num_mat_start:num_mat_end, evR_list)

evR_list = Float64.(evR_list)
p2 = plot(num_mat_start:num_mat_end, evR_list, xlabel="Nodes", ylabel="Re{EW_max_imag}")
display(plot(p1,p2,layout=(2,1)))

#print(all_ev)
#println(" ")
#println(" ")
#println(all_ev[1]["im"])
#print(all_ev[2]["im"])

x_ax = []
y_ax = []
for n=num_mat_start:num_mat_end
    global x_ax = vcat(x_ax, repeat([n], length(all_ev[n]["im"])))
    global y_ax = vcat(y_ax, all_ev[n]["im"])

end


#x = length(all_ev[3]["im"])
#x_ax = repeat([3], x)

#print(x_ax)
#print(y_ax)

p3 = scatter(x_ax, y_ax, xlabel="Nodes", ylabel="Im{EWs(A)} / (2* π)")

display(plot(p3))

####### plot step response
# use therefore the last (!) model
Ad = exp(A*ts)
Bd = A \ (Ad - C) * B

global sys_d = StateSpace(Ad, Bd, C, D, ts)

t = collect(0:ts:0.005)

ns = length(A[1,:])
na = length(B[1,:])

#x0 = real(eigenvec[:,2])

#global x0 = [0.0 for i = 1:length(t)]

global x0 = [0.0 for i = 1:ns]
#global u = rand(Float64, ( length(t) )) .*2 .-1
global u = [250.0 for i = 1:length(t)]
global uu = [u for i = 1:na ]
global uuu = mapreduce(permutedims, vcat, uu)
global ttt = t

#xout = lsim(sys_d, uuu, ttt, x0)
xout, _, _, _ = lsim(sys_d,uuu,ttt,x0=x0)

p4 = plot(t, xout[2,:], xlabel="time", ylabel="v_LC / V")
display(plot(p4))

p5 = plot(t, xout[ns,:], xlabel="time", ylabel="v_cable / V")
display(plot(p5))