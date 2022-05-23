using DrWatson
@quickactivate "MicroGridSimWithRL"

using LinearAlgebra

include(srcdir("nodeconstructor.jl"))

"""
investigate the occuring resonance frequencies for a given configuration
"""

CM = [ 0.  1.
        -1.  0.]

a = 34

parameters = Dict()
parameters["source"] = [Dict("fltr" => "LCL", "R" => 0.4, "L1" => 2.3e-3, "L2" => 2.3e-3,"C" => 10e-6)]
parameters["cable"] = [Dict("R" => 0.722, "L" => 0.955e-3, "C" => 8e-09)]
parameters["load"] = [Dict("impedance" => "R", "R" => 14)]

Grid_FC = NodeConstructor(num_source=1, num_loads=1, CM=CM, parameters=parameters)

#draw_graph(Grid_FC)   ---   not yet implemented

A, B, C, D = get_sys(Grid_FC)

evs = eigvals(A)

evs_re_list = []
evs_im_list = []

for ev in evs
    push!(evs_re_list, real(ev))
    push!(evs_im_list, imag(ev)/(2*pi))
end

sort!(evs_im_list)

println(evs)
println(" ")
println(" ")
println(evs_im_list)
println(" ")
println(" ")
println(count(i->(i == 0), evs_re_list))
println(" ")
println(" ")