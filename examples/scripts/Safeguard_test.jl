using ElectricGrid
using LinearAlgebra

CM = [ 0. 0. 1.
        0. 0. 2.
        -1. -2. 0.]

#R_load, L_load, _, _ = ParallelLoadImpedance(50e3, 0.95, 230)

parameters =
Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}(
                        "pwr" => 200e3,
                        "fltr" => "LCL",
                        "L1" => 70e-6,
                        "R1" => 1.1e-3,
                        "C" => 300e-6,
                        "R_C" => 7e-3,
                        "control_type" => "classic",
                        "mode" => "Droop",),
                    ],
    "load"   => Any[
        Dict{Any, Any}(
            "impedance" => "R",
            "R" => 100,
            "v_limit" => 1e4,
            "i_limit" => 1e4)
        ],
    "grid" => Dict{Any, Any}(
        "phase" => 1,
        "ramp_end" => 0.04,)
)



env = ElectricGridEnv(
    #CM =  CM,
    parameters = parameters,
    t_end = 1,
    action_delay = 0,
    verbosity = 0)


env.A

#= To exclude the cable C from matrix and conpare result with matlab:
 - LCL filter, only use LC (not 2nd L) for SG calculation, since we have no cable in matlab
    - Matrix look to be the same -> A_P10
 - Caculate condensed system A_N and B_N
    - Same?
 - Calcuate G and E,...
    - Same?
 - Calcuate Polyhedron
    - Same?
 - Use it to ensure safety
=#

A_P10 = env.A[1:2,1:2]
B_P10 = env.B[1:2]

W_u = [-1; 1]
omega_u = [
    env.nc.parameters["source"][1]["vdc"];
    env.nc.parameters["source"][1]["vdc"]
         ]


W_x = [
    -1 0;
    1 0;
    0 -1
    0 1
    ]

omega_x = [
    env.nc.parameters["source"][1]["i_limit"];
    env.nc.parameters["source"][1]["i_limit"];
    env.nc.parameters["source"][1]["v_limit"];
    env.nc.parameters["source"][1]["v_limit"]
]

A_sqare = A_P10^2

N = 3  # maximum number of iteration

m = 1 # number of inputs to the system,
n = size(B_P10)[1]
#TODO: get from B -> here has (2,) dim and not (2,1), why?
# better: (n, m) = size(B_P10)

global A_N = 1* Matrix(I, size(A_P10)[1], size(A_P10)[1])
global B_N = zeros(n, m*N)

using SparseArrays

W_x_cal = sparse(W_x)
W_x_sp = sparse(W_x)
W_u_cal = sparse(W_u)

W_x_cal = blockdiag(W_x_cal, W_x_cal)
W_u_cal = cat(W_u,  reverse(W_u); dims=(1,2))

for ii  = 1:N

    global A_N = [A_N; A_P10^ii]

    global W_x_cal = blockdiag(W_x_cal, W_x_sp)
    #Wglobal W_u_cal = blockdiag(W_u_cal, W_u_cal)



    B_N_newline = A_P10^(ii-1)*B_P10


    for jj in 2:ii

        B_N_newline = hcat(B_N_newline, A_P10^(ii-jj)*B_P10)


    end
    B_N_newline = hcat(B_N_newline, zeros(n, m*(N-ii)))




    global B_N = [B_N; B_N_newline]



end

println("")
println(A_N)
println("")
println(B_N)
println("")
println("")
println(W_x_cal)
println("")
println(W_u_cal)
println("")
