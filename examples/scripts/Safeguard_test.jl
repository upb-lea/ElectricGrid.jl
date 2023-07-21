using ElectricGrid
using LinearAlgebra
using Polyhedra

parameters =
Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}(
                        "pwr" => 10e3,
                        "fltr" => "LCL",
                        "L1" => 70e-6,
                        "R1" => 1.1e-3,
                        "C" => 300e-6,
                        "R_C" => 7e-3,
                        "control_type" => "classic",
                        "mode" => "Droop",
                        "v_limit" => 60,
                        "i_limit" => 40,
                        "vdc" => 50),
                    ],
    "load"   => Any[
        Dict{Any, Any}(
            "impedance" => "R",
            "R" => 100,
            "v_limit" => 1e2,
            "i_limit" => 1e2)
        ],
    "grid" => Dict{Any, Any}(
        "phase" => 1,
        "ramp_end" => 0.04,)
)



env = ElectricGridEnv(
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
#=
A = env.A[1:2,1:2]
B = env.B[1:2]
=#

# Use matlab vals to compare
A = [-115.705536  -14284.4644;
     3333.04169 -41.6630212]
B = [14285.7143; 0 ]

A_P10 = exp(A * env.ts) #fastExpm(A*ts) might be a better option
    #Bd = A \ (Ad - I) * B #This may be bad for large sizes, maybe QR factorise, then use ldiv!
B_P10 = A \ (A_P10 - I) * B #


i_dead = env.nc.parameters["source"][1]["i_limit"]
v_dead = env.nc.parameters["source"][1]["v_limit"]
v_dc = env.nc.parameters["source"][1]["vdc"]
A_P10[1,2] = A_P10[1,2]*v_dead/i_dead;
A_P10[2,1] = A_P10[2,1]*i_dead/v_dead;
B_P10[1,1] = B_P10[1,1]/i_dead*v_dc;
B_P10[2,1] = B_P10[2,1]*v_dc/v_dead;


W_u = [-1; 1]

omega_u = [1; 1]


W_x = [
    -1 0;
    1 0;
    0 -1
    0 1
    ]

v_ref = 40
v_lim_poly = 1.2 * v_ref
i_lim_poly = 30

#=
omega_x = [
    env.nc.parameters["source"][1]["i_limit"];
    env.nc.parameters["source"][1]["i_limit"];
    env.nc.parameters["source"][1]["v_limit"];
    env.nc.parameters["source"][1]["v_limit"]
]

omega_u = [
    env.nc.parameters["source"][1]["vdc"];
    env.nc.parameters["source"][1]["vdc"]
         ]


=#

omega_x = [
    i_lim_poly/env.nc.parameters["source"][1]["i_limit"];
    i_lim_poly/env.nc.parameters["source"][1]["i_limit"];
    v_lim_poly/env.nc.parameters["source"][1]["v_limit"];
    v_lim_poly/env.nc.parameters["source"][1]["v_limit"]
]

function my_blkdiag_from_vec(v, N)

    count = 1
    v_lec = length(v)
    W_u_cal = zeros(v_lec*N, N)

    for n in 1:N+1
        println(n)
        for i in 1:size(v)[1]+1
            if n == i
                W_u_cal[count, i] = v[1]
                W_u_cal[count+v_lec-1, i] = v[v_lec]
            end
        end
        count+=2
    end

    return W_u_cal
end

function my_blkdiag_matrix(m, N)

    m_lines, m_rows = size(m)
    W_x_cal = zeros(m_lines*(N+1), m_rows*(N+1))
    for n in 1:N+1
        for line in 1:(m_lines)
            for row in 1:(m_rows)
                    W_x_cal[line+m_lines*(n-1), row+m_rows*(n-1)] = m[line, row]
            end
        end
    end
    return W_x_cal
end

A_sqare = A_P10^2

N = 3  # maximum number of iteration

m = 1 # number of inputs to the system,
n = size(B_P10)[1]
#TODO: get from B -> here has (2,) dim and not (2,1), why?
# better: (n, m) = size(B_P10)

global A_N = 1* Matrix(I, size(A_P10)[1], size(A_P10)[1])
global B_N = zeros(n, m*N)

global Omega_x = omega_x
global Omega_u = []

using SparseArrays


W_x_cal = my_blkdiag_matrix(W_x, N)
#W_u_cal = cat(W_u,  reverse(W_u); dims=(1,2))
W_u_cal = my_blkdiag_from_vec(W_u, N)

for ii  = 1:N

    global A_N = [A_N; A_P10^ii]

    global Omega_u = [Omega_u; omega_u]
    global Omega_x = [Omega_x; omega_x]

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

E = [
    -W_x_cal*A_N;
    zeros(length(Omega_u), n)
    ]

G = [
    W_x_cal * B_N;
    W_u_cal
    ]

e = [
    Omega_x;
    Omega_u
    ]

A_poly = [
    hcat(W_x,zeros(size(W_x)[1], size(G)[2]));
    hcat(-E, G)
    ]

B_poly = Float64[omega_x; e]

hrep_poly = hrep(A_poly, B_poly, BitSet())

poly = polyhedron(hrep_poly)

volume(poly)

poly_2d = eliminate(poly, [3, 4, 5])

using Plots
plot(poly_2d)

poly_2d_i_u = eliminate(poly, [2, 4, 5])
plot(poly_2d_i_u)

poly_2d_v_u = eliminate(poly, [1, 4, 5])
plot(poly_2d_v_u)

poly_2d_v_u_hrep = doubledescription(poly_2d_v_u.vrep)
poly_2d_v_u_hrep.A

poly_3d = eliminate(poly, [4, 5])

poly_3D_Hrep = doubledescription(poly_3d.vrep)
poly_3D_Hrep.A
#=
m = Polyhedra.Mesh(poly_3d)
import Makie
Makie.mesh(m, color=:blue)
Makie.wireframe(m)
plot(poly_3d)
=#
using MeshCat
vis = Visualizer()
setobject!(vis, m)
open(vis)
