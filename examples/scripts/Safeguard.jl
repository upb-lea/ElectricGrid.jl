using ElectricGrid
using LinearAlgebra
using Polyhedra

############################################################################################
# New functions

"""
    my_blkdiag_from_vec(v, N)

# Description
Creates a blockdiagonal matrix with N-times vector v on the diagonal

# Arguments

# Keyword Arguments
- `v::Vector`: Vector that appears on diagonal.
- `N::Int`: number of diagonal elements (v).

# Return Values
- `M_blk::Matrix`: Blockdiagonal Matrix with vector v on diagonal

"""
function my_blkdiag_from_vec(v, N)

    count = 1
    v_lec = length(v)
    M_blk = zeros(v_lec*N, N)

    for n in 1:N+1
        for i in 1:size(v)[1]+1
            if n == i
                M_blk[count, i] = v[1]
                M_blk[count+v_lec-1, i] = v[v_lec]
            end
        end
        count+=v_lec
    end

    return M_blk
end

"""
    my_blkdiag_matrix(m, N)

# Description
Creates a blockdiagonal matrix with N-times matrix m on the diagonal

# Arguments

# Keyword Arguments
- `m::Matrix`: Vector that appears on diagonal.
- `N::Int`: number of diagonal elements (v).

# Return Values
- `M_blk::Matrix`: Blockdiagonal Matrix with vector v on diagonal
"""
function my_blkdiag_matrix(m, N)

    m_lines, m_rows = size(m)
    M_blk = zeros(m_lines*(N+1), m_rows*(N+1))
    for n in 1:N+1
        for line in 1:(m_lines)
            for row in 1:(m_rows)
                M_blk[line+m_lines*(n-1), row+m_rows*(n-1)] = m[line, row]
            end
        end
    end
    return M_blk
end

"""
    get_condensed_matrices(A, B, omega_x, omega_u, W_x, W_u, N)

# Description
Calulates the condesed state space and contraints matrices

# Arguments

# Keyword Arguments
- `A::Matrix`: .
- `B::Matrix`: .

# Return Values
-
"""
function get_condensed_matrices(A, B, omega_x, omega_u, W_x, W_u, N)

    m = 1
    n = size(B)[1]

    A_N = 1* Matrix(I, size(A)[1], size(A)[1])
    B_N = zeros(n, m*N)

    Omega_x = omega_x
    Omega_u = []


    W_x_cal = my_blkdiag_matrix(W_x, N)
    W_u_cal = my_blkdiag_from_vec(W_u, N)

    for ii  = 1:N

        A_N = [A_N; A_P10^ii]
        Omega_u = [Omega_u; omega_u]
        Omega_x = [Omega_x; omega_x]

        B_N_newline = A^(ii-1)*B


        for jj in 2:ii

            B_N_newline = hcat(B_N_newline, A^(ii-jj)*B)


        end
        B_N_newline = hcat(B_N_newline, zeros(n, m*(N-ii)))
        B_N = [B_N; B_N_newline]

    end

    return A_N, B_N, Omega_x, Omega_u, W_x_cal, W_u_cal
end

function get_poly(A_N, B_N, Omega_x, Omega_u, W_x_cal, W_u_cal, m, n)

    E = [-W_x_cal*A_N;
    zeros(length(Omega_u), n)]

    G = [W_x_cal * B_N;
        W_u_cal]

    e = [Omega_x;
        Omega_u]

    A_poly = [hcat(W_x,zeros(size(W_x)[1], size(G)[2]));
        hcat(-E, G)]
    B_poly = Float64[omega_x; e]

    hrep_poly = hrep(A_poly, B_poly, BitSet())
    poly = polyhedron(hrep_poly)

    return poly
end

############################################################################################

#=
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

A = env.A[1:2,1:2]
B = env.B[1:2]
=#
# Use matlab vals to compare
A = [-115.705536  -14284.4644;
     3333.04169 -41.6630212]
B = [14285.7143; 0 ]

A_P10 = exp(A * 1e-4)#env.ts)
B_P10 = A \ (A_P10 - I) * B #

i_dead = 40 #env.nc.parameters["source"][1]["i_limit"]
v_dead = 60 #env.nc.parameters["source"][1]["v_limit"]
v_dc = 50 #env.nc.parameters["source"][1]["vdc"]
v_ref = 40
v_lim_poly = 1.2 * v_ref
i_lim_poly = 30

# Work on normalized state space; #TODO: same if not normalized?
A_P10[1,2] = A_P10[1,2]*v_dead/i_dead;
A_P10[2,1] = A_P10[2,1]*i_dead/v_dead;
B_P10[1,1] = B_P10[1,1]/i_dead*v_dc;
B_P10[2,1] = B_P10[2,1]*v_dc/v_dead;


W_u = [-1; 1]
omega_u = [1; 1]
W_x = [-1 0; 1 0; 0 -1; 0 1]
omega_x = [
    i_lim_poly/i_dead # env.nc.parameters["source"][1]["i_limit"];
    i_lim_poly/i_dead;
    v_lim_poly/v_dead;
    v_lim_poly/v_dead
]

N = 3  # maximum number of iteration

m = 1 # number of inputs to the system,
n = size(B)[1]

A_N2, B_N2, Omega_x2, Omega_u2, W_x_cal2, W_u_cal2 = get_condensed_matrices(A_P10, B_P10, omega_x, omega_u, W_x, W_u, N)

poly = get_poly(A_N2, B_N2, Omega_x2, Omega_u2, W_x_cal2, W_u_cal2, m, n)


volume(poly)

poly_2d = eliminate(poly, [3, 4, 5])
volume(poly_2d)

using Plots
plot(poly_2d)

poly_2d_i_u = eliminate(poly, [2, 4, 5])
plot(poly_2d_i_u)

poly_2d_v_u = eliminate(poly, [1, 4, 5])
plot(poly_2d_v_u)

# geting matrices for safeguard
poly_3d = eliminate(poly, [4, 5])
poly_3D_Hrep = doubledescription(poly_3d.vrep)
feas_bA = poly_3D_Hrep.A



for i in 1:N

    A_N, B_N, Omega_x, Omega_u, W_x_cal, W_u_cal = get_condensed_matrices(A_P10, B_P10, omega_x, omega_u, W_x, W_u, i)
    poly_vol = get_poly(A_N, B_N, Omega_x, Omega_u, W_x_cal, W_u_cal, m, n)
    println("")
    println(n)
    println(volume(poly_vol))
    println("")
end
