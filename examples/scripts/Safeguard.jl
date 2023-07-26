using ElectricGrid
using LinearAlgebra
using Polyhedra

using Ipopt
using JuMP

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
    v_len = length(v)
    M_blk = zeros(v_len*N, N)

    for n in 1:N
        for i in 1:v_len*N
            if n == i
                for k in 1:v_len
                    M_blk[k+count-1, n] = v[k]
                end
            end
        end
        count+=v_len
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

A_N2, B_N2, Omega_x2, Omega_u2, W_x_cal2, W_u_cal2 = get_condensed_matrices(A_P10, B_P10, omega_x, omega_u, W_x, W_u, 2)

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

volume(poly_3d)  # TODO: Why is this 0? Problem -> when abort algorithm?

poly_3d_proj = project(poly, [1,2,3])
volume(poly_3d_proj)

poly_2d_v_u_prj = project(poly, [2,3])
plot(poly_2d_v_u_prj)

feas_bA = poly_3D_Hrep.A


#=
for i in 1:N

    A_N, B_N, Omega_x, Omega_u, W_x_cal, W_u_cal = get_condensed_matrices(A_P10, B_P10, omega_x, omega_u, W_x, W_u, i)
    poly_vol = get_poly(A_N, B_N, Omega_x, Omega_u, W_x_cal, W_u_cal, m, n)
    println("")
    println(i)
    println(volume(poly_vol))
    println("")
end
=#

# From Python - validated with quadprog in matlab
i = [0.7, 0.1, 0.5, -0.6, -0.4, -0.6, -0.7, 0.7]
v = [-0.5, -0.4, -0.2, 0.55, 0.3, -0.6, -0.5, 0.55]
a_rl = [0.95, 0.8, 0.9, -0.8, -0.6, -0.95, -0.8, 0.9]
a_sg = [-0.619, -0.123, -0.021, 0.644, 0.128, -0.19, -0.493, 0.275]

feas_b = feas_bA[:,1]
feas_A = feas_bA[:,2:end]

action_rl = a_rl[1]
state = [i[1]; v[1]]
safeguard = Model(Ipopt.Optimizer)

@variable(safeguard, action)

@objective(safeguard, Min, (action - action_rl)^2)

A_x = feas_A[:,1:2]
A_u = feas_A[:,3]
A_x *state
@constraint(safeguard, A_u * action <= feas_b - A_x*state)

print(safeguard)

optimize!(safeguard)


p1 = plot(poly_2d_v_u, label="feas_Set")
scatter!(p1, (v[1], value(action)), mc=:green, ms=10, ma=1, label="u_SG")
plot!(p1, legend=:bottomright)
plot!(p1, [(v[1], action_rl), (v[1], value(action))], arrow = arrow(:closed, 0.1), color = :black)
scatter!(p1, (v[1], action_rl), mc=:red, ms=10, ma=1, label="u_RL")

p = plot(poly_2d_v_u, label="feas_Set")
p_i = plot(poly_2d_i_u)

a_sg_jump = []

for ii in 1:length(v)
    action_rl = a_rl[ii]
    state = [i[ii]; v[ii]]
    safeguard = Model(Ipopt.Optimizer)

    @variable(safeguard, action)

    @objective(safeguard, Min, (action - action_rl)^2)

    A_x = feas_A[:,1:2]
    A_u = feas_A[:,3]
    A_x *state
    @constraint(safeguard, A_u * action <= feas_b - A_x*state)

    print(safeguard)

    optimize!(safeguard)


    scatter!(p, (v[ii], value(action)), mc=:green, ms=10, ma=1, label="u_SG")
    plot!([(v[ii], action_rl), (v[ii], value(action))], arrow = arrow(:closed, 0.1), color = :black)
    scatter!(p, (v[ii], action_rl), mc=:red, ms=10, ma=1, label="u_RL")
    push!(a_sg_jump, value(action))

    scatter!(p_i, (i[ii], value(action)), mc=:green, ms=10, ma=1, label="u_SG")
    plot!([(i[ii], action_rl), (i[ii], value(action))], arrow = arrow(:closed, 0.1), color = :black)
    scatter!(p_i, (i[ii], action_rl), mc=:red, ms=10, ma=1, label="u_RL")
    push!(a_sg_jump, value(action))

end
xlabel!(p, "v / v_lim")
ylabel!(p, "action")

xlabel!(p_i, "i / i_lim")
ylabel!(p_i, "action")

display(p)
display(p_i)
