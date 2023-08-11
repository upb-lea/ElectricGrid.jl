using ElectricGrid
using LinearAlgebra
using ControlSystemsBase
using BenchmarkTools
using DifferentialEquations
using PlotlyJS

function Get_A_x(AAA)
    (rows, columns) = size(AAA)
    AA = Matrix{Any}(undef, (rows, columns))
    for row in 1:rows
        for column in 1:columns
            h = AAA[row, column]
            if isa(h, Number)
                AA[row, column] = x -> h
            else
                AA[row, column] = h
            end
        end
    end

    A(x) = (|>).(x, AA)
    return x -> A(x)
end

CM = [0 1
    -1 0]

# start: Setting up parameters
# Source
R = 1.1e-3
L = 70e-6
C = 250e-6

# Cable
C_b = 1e-4 / 2
L_b = x -> 1e-4
R_b = 1e-3

# Load
R_l = 100
C_l = 1e-2
L_l = 1e-2;

parameters = Dict()

grid_properties = Dict()
grid_properties["fs"] = 10e3
grid_properties["v_rms"] = 230
grid_properties["phase"] = 1;
parameters["grid"] = grid_properties

source1 = Dict()
source_list = []

source1["fltr"] = "L"
source1["R1"] = R
source1["L1"] = L

push!(source_list, source1)
parameters["source"] = source_list

cable = Dict()
cable["R"] = R_b
cable["L"] = L_b
cable["C"] = C_b
cable_list = []

push!(cable_list, cable);
parameters["cable"] = cable_list

load1 = Dict()
load_list = []

load1["impedance"] = "RLC"
load1["R"] = R_l;
load1["L"] = L_l;
load1["C"] = C_l;

push!(load_list, load1);
parameters["load"] = load_list;

# end: Setting up parameters

function bench_nodes_nonlinear(A_nonlinear,A,B,C,D)
    (row, col) = size(A)
    x0 = zeros(row, 1)
    u = [230 for i = 1:n]
    tspan = (0, 0.1)
    function f(dx, x, p, t)
        dx .= A_nonlinear(x) * x + B * p
    end
    prob = ODEProblem(f, x0, tspan, u)
    alg = SSPRK22()
    alg = Euler()
    Δt = 0.01
    sol = solve(prob, alg, dt=Δt)

    return sol.u
end

function bench_nodes_linear(A_nonlinear,A,B,C,D)
    (row, col) = size(A)
    x0 = zeros(row, 1)
    A_linear = A_nonlinear(x0)

    A′ = Float64.(A_linear)
    ts = 0.01
    Ad = exp(A′ * ts)
    Bd = A′ \ (Ad - C) * B
    sys_d = StateSpace(Ad, Bd, C, D, ts)
    ns = length(A[1, :]) # get num of states
    ni = length(B[1, :]) # get num of inputs
    t = collect(0:ts:0.1)
    x0 = [0.0 for i = 1:ns]
    u = [230.0 for i = 1:length(t)]
    uu = [u for i = 1:ni]
    uuu = mapreduce(permutedims, vcat, uu)
    xout, _, _, _ = lsim(sys_d, uuu, t, x0=x0)
    return xout
end

times = []
n = 0
range = 5:5:20
for i = range
    n = i
    source_list = []
    load_list = []
    for j = 1:n
        push!(source_list, source1)
        push!(load_list, load1)
    end
    parameters["source"] = source_list
    parameters["load"] = load_list
    nc = NodeConstructor(num_sources=n, num_loads=n, parameters=parameters)
    global A, B, C, D = GetSystem(nc)

    global A_x = Get_A_x(A)

    b_time = median(@benchmark bench_nodes_linear(A_x,A,B,C,D)).time
    safe = [b_time]
    @show b_time

    result_linear = bench_nodes_linear(A_x,A,B,C,D)

    b_time = median(@benchmark bench_nodes_nonlinear(A_x,A,B,C,D)).time
    @show b_time

    result_nonlinear = bench_nodes_nonlinear(A_x,A,B,C,D)
    append!(safe, b_time)
    append!(times, [safe])

    # @show max_error = norm(result_nonlinear-result_linear,Inf)
end