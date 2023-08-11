using ElectricGrid
using LinearAlgebra
using ControlSystemsBase
using BenchmarkTools

CM = [0 1
    -1 0]

# start: Setting up parameters
# Source
R = 1.1e-3
L = 70e-6
C = 250e-6

LL(i) = L_low_1 + (L_high_1-L_low_1)/2 * (1-2/π*atan(i)) + L_low_2 + (L_high_2-L_low_2)/2*(1-2/π*atan(i))
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

nc = NodeConstructor(num_sources=1, num_loads=1, parameters=parameters, CM=CM);

A, B, C, D = GetSystem(nc);

# (rows, columns) = size(A)
# AA = Matrix{Any}(undef, (rows, columns))
# for row in 1:rows
#     for column in 1:columns
#         h = A[row, column]
#         if isa(h, Number)
#             AA[row, column] = x -> h
#         else
#             AA[row, column] = h
#         end
#     end
# end

# # Here we could try to make it much faster for example
# # when we remember in which line functions are and in which line not. (Pipe operator)
# AAA(x) = (|>).(x, AA)

# x_state_space = [0 0 0 0 0]
# AAAA = AAA(x_state_space)


function get_b!(A)
    b = []
    (col,row) = size(A)
    for i = 1:col
        a = A[i,:]
        if any(x->isa(x,Function)==true,a)
            append!(b,true)
            for j = 1:row
                a_j = A[i,j]
                if isa(a_j,Number)
                    A[i,j] = x->a_j
                end
            end
        else append!(b,false)
        end
    end
    return b
end

b = get_b!(A)
E = Matrix{Int}(I,size(A))
idx_list = findall(x->x==true,b)
# .! negates all entries of b
A_linear = Diagonal(.!b)*A
siz = size(A)
function AA(x)
    # Here I use smart Matrix multiplikation Rules. At least in my eyes
    for i in idx_list
        h = x[i].|>A[i,:]
        result += E[:,i]*h'
    end
    return A_linear+result
end

function benchtestA_x(n)
    x = [0 0 0 0 0]'
    for i = 1:n
        x = AA(x) * x
    end
end
x_state = [0,0,0,0,0]
AAA = AA(x_state);
function benchtestA(n)
    x = [0 0 0 0 0]'
    for i = 1:n
        x = AAA * x
    end
end
n = 1000;
@btime benchtestA_x(n)
@btime benchtestA(n)