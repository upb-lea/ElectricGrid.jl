# f(x) = x^2
# g(x) = x

# A(x) = [f(x),g(x)]

# A(2)
function NonlinearInductance(value, drop, length)
    return x -> value * (1 - drop / 100 * 1 / (1 + exp(-(abs(x) - length / 2) / 2)))
end


using ElectricGrid
using ControlSystemsBase
using PlotlyJS

CM = [0 1
    -1 0]

# Source
R = 1.1e-3
L = 70e-6
C = 250e-6

# Cable
C_b = 1e-4 / 2
value = 1e-4
drop = 20
len = 5000
# L_b = 1e-4
L_b = NonlinearInductance(value, drop, len)
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

# source1["fltr"] = "LC"
source1["fltr"] = "L"
source1["R1"] = R
source1["L1"] = L
# source1["C"] = C

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

nc = NodeConstructor(num_sources=1, num_loads=1, parameters=parameters, CM=CM);
A,B,C,D = GetSystem(nc)