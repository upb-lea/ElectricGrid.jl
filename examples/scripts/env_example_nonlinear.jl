using ElectricGrid
using DifferentialEquations
using PlotlyJS

CM = [0 1
    -1 0]

# Source
R = 1.1e-3
L = 70e-6
C = 250e-6

# Cable
C_b = 1e-4/2
L_b = 1e-4
R_b = 1e-3

# Load
R_l = 100
C_l = 1e-2
L_l = 1e-2;
parameters = Dict()

grid_properties = Dict()
grid_properties["fs"] =  10e3
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

# @show parameters
env1 = ElectricGridEnv(num_sources=1, num_loads=1, CM = CM, parameters = parameters, verbosity = 2)
value = 5
sol = []

for i = 1:1000
    env1([230])
    append!(sol,env1.x[value])
end

parameters = Dict()

grid_properties = Dict()
grid_properties["fs"] =  10e3
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
cable["L"] = x->L_b
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

# @show parameters
env2 = ElectricGridEnv(num_sources=1, num_loads=1, CM = CM, parameters = parameters, verbosity = 2)

t_t = collect(env1.t0:env1.ts:env1.t)
p1 = scatter(x=t_t,y=sol,mode="lines",name="Nonlinear")
sol = []

for i = 1:1000
    env2([230])
    append!(sol,env2.x[value])
end

t_t = collect(env2.t0:env2.ts:env2.t)
p2 = scatter(x=t_t,y=sol,mode="lines",name="Linear")

plot([p1,p2])