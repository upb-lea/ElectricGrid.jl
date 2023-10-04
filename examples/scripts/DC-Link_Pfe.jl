using ElectricGrid
using JuMP
using Ipopt
#using ReinforcementLearning

num_sources = 2
num_loads = 1

CM = [0.0 1.0 2.0
    -1.0 0.0 3.0
    -2.0 -3.0 0.0]

S_source = 1e4

S_load = 22e3
pf_load = 0.9
v_rms = 230
R_load, L_load, X, Z = ParallelLoadImpedance(S_load, pf_load, v_rms)



parameters = Dict{Any,Any}(
    "source" => Any[
        Dict{Any,Any}("pwr" => S_source, "control_type" => "classic", "mode" => "Swing", "fltr" => "LC", "i_limit" => 1e4, "v_limit" => 1e4,),
    ],
    "load" => Any[
        Dict{Any,Any}("impedance" => "RLC", "R" => R_load, "v_limit" => 1e4, "i_limit" => 1e4)
    ],
    #"cable"   => Any[
    #                Dict{Any, Any}("R" => 1e-3, "L" => 1e-4, "C" => 1e-4),
    #                ],
    "grid" => Dict{Any,Any}("fs" => 1e4, "phase" => 3, "v_rms" => 230, "f_grid" => 50, "ramp_end" => 0.00)
)

cable_list = []
cable = Dict()
cable["R"] = 1e-7#0.00001    # Ω, line resistance #0.208
cable["L"] = 1e-7#0.25e-5 # H, line inductance
cable["C"] = 3e-7#5e-12  # F, line capacitance
cable["i_limit"] = 10e7   # A, line current limit

push!(cable_list, cable, cable, cable)
parameters["cable"] = cable_list

# S2_L1 = NodeConstructor(num_sources=10, num_loads=10)#, CM=CM, parameters=parameters)
S2_L1 = NodeConstructor(num_sources=num_sources, num_loads=num_loads, CM=CM, parameters=parameters, invoke_pfe=false)

num_nodes = num_sources + num_loads

Y = GetYBus(S2_L1)
G = real(Y)
B = imag(Y)

model = Model(Ipopt.Optimizer)

@variable(model, nodes[1:num_nodes, ["v", "theta", "P", "Q"]])

P1_max = 50e3
Q1_max = 10e3
# Bus 1: Reference Bus -> Swing
fix(nodes[1, "v"], 230)
fix(nodes[1, "theta"], 0)
SetBounds(nodes[1, "P"], 10e3, -P1_max, P1_max)
SetBounds(nodes[1, "Q"], 500, -Q1_max, Q1_max)

# Bus 3: Load Bus, draws P and Q -> PQ --> LOAD
fix(nodes[3, "P"], -20e3)
fix(nodes[3, "Q"], -1e3)
SetBounds(nodes[3, "theta"], 0.0, -0.25 * pi / 2, 0.25 * pi / 2)
SetBounds(nodes[3, "v"], 230, 0.95 * 230, 1.05 * 230)

Q2_max = 50e3
# Bus 2: PV Bus -> PV
fix(nodes[2, "P"], 5e3)
#fix(nodes[2, "Q"], 1e3)
fix(nodes[2, "v"], 230)
SetBounds(nodes[2, "theta"], 0.0, -0.25 * pi / 2, 0.25 * pi / 2)
#SetBounds(nodes[2, "v"], 230, 0.95 * 230, 1.05 * 230)
SetBounds(nodes[2, "Q"], 0, -Q2_max, Q2_max)


node_cons = GetNodeConnections(CM)
P_node = Array{NonlinearConstraintRef,1}(undef, num_nodes)
Q_node = Array{NonlinearConstraintRef,1}(undef, num_nodes)

for i in 1:num_nodes
    P_node[i] = @NLconstraint(model,
        nodes[i, "P"] == (nodes[i, "v"] * sum(
            nodes[j, "v"] * (
                G[i, j] * cos(nodes[i, "theta"] - nodes[j, "theta"])
                + B[i, j] * sin(nodes[i, "theta"] - nodes[j, "theta"])
            ) for j in node_cons[i]))
    )

    Q_node[i] = @NLconstraint(model, nodes[i, "Q"] == nodes[i, "v"] * sum(
        nodes[j, "v"] * (
            G[i, j] * sin(nodes[i, "theta"] - nodes[j, "theta"])
            - B[i, j] * cos(nodes[i, "theta"] - nodes[j, "theta"])
        ) for j in node_cons[i]))
end

set_silent(model)    # turns off output of the solver #TODO put to verbosity?
optimize!(model)

println("""
        termination_status = $(termination_status(model))
        primal_status      = $(primal_status(model))
        objective_value    = $(objective_value(model))
        """)
global optimizer_status = GetOptimizerStatus(model)
@show optimizer_status["termination_status"]
println(value.(nodes))

println("""
        -----------------POWERSYSTEM THEORY-----------------
        """)

parameters2 = Dict{Any,Any}(
    "source" => Any[
        Dict{Any,Any}("pwr" => sqrt(P1_max^2 + Q1_max^2)*3, "control_type" => "classic", "mode" => "Swing", "fltr" => "LC", "i_limit" => 1e4, "v_limit" => 1e4,),
        Dict{Any,Any}("pwr" => 10e3*3, "control_type" => "classic", "mode" => "PV", "fltr" => "LC", "i_limit" => 1e4, "v_limit" => 1e4, "p_set" =>5e3*3),
        Dict{Any,Any}("pwr" => 21e3*3, "control_type" => "classic", "mode" => "PQ", "fltr" => "LC", "i_limit" => 1e4, "v_limit" => 1e4, "p_set" =>-20e3*3, "q_set" => -1e3*3), #TODO: use that instead of "load"; - correct? -> load-flag
    ],
    #"load" => Any[
    #    Dict{Any,Any}("impedance" => "RL", "R" => R_load, "L" => L_load, "v_limit" => 1e4, "i_limit" => 1e4)
    #],   # Works!
    #
    #"cable"   => Any[
    #                Dict{Any, Any}("R" => 1e-3, "L" => 1e-4, "C" => 1e-4),
    #                ],
    "grid" => Dict{Any,Any}("fs" => 1e4, "phase" => 3, "v_rms" => 230, "f_grid" => 50, "ramp_end" => 0.00)
)

S2_L1 = NodeConstructor(num_sources=3, num_loads=0, CM=CM, parameters=parameters2, invoke_pfe=true, verbosity=2);
