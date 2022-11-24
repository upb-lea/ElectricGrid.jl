using JuMP
import Ipopt

model = Model(Ipopt.Optimizer)

zero_expression = @NLexpression(model, 0.0)

# Constant values
f = 50
omega = 2π*f

# for every Source: v is fixed 230
# for one Source: theta is fixed 0

CM = [  0   0   0   1
        0   0   0   2
        0   0   0   3
        -1  -2  -3  0  ]

num_source = 3
num_load = 1
num_nodes = num_source + num_load

@variable(model, nodes[1 : num_source + num_load, ["v", "theta", "P", "Q"]])
fix(nodes[1, "theta"], 0.0)
fix(nodes[1, "v"], 230.0)

function set_bounds(variable, start_value, low_bound, up_bound)
    if !is_fixed(variable)
        set_lower_bound(variable, low_bound)
        set_upper_bound(variable, up_bound)
        set_start_value(variable, start_value)
    end
end

for i = 1:num_source+num_load
    if i <= num_source
        fix(nodes[i, "v"], 230.0)

        set_bounds(nodes[i, "theta"], 0.0, -pi/2, pi/2)
        set_bounds(nodes[i, "P"], (num_load * 1000.0) / num_source, -1000.0, 1000.0)
        set_bounds(nodes[i, "Q"], (num_load * -100.0) / num_source, -1000.0, 1000.0)
    else
        fix(nodes[i, "P"], -1000.0)
        fix(nodes[i, "Q"], -100.0)
    end
end

function get_degree(CM = CM)
    result = zeros(Int, size(CM)[1])

    for i=1:size(CM)[1]
        result[i] = count(x -> x != 0, CM[i,:])
    end

    result
end

function get_cable_connections(CM = CM)

    result = Vector{Vector{Int64}}(undef, size(CM)[1])

    for i=1:size(CM)[1]
        result[i] = filter(x -> x != 0, abs.(CM[i,:]))
    end

    return result
end

@variable(model, cables[1 : maximum(CM), ["L", "X_R", "C_L"]])
cable_conductance = Array{NonlinearExpression, 1}(undef, maximum(CM))
cable_susceptance_0 = Array{NonlinearExpression, 1}(undef, maximum(CM))
cable_susceptance_1 = Array{NonlinearExpression, 1}(undef, maximum(CM))

for i=1:maximum(CM)
    set_bounds(cables[i, "L"], 0.00025, 0.0002, 0.0003)
    set_bounds(cables[i, "X_R"], 0.38, 0.1, 0.5)
    set_bounds(cables[i, "C_L"], 0.0016, 0.0001, 0.01)

    #R = (omega*L)/X_R
    #C = C_L*L
    cable_conductance[i] = @NLexpression(model, ((omega * cables[i, "L"]) / cables[i, "X_R"]) / (((omega*cables[i, "L"]) / cables[i, "X_R"])^2 + omega^2 * cables[i, "L"]^2))
    cable_susceptance_1[i] = @NLexpression(model, (-omega * cables[i, "L"] / (((omega*cables[i, "L"])/cables[i, "X_R"])^2 + omega^2 * cables[i, "L"]^2)))
    cable_susceptance_0[i] = @NLexpression(model, (-omega * cables[i, "L"] / (((omega*cables[i, "L"])/cables[i, "X_R"])^2 + omega^2 * cables[i, "L"]^2)) + omega*cables[i, "C_L"]*cables[i, "L"]/2)
end

cable_connections = get_cable_connections()
G = Array{NonlinearExpression, 2}(undef, num_nodes, num_nodes) # should be symmetric
B = Array{NonlinearExpression, 2}(undef, num_nodes, num_nodes) # should be symmetric

for i in 1:num_nodes

    # diagonal terms
    G[i, i] = @NLexpression(model, sum( cable_conductance[cable_connections[i]][j] for j in 1:length(cable_connections[i])))
    B[i, i] = @NLexpression(model, sum( cable_susceptance_0[cable_connections[i]][j] for j in 1:length(cable_connections[i])))

    # off diagonal terms
    for k in (i+1):num_nodes # this is over the upper triangle

        if CM[i, k] != 0

            cable_num = abs(CM[i, k])

            G[i, k] = @NLexpression(model, -1*cable_conductance[cable_num])
            B[i, k] = @NLexpression(model, -1*cable_susceptance_1[cable_num])
            G[k, i] = G[i, k]
            B[k, i] = B[i, k]
            
        else

            G[i, k] = zero_expression
            B[i, k] = zero_expression
            G[k, i] = zero_expression
            B[k, i] = zero_expression
        end
    end
end

#= Goal:

    At the end the optimisation should solve for L, X_R, and C_L.
    Constraints (i.e. minimum and maximum) should be enforced on these values

=#
#= 
    another constraint to add - the active power, P, should be less than the Surge Impedance Loading, SIL)

    SIL = 2*v^2*sqrt(cables[i, "C_L"])

    nodes[i, "P"] < 2*v^2*sqrt(cables[i, "C_L"])
=#

# non-linear objectives
#@NLobjective(model, Min, abs(sum(nodes[1:num_source,"P"]))/1000 + abs(sum(nodes[1:num_source,"Q"]))/1000 + sum(nodes[num_source+1:end,"v"])/230 + abs(sum(nodes[2:end,"theta"]))/(2π))

function get_node_connections(CM = CM)

    result = Vector{Vector{Int64}}(undef, size(CM)[1])

    for i=1:size(CM)[1]
        result[i] = findall(x -> x != 0, CM[i,:])
    end

    return result
end

for i in 1:num_nodes

    @NLconstraint(model, P_Bus1,

        #= for j in num_nodes
            if CM[i,j] != 0

            end
        end =#

        nodes[1,"v"] * nodes[1,"v"] * (G[1,1] * cos(theta1 - theta1) + B[1,1] * sin(theta1 - theta1)) + 
        nodes[1,"v"]  * nodes[4,"v"]  * (G[1,4] * cos(theta1 - theta2) + B[1,4] * sin(theta1 - theta2))  
        == nodes[1,"P"])

    @NLconstraint(model, P_Bus2,
        v2 * v1 * (G[2,1] * cos(theta2 - theta1) + B[2,1] * sin(theta2 - theta1)) + 
        v2 * v2 * (G[2,2] * cos(theta2 - theta2) + B[2,2] * sin(theta2 - theta2)) 
        == P)

end
    

# Linear constraints:
#@constraint(model, lc1, PG_1 + QG_1 <= Smax_1)

optimize!(model)
println("""
termination_status = $(termination_status(model))
primal_status      = $(primal_status(model))
objective_value    = $(objective_value(model))
""")
println("PG1, QG1, v2, theta2, L, X_R, C_L")
println("$(value(PG_1)), $(value(QG_1)), $(value(v2)), $(value(theta2)), $(value(L)), $(value(X_R)), $(value(C_L))")

