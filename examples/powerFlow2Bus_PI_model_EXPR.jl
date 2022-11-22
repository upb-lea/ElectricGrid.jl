using JuMP
import Ipopt

model = Model(Ipopt.Optimizer)

# Constant values
f = 50
omega = 2*π*f

# Reference Bus
v1 = 230
theta1 = 0
PL_1 = 0
QL_1 = 0
Pmax_1 = 5000
Qmax_1 = 5000

# Load Bus
PL_2 = 1000
PG_2 = 0
QL_2 = 500
QG_2 = 0

# Concuctance and Susceptance matrices -> get from nc.get_Ybus() -> real()/imag()
#L = 0.00025
#X_R = 0.38 # ratio of omega*L/R
#C_L = 0.0016 # ratio of C/L
@variable(model, 0.0002 <= L <= 0.0003)
@variable(model, 0.1 <= X_R <= 0.5)
@variable(model, 0.0001 <= C_L <= 0.01)

set_start_value(L, 0.00025)
set_start_value(X_R, 0.38)
set_start_value(C_L, 0.0016)

#C = C_L*L = 0.4e-6 
#R = (omega*L)/X_R = 0.208 
#B1 = -omega * L /(((omega*L)/X_R)^2 + omega^2 * L^2)
#G1 = ((omega*L)/X_R) /(((omega*L)/X_R)^2 + omega^2 * L^2)
#B0 = (-omega * L /(((omega*L)/X_R)^2 + omega^2 * L^2)) + omega*C_L*L
@NLexpression(model, B1, -omega * L /(((omega*L)/X_R)^2 + omega^2 * L^2))
@NLexpression(model, G1, ((omega*L)/X_R) /(((omega*L)/X_R)^2 + omega^2 * L^2))
@NLexpression(model, B0, (-omega * L /(((omega*L)/X_R)^2 + omega^2 * L^2)) + omega*C_L*L)

#= Goal:

    At the end the optimisation should solve for L, X_R, and C_L.
    Constraints (i.e. minimum and maximum) should be enforced on these values

=#

#G = [[G1, -G1],
#     [-G1, G1]]

#B = [[B0, -B1],
#     [-B1, B0]]

@NLexpression(model, G[1:2, 1:2], 0)
@NLexpression(model, B[1:2, 1:2], 0)

G[1,1] = @NLexpression(model, G1)
G[1,2] = @NLexpression(model, -G1)
G[2,1] = @NLexpression(model, -G1)
G[2,2] = @NLexpression(model, G1)
B[1,1] = @NLexpression(model, B0)
B[1,2] = @NLexpression(model, -B1)
B[2,1] = @NLexpression(model, -B1)
B[2,2] = @NLexpression(model, B0)


# Variables including var_constraints
@variable(model, 0 <= PG_1 <= Pmax_1)
@variable(model, 0 <= QG_1 <= Qmax_1)
@variable(model, 0 <= v2 <= 400)
@variable(model, -2 <= theta2 <= 2)

set_start_value(PG_1, 1004.9)
set_start_value(QG_1, 488.6)
set_start_value(v2, 230.0)
set_start_value(theta2, 0)

# non-linear objectives
@NLobjective(model, Min, abs(PG_1)/200 + abs(QG_1)/500 + abs(v2)/240 + abs(theta2)/0.1)


@NLconstraint(model, P_Bus1,
    v1 * v1 * (G[1,1] * cos(theta1 - theta1) + B[1,1] * sin(theta1 - theta1)) + 
    v1 * v2 * (G[1,2] * cos(theta1 - theta2) + B[1,2] * sin(theta1 - theta2))  
    == PG_1 - PL_1)

@NLconstraint(model, P_Bus2,
    v2 * v1 * (G[2,1] * cos(theta2 - theta1) + B[2,1] * sin(theta2 - theta1)) + 
    v2 * v2 * (G[2,2] * cos(theta2 - theta2) + B[2,2] * sin(theta2 - theta2)) 
    == PG_2 - PL_2)

@NLconstraint(model, Q_Bus1,
    v1 * v1 * (G[1,1] * sin(theta1 - theta1) - B[1,1] * cos(theta1 - theta1)) + 
    v1 * v2 * (G[1,2] * sin(theta1 - theta2) - B[1,2] * cos(theta1 - theta2)) 
    == QG_1 - QL_1)

@NLconstraint(model, Q_Bus2,
    v2 * v1 * (G[2,1] * sin(theta2 - theta1) - B[2,1] * cos(theta2 - theta1)) + 
    v2 * v2 * (G[2,2] * sin(theta2 - theta2) - B[2,2] * cos(theta2 - theta2))  
    == QG_2 - QL_2)

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

