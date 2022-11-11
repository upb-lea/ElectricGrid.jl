using JuMP
import Ipopt


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
L = 0.00025
R = 0.208
B1 = -omega * L /(R^2 + omega^2 * L^2)
G1 = R /(R^2 + omega^2 * L^2)

G = [[G1, -G1],
     [-G1, G1]]

B = [[B1, -B1],
     [-B1, B1]]

model = Model(Ipopt.Optimizer)

# Variables including var_constraints
@variable(model, 0 <= PG_1 <= Pmax_1)
@variable(model, 0 <= QG_1 <= Qmax_1)
@variable(model, 0 <= v2 <= 400)
@variable(model, -2 <= theta2 <= 2)

# non-linear objectives
@NLobjective(model, Min, PG_1 + QG_1 + v2 + theta2)


@NLconstraint(model, P_Bus1,
    v1 * v1 * (G[1][1] * cos(theta1 - theta1) + B[1][1] * sin(theta1 - theta1)) + 
    v1 * v2 * (G[1][2] * cos(theta1 - theta2) + B[1][2] * sin(theta1 - theta2))  
    == PG_1 - PL_1)

@NLconstraint(model, P_Bus2,
    v2 * v1 * (G[2][1] * cos(theta2 - theta1) + B[2][1] * sin(theta2 - theta1)) + 
    v2 * v2 * (G[2][2] * cos(theta2 - theta2) + B[2][2] * sin(theta2 - theta2)) 
    == PG_2 - PL_2)

@NLconstraint(model, Q_Bus1,
    v1 * v1 * (G[1][1] * sin(theta1 - theta1) - B[1][1] * cos(theta1 - theta1)) + 
    v1 * v2 * (G[1][2] * sin(theta1 - theta2) - B[1][2] * cos(theta1 - theta2)) 
    == QG_1 - QL_1)

@NLconstraint(model, Q_Bus2,
    v2 * v1 * (G[2][1] * sin(theta2 - theta1) - B[2][1] * cos(theta2 - theta1)) + 
    v2 * v2 * (G[2][2] * sin(theta2 - theta2) - B[2][2] * cos(theta2 - theta2))  
    == QG_2 - QL_2)

# Linear constraints:
#@constraint(model, lc1, PG_1 + QG_1 <= Smax_1)

optimize!(model)
println("""
termination_status = $(termination_status(model))
primal_status      = $(primal_status(model))
objective_value    = $(objective_value(model))
""")
println("PG1, QG1, v2, theta2")
println("$(value(PG_1)), $(value(QG_1)), $(value(v2)), $(value(theta2))")

