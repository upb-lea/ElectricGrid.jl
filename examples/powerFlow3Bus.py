import gekko

from gekko import GEKKO
import math
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
from scipy.optimize import minimize
from scipy.optimize import Bounds

f =50
omega = 2*np.pi*f

n = GEKKO(remote=False)
PG_1 = n.Var(value=500)  # Power at Bus 1/W
QG_1 = n.Var(value=100)  # Power at Bus 1/W
v2 = n.Var(value=200)  # Load voltage/Volt
theta2 = n.Var(value=0)  # Phase angle at Bus 2/rad
theta3 = n.Var(value=0)  # Phase angle at Bus 2/rad
QG_3 = n.Var(value=100)  # Phase angle at Bus 2/rad
#b_var = n.Var(value=0)


v1 = n.Const(value=230)  # Output voltage first Inverter/Volt
theta1 = n.Const(value=0)  # Phase angle at Bus 1/rad
PL_1 = n.Const(value=0)  # Input power first inverter/W
QL_1 = n.Const(value=0)  # Input power first inverter/W

PL_2 = n.Const(value=100)  # Input power second Inverter7W
PG_2 = n.Const(value=0)  # Power at Bus 2/W
QL_2 = n.Const(value=500)  # Input power second Inverter7W
QG_2 = n.Const(value=0)  # Power at Bus 2/W

v3 = n.Const(value=235)  # Output voltage first Inverter/Volt
PG_3 = n.Const(value=100)  # Power at Bus 2/W
PL_3 = n.Const(value=100)  # Input power second Inverter7W
QL_3 = n.Const(value=100)  # Input power first inverter/W


L = 0.0025
#L = 0.00025
R = 0.208
B1 = -omega * L /(R**2 + omega**2 * L**2)

#I = B *V

G1 = R /(R**2 + omega**2 * L**2)

#B = np.array([[B1 + b_var, -b_var, -B1],
#              [-b_var, 2*B1, -B1],
#              [-B1, -B1, 2*B1]])

B = np.array([[2*B1, -B1, -B1],
              [-B1, 2*B1, -B1],
              [-B1, -B1, 2*B1]])

G = np.array([[2*G1, -G1, -G1],
              [-G1, 2*G1, -G1],
              [-G1, -G1, 2*G1]])

# Power flow equations
n.Equation(v1 * v1 * (G[0][0] * n.cos(theta1 - theta1) + B[0][0] * n.sin(theta1 - theta1)) + \
           v1 * v2 * (G[0][1] * n.cos(theta1 - theta2) + B[0][1] * n.sin(theta1 - theta2)) + \
           v1 * v3 * (G[0][2] * n.cos(theta1 - theta3) + B[0][2] * n.sin(theta1 - theta3))
           == PG_1 - PL_1)
n.Equation(v2 * v1 * (G[1][0] * n.cos(theta2 - theta1) + B[1][0] * n.sin(theta2 - theta1)) + \
           v2 * v2 * (G[1][1] * n.cos(theta2 - theta2) + B[1][1] * n.sin(theta2 - theta2)) +\
           v2 * v3 * (G[1][2] * n.cos(theta2 - theta3) + B[1][2] * n.sin(theta2 - theta3)) == PG_2 - PL_2)

n.Equation(v3 * v1 * (G[2][0] * n.cos(theta3 - theta1) + B[2][0] * n.sin(theta3 - theta1)) + \
           v3 * v2 * (G[2][1] * n.cos(theta3 - theta2) + B[2][1] * n.sin(theta3 - theta2)) +\
           v3 * v3 * (G[2][2] * n.cos(theta3 - theta3) + B[2][2] * n.sin(theta3 - theta3)) == PG_3 - PL_3)


n.Equation(v1 * v1 * (G[0][0] * n.sin(theta1 - theta1) - B[0][0] * n.cos(theta1 - theta1)) + \
           v1 * v2 * (G[0][1] * n.sin(theta1 - theta2) - B[0][1] * n.cos(theta1 - theta2))+ \
           v1 * v3 * (G[0][2] * n.sin(theta1 - theta3) - B[0][2] * n.cos(theta1 - theta3))== QG_1 - QL_1)

n.Equation(v2 * v1 * (G[1][0] * n.sin(theta2 - theta1) - B[1][0] * n.cos(theta2 - theta1)) + \
           v2 * v2 * (G[1][1] * n.sin(theta2 - theta2) - B[1][1] * n.cos(theta2 - theta2)) + \
           v2 * v3 * (G[1][2] * n.sin(theta2 - theta3) - B[1][2] * n.cos(theta2 - theta3)) == QG_2 - QL_2)

n.Equation(v3 * v1 * (G[2][0] * n.sin(theta3 - theta1) - B[2][0] * n.cos(theta3 - theta1)) + \
           v3 * v2 * (G[2][1] * n.sin(theta3 - theta2) - B[2][1] * n.cos(theta3 - theta2))+ \
           v3 * v3 * (G[2][2] * n.sin(theta3 - theta3) - B[2][2] * n.cos(theta3 - theta3))== QG_3 - QL_3)



#n.options.IMODE = 7


# Solve simulation
n.solve()
print("PG1, QG1, v2, theta2, theat3, QG3")
print(n._variables)
n.cleanup()


