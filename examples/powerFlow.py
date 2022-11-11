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

v1 = n.Const(value=230)  # Output voltage first Inverter/Volt
PG_2 = n.Const(value=0)  # Power at Bus 2/W
QG_2 = n.Const(value=0)  # Power at Bus 2/W

PL_1 = n.Const(value=0)  # Input power first inverter/W
PL_2 = n.Const(value=1000)  # Input power second Inverter7W
QL_1 = n.Const(value=0)  # Input power first inverter/W
QL_2 = n.Const(value=500)  # Input power second Inverter7W

theta1 = n.Const(value=0)  # Phase angle at Bus 1/rad

L = 0.00025
R = 0.208
B1 = -omega * L /(R**2 + omega**2 * L**2)
G1 = R /(R**2 + omega**2 * L**2)

B = np.array([[B1, -B1],
              [-B1, B1]])

G = np.array([[G1, -G1],
              [-G1, G1]])

# Power flow equations
n.Equation(v1 * v1 * (G[0][0] * n.cos(theta1 - theta1) + B[0][0] * n.sin(theta1 - theta1)) + \
           v1 * v2 * (G[0][1] * n.cos(theta1 - theta2) + B[0][1] * n.sin(theta1 - theta2))  == PG_1 - PL_1)
n.Equation(v2 * v1 * (G[1][0] * n.cos(theta2 - theta1) + B[1][0] * n.sin(theta2 - theta1)) + \
           v2 * v2 * (G[1][1] * n.cos(theta2 - theta2) + B[1][1] * n.sin(theta2 - theta2)) == PG_2 - PL_2)


n.Equation(v1 * v1 * (G[0][0] * n.sin(theta1 - theta1) - B[0][0] * n.cos(theta1 - theta1)) + \
           v1 * v2 * (G[0][1] * n.sin(theta1 - theta2) - B[0][1] * n.cos(theta1 - theta2)) == QG_1 - QL_1)
n.Equation(v2 * v1 * (G[1][0] * n.sin(theta2 - theta1) - B[1][0] * n.cos(theta2 - theta1)) + \
           v2 * v2 * (G[1][1] * n.sin(theta2 - theta2) - B[1][1] * n.cos(theta2 - theta2))  == QG_2 - QL_2)



#n.options.IMODE = 7


# Solve simulation
n.solve()
print("PG1, QG1, v2, theta2, theat3, QG3")
print(n._variables)
#n.cleanup()


