import json
from os import makedirs

import numpy as np
from pre_investigations.python.dare.utils.nodeconstructor import NodeConstructor
from pre_investigations.python.dare.utils.nodeconstructorcable import NodeConstructorCable


"""
asd
"""

makedirs('CM_matrices', exist_ok=True)

# define parameters
R = 0.4
L = 2.3e-3
C = 10e-6
LT = 2.3e-3
R_load = 14

parameter = dict()
parameter['R_source'] = R
parameter['L_source'] = L
parameter['C_source'] = C
parameter['L_cable'] = LT
parameter['R_cable'] = R
parameter['R_load'] = R_load

# set up CM
#An entry in the CM reads as follows: Object x (row) is connected to object y (column) via connection line z:
# x ---z---> y
CM = np.array([[0, 0, 1],
               [0, 0, 2],
               [-1, -2, 0]])

num_nodes = [2, 4, 6, 8, 10]
loops = 1
for nodes in num_nodes:
    CM_array = []
    parameter_array = []
    for loop in range(loops):
        # grid with 2 sources, 1 load and a total of 2 connections
        power_grid = NodeConstructorCable(nodes, nodes,  CM = None)
        #power_grid.draw_graph()

        #print(power_grid.CM)

        CM_array.append(power_grid.CM.tolist())
        parameter_array.append(power_grid.parameter)

    with open('CM_matrices_lcl/CM_nodes'+str(nodes)+'.json', 'w') as \
            outfile: json.dump(CM_array, outfile)

    with open('CM_matrices_lcl/parameter_nodes'+str(nodes)+'.json', 'w') as \
            outfile: json.dump(parameter_array, outfile)


