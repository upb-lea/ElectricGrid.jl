import json
from os import makedirs

import numpy as np
from pre_investigations.python.dare.utils.nodeconstructor import NodeConstructor


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
parameter['L_cabel'] = LT
parameter['R_cabel'] = R
parameter['R_load'] = R_load

# set up CM
#An entry in the CM reads as follows: Object x (row) is connected to object y (column) via connection line z:
# x ---z---> y
CM = np.array([[0, 0, 1],
               [0, 0, 2],
               [-1, -2, 0]])

num_nodes = [2, 4, 6, 8, 10]
loops = 30
for nodes in num_nodes:
    CM_array = []
    for loop in range(loops):
        # grid with 2 sources, 1 load and a total of 2 connections
        power_grid = NodeConstructor(nodes, nodes, parameter=parameter, CM = None)
        #power_grid.draw_graph()

        #print(power_grid.CM)

        CM_array.append(power_grid.CM.tolist())

    with open('CM_matrices/CM_nodes'+str(nodes)+'.json', 'w') as \
            outfile: json.dump(CM_array, outfile)
    """
    f = open('CM_matrices/CM_nodes'+str(nodes)+'.json', "r")
    # Reading from file
    CM_load = json.loads(f.read())

    for CM_load_use in CM_array:
        power_grid_2 = NodeConstructor(nodes, nodes, parameter=parameter, CM=np.array(CM_load_use))
        power_grid_2.draw_graph()
    """

