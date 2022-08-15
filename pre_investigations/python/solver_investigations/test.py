import json
import numpy as np
import scipy
import copy
import control

from pyinstrument import Profiler

from pre_investigations.python.dare.utils.nodeconstructor import NodeConstructor
import pre_investigations.python.solver_investigations.custom_control as cc


num_nodes = 30
numcm = 0
bit32 = False
offline_expm = True
discrete = False
profiling = False
create_tv = True
ts = 1e-4

parameter = dict()
parameter['R_source'] = 0.4
parameter['L_source'] = 2.3e-3
parameter['C_source'] = 10e-6
parameter['L_cable'] = 2.3e-3
parameter['R_cable'] = 0.4
parameter['R_load'] = 14
parameter['V_dc'] = 300

t = np.arange(0, 0.03 + ts, ts)

f = open('pre_investigations/python/solver_investigations/CM_matrices/CM_nodes' + str(num_nodes) + '.json', "r")
CM_list = json.loads(f.read())
CM = np.array(CM_list[numcm])

power_grid = NodeConstructor(num_nodes, num_nodes, parameter, CM=CM)

A_sys, B_sys, C_sys, D_sys = power_grid.get_sys()

A_d = scipy.linalg.expm(A_sys * ts)
A_inv = scipy.linalg.inv(A_sys)
B_d = A_inv @ (A_d - np.eye(A_sys.shape[0])) @ B_sys
C_d = copy.copy(C_sys)

if discrete:
    dt = ts
else:
    dt = 0

sys = cc.ss(A_d, B_d, C_d, 0, dt=dt, bit32=bit32, offline_expm=offline_expm)

u_fix = np.array([230] * power_grid.num_source)[:, None] * np.ones(
                                (power_grid.num_source, len(t)))
x0 = np.zeros((A_sys.shape[0],))
if bit32:
    u_fix = u_fix.astype(np.float32)
    x0 = x0.astype(np.float32)

if profiling:
    profiler = Profiler(interval=0.0001)
    profiler.start()

T, yout_d, xout_d = cc.forced_response(sys, T=t, U=u_fix, X0=x0, return_x=True, squeeze=True)
#T, yout_d, xout_d = control.forced_response(sys, T=t, U=u_fix, X0=x0, return_x=True, squeeze=True)

if profiling:
    profiler.stop()
    #profiler.print()
    profiler.open_in_browser()

if create_tv:
    test_variables = {}
    test_variables['A_d'] = A_d.tolist()
    test_variables['A_sys'] = A_sys.tolist()
    test_variables['B_sys'] = B_sys.tolist()
    test_variables['B_d'] = B_d.tolist()
    test_variables['C'] = C_d.tolist()
    test_variables['D'] = np.zeros((C_d.shape[0], B_d.shape[1])).tolist()
    test_variables['u_fix'] = u_fix.tolist()
    test_variables['A_d_foh'] = sys.Ad.tolist()
    test_variables['B1_d_foh'] = sys.Bd1.tolist()
    test_variables['B0_d_foh'] = sys.Bd0.tolist()
    with open('pre_investigations/python/solver_investigations/test_variables.json', 'w') as tv_outfile: json.dump(test_variables, tv_outfile)