import json
import timeit
from threadpoolctl import threadpool_limits
import numpy as np
import cupy as cp
import os

numpydot = True
discrete = True
threads = 8
use_cuda = False

#os.environ["MKL_DYNAMIC"] = "FALSE"
#os.environ["OMP_DYNAMIC"] = "FALSE"
#os.environ["MKL_NUM_THREADS"] = str(threads)
#os.environ["OMP_NUM_THREADS"] = str(threads)
#os.environ["MKL_DOMAIN_NUM_THREADS"] = "MKL_BLAS=4"

f = open('./test_variables.json', "r")
test_variables = json.loads(f.read())
if use_cuda:
    A = cp.array(test_variables['A_d'])
    B = cp.array(test_variables['B_d'])
    C = cp.array(test_variables['C'])
    D = cp.array(test_variables['D'])
    u_fix = cp.array(test_variables['u_fix'])
    Ad = cp.array(test_variables['A_d_foh'])
    Bd1 = cp.array(test_variables['B1_d_foh'])
    Bd0 = cp.array(test_variables['B0_d_foh'])
else:
    A = np.array(test_variables['A_d'])
    B = np.array(test_variables['B_d'])
    C = np.array(test_variables['C'])
    D = np.array(test_variables['D'])
    u_fix = np.array(test_variables['u_fix'])
    Ad = np.array(test_variables['A_d_foh'])
    Bd1 = np.array(test_variables['B1_d_foh'])
    Bd0 = np.array(test_variables['B0_d_foh'])

n_states = A.shape[0]
n_inputs = B.shape[1]
n_outputs = C.shape[0]
n_steps = len(u_fix[0, :])

if use_cuda:
    xout = cp.zeros((n_states, n_steps))
    yout = cp.zeros((n_outputs, n_steps))
else:
    xout = np.zeros((n_states, n_steps))
    yout = np.zeros((n_outputs, n_steps))


limit = 10

def sim():
    for n in range(0, limit):
        # print(str(100*n/limit) + '%')
        if discrete:
            if use_cuda:
                for i in range(0, n_steps - 1):
                    global yout
                    xout[:, i + 1] = cp.dot(A, xout[:, i]) + cp.dot(B, u_fix[:, i])
                    yout[:, i] = cp.dot(C, xout[:, i]) + cp.dot(D, u_fix[:, i])
            elif numpydot:
                for i in range(0, n_steps - 1):
                    xout[:, i + 1] = np.dot(A, xout[:, i]) + np.dot(B, u_fix[:, i])
                    yout[:, i] = np.dot(C, xout[:, i]) + np.dot(D, u_fix[:, i])
            else:
                for i in range(0, len(u_fix[0, :]) - 1):
                    xout[:, i + 1] = A @ xout[:, i] + B @ u_fix[:, i]
                    yout[:, i] = C @ xout[:, i] + D @ u_fix[:, i]
        else:
            if use_cuda:
                for i in range(1, n_steps):
                    xout[:, i] = cp.dot(Ad, xout[:, i - 1]) + cp.dot(Bd0, u_fix[:, i - 1]) + cp.dot(Bd1, u_fix[:, i])
                yout = cp.dot(C, xout) + cp.dot(D, u_fix)
            elif numpydot:
                for i in range(1, n_steps):
                    xout[:, i] = np.dot(Ad, xout[:, i - 1]) + np.dot(Bd0, u_fix[:, i - 1]) + np.dot(Bd1, u_fix[:, i])
                yout = np.dot(C, xout) + np.dot(D, u_fix)
            else:
                for i in range(1, n_steps):
                    xout[:, i] = (Ad @ xout[:, i - 1] + Bd0 @ u_fix[:, i - 1] + Bd1 @ u_fix[:, i])
                yout = C @ xout + D @ u_fix
    if use_cuda:
        cp.cuda.Device().synchronize()


with threadpool_limits(limits=threads):
    res = timeit.timeit(lambda: sim(), number=5)
    #sim()

print(res)
