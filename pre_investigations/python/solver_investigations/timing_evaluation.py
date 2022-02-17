import copy
import math
import timeit

import control
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import ode, odeint

from pre_investigations.python.dare.utils.nodeconstructor import NodeConstructor


# time measure based on https://note.nkmk.me/en/python-timeit-measure/

def plot_result(results: dict, nodes: list, t_end_vec: list):
    """
    Function to plot the results of time comparison for scaling issue of the power grid for different number of nodes
    and runtimes.
    :param results: Includes measured data:
            results['methode'] -> Used methode (e.g. control.py; scipy;... - more to be defined)
            results['methode_args'] -> Additional arguments of the methode: e.g. control.py -> continuous/discrete
                                                                                     scipy -> solver like LSODA,...
            results['times'] -> list of time results of the experiments in numpy.array(len(nodes),len(t_end_vec))
    :param nodes: Number of grid nodes (doubles, if =5 -> 5 loads and 5 sources) - multiplied by 2 for plotting
    :param t_end_vec: time till simulation is run
    :return:
    """
    for m in range(len(results['methode'])):
        time_result = results['times_mean'][m]

        num_plt_rows = int(math.ceil(len(t_end_vec) / 2))
        fig, ax = plt.subplots(num_plt_rows, 2, figsize=(12, 10))

        col = 0
        row = 0

        jump = False

        for l in range(len(t_end_vec)):

            if l >= round(len(t_end_vec) / 2) and not jump and not num_plt_rows == 1:
                row += 1
                col = 0
                jump = True
            if num_plt_rows == 1:
                ax[col].step(nodes, time_result[:, l])
                # lt.plot(t,result[:steps,0], label = 'i1')
                if col == round(len(t_end_vec) / 2) - 1:
                    ax[col].set_xlabel('nodes')
                if row == 0:
                    ax[col].set_ylabel('$mean(execution-time)\,/\,\mathrm{s}$')
                ax[col].title.set_text(f't_end =  {t_end_vec[l]} s')
                ax[col].grid()
            else:
                ax[col, row].step(nodes, time_result[:, l])
                # lt.plot(t,result[:steps,0], label = 'i1')
                if col == round(len(t_end_vec) / 2) - 1:
                    ax[col, row].set_xlabel('nodes')
                if row == 0:
                    ax[col, row].set_ylabel('$mean(execution-time)\,/\,\mathrm{s}$')
                ax[col, row].title.set_text(f't_end =  {t_end_vec[l]} s')
                ax[col, row].grid()

            col += 1
        fig.suptitle(results['methode'][m] + " using " + results['methode_args'][m])
        plt.show()


def timing_experiment_simulation(repeat: int = 5, loops: int = 10, num_nodes: list[int] = np.arange(2, 12, 2).tolist(),
                                 t_end: list[float] = np.arange(2, 12, 2).tolist(), ts: float = 1e-4,
                                 methode: list[str] = ['control.py', 'control.py', 'scipy_ode', 'scipy_odeint'],
                                 methode_args: list[str] = ['continuous', 'discrete', 'LSODA', 'LSODA'],
                                 parameter: dict = None):
    """

    :param repeat:
    :param loops:
    :param num_nodes:
    :param t_end:
    :param ts:
    :param methode:
    :param methode_args:
    :param parameter:
    :return:
    """
    steps = int(1 / ts)

    if parameter is None:
        parameter = dict()
        parameter['R_source'] = 0.4
        parameter['L_source'] = 2.3e-3
        parameter['C_source'] = 10e-6
        parameter['L_cabel'] = 2.3e-3
        parameter['R_cabel'] = 0.4
        parameter['R_load'] = 14
        parameter['V_dc'] = 300

    result_per_execution = []
    # methode = ['control.py', 'control.py', 'scipy_ode', 'scipy_odeint']
    # methode_args = ['continuous', 'discrete', 'LSODA', 'LSODA']

    # methode = ['scipy_ode', 'scipy_odeint']
    # methode_args = ['LSODA', 'LSODA']

    t_result_mean = np.zeros([len(methode), len(num_nodes), len(t_end)])
    t_result_std = np.zeros([len(methode), len(num_nodes), len(t_end)])

    res_mean_list = []
    res_std_list = []

    def env_model_ode(t, x):  # , arg):
        u_ode = np.array([230] * power_grid.num_source)
        u_random = np.random.uniform(-1, 1, power_grid.num_source) * parameter['V_dc']
        return A_sys @ x + B_sys @ u_ode

    def run_scipy_ode(t_end_ode, ode_solver_):
        while ode_solver_.successful() and ode_solver_.t <= steps * ts * t_end_ode:
            ode_solver_.integrate(ode_solver_.t + ts)

    for n in range(len(methode)):

        for k in range(len(num_nodes)):
            # Define dummy grid
            power_grid = NodeConstructor(num_nodes[k], num_nodes[k], parameter)
            A_sys, B_sys, C_sys, D_sys = power_grid.get_sys()

            if methode in ['control.py'] and methode_args == 'discrete':
                A_d = scipy.linalg.expm(A_sys * ts)
                A_inv = scipy.linalg.inv(A_sys)
                B_d = A_inv @ (A_d - np.eye(A_sys.shape[0])) @ B_sys
                C_d = copy.copy(C_sys)

                sys = control.ss(A_d, B_d, C_d, 0, dt=True)

            else:
                sys = control.ss(A_sys, B_sys, C_sys, D_sys)

            # generate init state
            x0 = np.zeros((A_sys.shape[0],))

            if methode[n] in ['scipy_ode']:
                ode_solver = ode(env_model_ode)
                ode_solver.set_initial_value(x0, 0)

            for l in range(len(t_end)):
                # define time vector
                t = np.arange(0, t_end[l] + ts, ts)

                # simple input signal of constant 230V from all sources
                u_vec = np.array([230]).repeat(power_grid.num_source)[:, None] * np.ones(
                    (power_grid.num_source, len(t)))

                u_random = np.random.uniform(-1, 1, (power_grid.num_source, len(t))) * parameter['V_dc']

                if methode[n] in ['control.py']:
                    res_list = timeit.repeat(
                        lambda: control.forced_response(sys, T=t, U=u_vec, X0=x0, return_x=True, squeeze=True)
                        , repeat=repeat, number=loops)

                if methode[n] in ['scipy_ode']:
                    res_list = timeit.repeat(lambda: run_scipy_ode(t_end[l], ode_solver), repeat=repeat, number=loops)

                if methode[n] in ['scipy_odeint']:
                    res_list = timeit.repeat(lambda: odeint(env_model_ode, x0, t, tfirst=True), repeat=repeat,
                                             number=loops)

                # correct for loops-time execution:
                result_per_loop = ([x / loops for x in res_list])

                # print(result_per_loop)

                t_result_mean[n][k][l] = np.mean(result_per_loop)
                t_result_std[n][k][l] = np.std(result_per_loop)

    result_dict = dict()
    result_dict['methode'] = methode
    result_dict['methode_args'] = methode_args
    result_dict['times_mean'] = t_result_mean
    result_dict['times_std'] = t_result_std
    result_dict['t_end'] = t_end
    result_dict['num_grid_nodes'] = num_nodes
    # result_dict['t_s'] = ts
    # result_dict['repeats'] = repeat
    # result_dict['loops'] = loops
    result_dict['info'] = 'Logs the mean and std of the execution time for all defined methods to simulate the power ' \
                          'grid for different simulation times (t_end) and grid size (num_grid_nodes). num_grid_nodes '\
                          'thereby defines the number of sources and the number of loads ' \
                          '(grid size = 2*num_grid_nodes). '\
                          'Each experiment is executed loops*repeats-times ' \
                          'while the mean and std is calculated based' \
                          'on repeats' \

    return result_dict


if __name__ == '__main__':
    """
    Small execution exxample of the two functions
    """
    parameter_used = dict()
    parameter_used['R_source'] = 0.4
    parameter_used['L_source'] = 2.3e-3
    parameter_used['C_source'] = 10e-6
    parameter_used['L_cabel'] = 2.3e-3
    parameter_used['R_cabel'] = 0.4
    parameter_used['R_load'] = 14
    parameter_used['V_dc'] = 300

    repeat_used = 2
    loops_used = 3

    t_s_used = 1e-4
    t_end_used = [0.001, 0.002, 0.004]  # ,0.005, 0.01, 0.05]

    num_nodes_used = [2, 4, 6, 7, 10]

    time_result = timing_experiment_simulation(repeat=repeat_used, loops=loops_used, num_nodes=num_nodes_used,
                                               t_end=t_end_used, ts=t_s_used,
                                               methode=['control.py', 'control.py', 'scipy_ode', 'scipy_odeint'],
                                               methode_args=['continuous', 'discrete', 'LSODA', 'LSODA'],
                                               parameter=parameter_used)

    plot_result(time_result, num_nodes_used, t_end_used)
