import json
import math
import time

import matplotlib.pyplot as plt
import numpy as np


def run_plt(x_py, mean_py, std_py, x_jl, mean_jl, std_jl, x_label, title):

    #plt.figure(dpi=200, figsize=(15,9))
    plt.plot(x_py, mean_py, 'xb', label='python_mean')
    plt.plot(x_py, mean_py, 'b')

    plt.plot(x_py, (mean_py + std_py), '--b', linewidth=0.5, label='python_std')
    plt.plot(x_py, (mean_py - std_py), '--b', linewidth=0.5)
    #plt.xticks(x_py)
    plt.fill_between(x_py, (mean_py + std_py), (mean_py - std_py), facecolor='b', alpha=0.25)

    plt.plot(x_jl, mean_jl, 'r')
    plt.plot(x_jl, mean_jl, 'xr', label='julia_mean')
    plt.plot(x_jl, (mean_jl + std_jl), '--r', linewidth=0.5, label='julia_std')
    plt.plot(x_jl, (mean_jl - std_jl), '--r', linewidth=0.5)

    plt.fill_between(x_jl, (mean_jl + std_jl), (mean_jl - std_jl), facecolor='r', alpha=0.25)
    plt.xlabel(x_label)
    plt.ylabel('$execution time\,/\,\mathrm{s}$')
    plt.legend()
    plt.grid()
    plt.title(title)
    plt.show()
    time.sleep(1)


def plot_python_vs_julia(results_python: dict = None, results_julia: dict = None, plt_settings: dict = None):
    """
    Function to plot the results of time comparison for scaling issue of the power grid for different number of nodes
    and runtimes.
    Will loop over the in plt_Setting defined x-axis and mehtode-pairs
    :param results_python: Includes python measured data:
            results['methode'] -> Used methode (e.g. control.py; scipy;... - more to be defined)
            results['methode_args'] -> Additional arguments of the methode: e.g. control.py -> continuous/discrete
                                                                                     scipy -> solver like LSODA,...
            results['times'] -> list of time results of the experiments in
                                list_of_methods[numpy.array(len(nodes),len(t_end_vec))]
    :param results_julia: Includes python measured data: (similar to python)
    :param plt_settings: plt_settings[x_axis] -> num_grid_nodes and/or t_end
                         plt_settings[methods_python] -> list of methods to plot from python
                         plt_settings[methods_python] -> list of methods to plot from julia
                         plt_settings[t_end_show] -> list of t_end times to be plotted to x_label num_grid_nodes (if
                                                     None, all t_End are shown in different plots)
                         plt_settings[num_grid_nodes_show] -> list of num_grid_nodes times to be plotted to
                                                              x_label t_end (if None, all num_nodes are shown in
                                                              different plots)

                         Will compare depending on the order of defined mehtods and throw error if not the same length

                         Example:
                         plt_setting = dict({'x_axis': ['num_grid_nodes', 't_end'],
                        'methods_python': ['control.py', 'env_standalone'],
                        'methods_julia': ['control', 'lsoda'],
                        't_end_show': [0.01],
                        'num_grid_nodes_show': [10]})

                        -> Will compare control python and julia AND env_standalone from python to lsoda from julia
                        for num_grid_nodes and t_end.
                        For num_Grid_nodes as x-axis only t_end = 0.01 is shown and for t_end as x-axis only
                        num_Grid_nodes = 10 is shown
                        --> 4 plots
    """

    if len(plt_setting['methods_julia']) != len(plt_setting['methods_julia']):
        raise ValueError('Please define methode list of same length for python and julia')

    for m in range(len(plt_setting['methods_julia'])):

        time_result_python = np.array(results_python['times_mean'])[
            results_python['methode'].index(plt_setting['methods_python'][m])]
        time_result_julia = np.array(results_julia['times_mean'])[
            results_julia['methods'].index(plt_setting['methods_julia'][m])]
        sdt_result_python = np.array(results_python['times_std'])[
            results_python['methode'].index(plt_setting['methods_python'][m])]
        std_result_julia = np.array(results_julia['times_std'])[
            results_julia['methods'].index(plt_setting['methods_julia'][m])]
        # iterates different methode pairs
        for x in range(len(plt_setting['x_axis'])):
            if plt_setting['x_axis'][x] == 'num_grid_nodes':
                if plt_settings['t_end_show'] is None:
                    for l in range(len(results_python['t_end'])):
                        title = 'python:' + plt_setting['methods_python'][m] + ';  julia:' + \
                                plt_setting['methods_julia'][
                                    m] + ';  t_end:' + str(results_python['t_end'][l])

                        run_plt(results_python[plt_setting['x_axis'][x]], time_result_python[:, l],
                                sdt_result_python[:, l], results_julia[plt_setting['x_axis'][x]],
                                time_result_julia[:, l], std_result_julia[:, l],
                                plt_setting['x_axis'][x], title)

                else:
                    for s in range(len(plt_settings['t_end_show'])):
                        l = results_python['t_end'].index(plt_settings['t_end_show'][s])

                        title = 'python:' + plt_setting['methods_python'][m] + ';  julia:' + \
                                plt_setting['methods_julia'][
                                    m] + ';  t_end:' + str(results_python['t_end'][l])

                        run_plt(results_python[plt_setting['x_axis'][x]], time_result_python[:, l],
                                sdt_result_python[:, l], results_julia[plt_setting['x_axis'][x]],
                                time_result_julia[:, l], std_result_julia[:, l],
                                plt_setting['x_axis'][x], title)

            if plt_setting['x_axis'][x] == 't_end':
                if plt_settings['num_grid_nodes_show'] is None:
                    for l in range(len(results_python['num_grid_nodes'])):
                        title = 'python:' + plt_setting['methods_python'][m] + ';  julia:' + \
                                plt_setting['methods_julia'][m] + ';  grid_nodes:' + \
                                str(results_python['num_grid_nodes'][l])

                        run_plt(results_python[plt_setting['x_axis'][x]], time_result_python[l, :],
                                sdt_result_python[l, :], results_julia[plt_setting['x_axis'][x]],
                                time_result_julia[l, :], std_result_julia[l, :],
                                plt_setting['x_axis'][x], title)


                else:
                    for s in range(len(plt_settings['num_grid_nodes_show'])):
                        l = results_python['num_grid_nodes'].index(plt_settings['num_grid_nodes_show'][s])

                        title = 'python:' + plt_setting['methods_python'][m] + ';  julia:' + \
                                plt_setting['methods_julia'][m] + ';  grid_nodes:' + \
                                str(results_python['num_grid_nodes'][l])

                        run_plt(results_python[plt_setting['x_axis'][x]], time_result_python[l, :],
                                sdt_result_python[l, :], results_julia[plt_setting['x_axis'][x]],
                                time_result_julia[l, :], std_result_julia[l, :],
                                plt_setting['x_axis'][x], title)


if __name__ == '__main__':
    f = open('saves/time_result_all_python_1cm.json', "r")
    # Reading from file
    time_result_python = json.loads(f.read())

    f = open('saves/time_results_all_julia_1cm.json', "r")
    # Reading from file
    time_result_julia = json.loads(f.read())

    plt_setting = dict({'x_axis': ['num_grid_nodes'],  # alternative: 't_end'
                        'methods_python': ['scipy_odeint'],#, 'env_standalone'],#, 'scipy_odeint'],
                        #                   'scipy_ode', 'scipy_solve_ivp'],
                        'methods_julia': ['lsoda'],#, 'env_without_agent'],#, 'lsoda'],
                        #                  'lsoda', 'lsoda', 'lsoda'],
                        #'methods_python': ['scipy_solve_ivp'],
                        #'methods_julia': ['lsoda'],
                        't_end_show': [0.03],  # ,
                        'num_grid_nodes_show': None#[18, 20]
                         })

    plot_python_vs_julia(time_result_python, time_result_julia, plt_setting)
