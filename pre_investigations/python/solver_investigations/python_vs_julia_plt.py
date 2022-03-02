import json
import math
import time

import matplotlib.pyplot as plt
import numpy as np


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
        # iterates different methode pairs
        for x in range(len(plt_setting['x_axis'])):
            if plt_setting['x_axis'][x] == 'num_grid_nodes':
                if plt_settings['t_end_show'] is None:
                    for l in range(len(results_python['t_end'])):
                        plt.plot(results_python[plt_setting['x_axis'][x]], time_result_python[:, l], 'o',
                                 label='python')
                        plt.plot(results_julia[plt_setting['x_axis'][x]], time_result_julia[:, l], 'or', label='julia')
                        plt.xlabel(plt_setting['x_axis'][x])
                        plt.ylabel('$mean(execution-time)\,/\,\mathrm{s}$')
                        plt.legend()
                        plt.grid()
                        plt.title(
                            'python:' + plt_setting['methods_python'][m] + ';  julia:' + plt_setting['methods_julia'][
                                m] +
                            ';  t_end:' + str(results_python['t_end'][l]))
                        plt.show()
                        time.sleep(0.5)
                else:
                    for s in range(len(plt_settings['t_end_show'])):
                        l = results_python['t_end'].index(plt_settings['t_end_show'][s])
                        plt.plot(results_python[plt_setting['x_axis'][x]], time_result_python[:, l], 'o',
                                 label='python')
                        plt.plot(results_julia[plt_setting['x_axis'][x]], time_result_julia[:, l], 'or', label='julia')
                        plt.xlabel(plt_setting['x_axis'][x])
                        plt.ylabel('$mean(execution-time)\,/\,\mathrm{s}$')
                        plt.legend()
                        plt.grid()
                        plt.title(
                            'python:' + plt_setting['methods_python'][m] + ';  julia:' + plt_setting['methods_julia'][
                                m] +
                            ';  t_end:' + str(results_python['t_end'][l]))
                        plt.show()
                        time.sleep(0.5)

            if plt_setting['x_axis'][x] == 't_end':
                if plt_settings['num_grid_nodes_show'] is None:
                    for l in range(len(results_python['num_grid_nodes'])):
                        plt.plot(results_python[plt_setting['x_axis'][x]], time_result_python[l, :], 'o',
                                 label='python')
                        plt.plot(results_julia[plt_setting['x_axis'][x]], time_result_julia[l, :], 'or', label='julia')
                        plt.xlabel(plt_setting['x_axis'][x])
                        plt.ylabel('$mean(execution-time)\,/\,\mathrm{s}$')
                        plt.legend()
                        plt.grid()
                        plt.title(
                            'python:' + plt_setting['methods_python'][m] + ';  julia:' + plt_setting['methods_julia'][
                                m] +
                            ';  grid_nodes:' + str(results_python['num_grid_nodes'][l]))
                        plt.show()
                        time.sleep(0.5)
                else:
                    for s in range(len(plt_settings['t_end_show'])):
                        l = results_python['t_end'].index(plt_settings['t_end_show'][s])
                        plt.plot(results_python[plt_setting['x_axis'][x]], time_result_python[:, l], 'o',
                                 label='python')
                        plt.plot(results_julia[plt_setting['x_axis'][x]], time_result_julia[:, l], 'or', label='julia')
                        plt.xlabel(plt_setting['x_axis'][x])
                        plt.ylabel('$mean(execution-time)\,/\,\mathrm{s}$')
                        plt.legend()
                        plt.grid()
                        plt.title(
                            'python:' + plt_setting['methods_python'][m] + ';  julia:' + plt_setting['methods_julia'][
                                m] +
                            ';  grid_nodes:' + str(results_python['num_grid_nodes'][l]))
                        plt.show()
                        time.sleep(0.5)


if __name__ == '__main__':
    f = open('saves/time_result_all_python.json', "r")
    # Reading from file
    time_result_python = json.loads(f.read())

    f = open('saves/time_result_all_julia.json', "r")
    # Reading from file
    time_result_julia = json.loads(f.read())

    plt_setting = dict({'x_axis': ['num_grid_nodes', 't_end'],  # alternative: 't_End'
                        'methods_python': ['control.py', 'env_standalone'],
                        'methods_julia': ['control', 'lsoda'],
                        't_end_show': [0.01],
                        'num_grid_nodes_show': [10]})

    plot_python_vs_julia(time_result_python, time_result_julia, plt_setting)
