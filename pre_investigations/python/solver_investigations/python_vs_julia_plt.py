import json
import math

import matplotlib.pyplot as plt
import numpy as np


def plot_python_vs_julia(results_python: dict = None, results_julia: dict = None):
    """
    Function to plot the results of time comparison for scaling issue of the power grid for different number of nodes
    and runtimes.
    :param results_python: Includes measured data:
            results['methode'] -> Used methode (e.g. control.py; scipy;... - more to be defined)
            results['methode_args'] -> Additional arguments of the methode: e.g. control.py -> continuous/discrete
                                                                                     scipy -> solver like LSODA,...
            results['times'] -> list of time results of the experiments in numpy.array(len(nodes),len(t_end_vec))
    """

    nodes = results_python['num_grid_nodes']

    t_end_vec = results_python['t_end']
    for m in range(len(results_python['methode'])):
        time_result_python = np.array(results_python['times_mean'])[m]
        if m > 2:
            time_result_julia = np.array(results_julia['times_mean'])[2]
        else:
            time_result_julia = np.array(results_julia['times_mean'])[m]

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
                ax[col].step(nodes, time_result_python[:, l], label='python')
                ax[col].step(nodes, time_result_julia[:, l],'r', label='julia')
                # lt.plot(t,result[:steps,0], label = 'i1')
                if col == round(len(t_end_vec) / 2) - 1:
                    ax[col].set_xlabel('nodes')
                if row == 0:
                    ax[col].set_ylabel('$mean(execution-time)\,/\,\mathrm{s}$')
                    ax[col].legend()
                ax[col].title.set_text(f't_end =  {t_end_vec[l]} s')
                ax[col].grid()
            else:
                ax[col, row].step(nodes, time_result_python[:, l], label='P')
                ax[col].step(nodes, time_result_julia[:, l], 'r', label='P')
                # lt.plot(t,result[:steps,0], label = 'i1')
                if col == round(len(t_end_vec) / 2) - 1:
                    ax[col, row].set_xlabel('nodes')
                if row == 0:
                    ax[col, row].set_ylabel('$mean(execution-time)\,/\,\mathrm{s}$')
                ax[col, row].title.set_text(f't_end =  {t_end_vec[l]} s')
                ax[col, row].grid()

            col += 1
        if isinstance(results_python['methode_args'][m], dict):
            fig.suptitle(results_python['methode'][m] + " using " + results_python['methode_args'][m]['integrator'] +
                         " with solver  " + results_python['methode_args'][m]['methode'])
        else:
            fig.suptitle(results_python['methode'][m] + " using " + results_python['methode_args'][m])
        plt.show()

if __name__ == '__main__':
    f = open('saves/time_result_all_python.json', "r")
    # Reading from file
    time_result_python = json.loads(f.read())

    f = open('saves/time_result_all_julia.json', "r")
    # Reading from file
    time_result_julia = json.loads(f.read())

    plot_python_vs_julia(time_result_python, time_result_julia)
