import matplotlib.pyplot as plt
import numpy as np
import math

def plot_result(results: dict = None, nodes: list = None, t_end_vec: list = None):
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
    if nodes is None:
        nodes = results['num_grid_nodes']

    if t_end_vec is None:
        t_end_vec = results['t_end']
    for m in range(len(results['methode'])):
        time_result = np.array(results['times_mean'])[m]

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
        if isinstance(results['methode_args'][m], dict):
            fig.suptitle(results['methode'][m] + " using " + results['methode_args'][m]['integrator'] +
                         " with solver  " + results['methode_args'][m]['methode'])
        else:
            fig.suptitle(results['methode'][m] + " using " + results['methode_args'][m])
        plt.show()