import json
import math
import time

import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':

    num = [1, 8, 16]
    f = open(
        'saves/time_result_all_python_nodes[2-60]t_end0.03control.py[1,8,16].json', #Jan
        #'saves/time_result_all_python_nodes[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 35, 40]t_end0.03control.py[1, 4, 8, 16].json', #DW
        #'saves/julia_control,control1,control16,controlCUDA_2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,35,40,45,50,55,60_0.03_cm=3_loops=3_repeat=3.json',
        #'saves/julia_control,control1,control16_80,85,90,95,100,105,110_0.03_cm=3_loops=3_repeat=3.json',
        "r")

    # Reading from file
    time_result_python = json.loads(f.read())
    #time_result_julia = json.loads(f.read())

    for cc in range(len(time_result_python['methode_args'])):
    #for cc in range(len(time_result_julia['methods'])):
        y_debug_python = np.array(time_result_python['times_mean'][cc])
        #y_debug_julia = np.array(time_result_julia['times_mean'][cc])

        #plt.plot(time_result_python['num_grid_nodes'], y_debug_python, label=str(num[cc]) + 'threads')
        if cc == 2:
            plt.plot(time_result_python['num_grid_nodes'], y_debug_python,'--', label=str(num[cc]) + 'threads')
        else:
            plt.plot(time_result_python['num_grid_nodes'], y_debug_python, label=str(num[cc]) + 'threads')
        #plt.plot(time_result_julia['num_grid_nodes'], y_debug_julia, label=time_result_julia['methods'][cc]+'.jl')

    """
    f = open('saves/julia_controlCUDA,controlCUDA32_2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,35,40,45,50,55,60_0.03_cm=3_loops=3_repeat=3.json',
        "r")
    time_result_julia = json.loads(f.read())

    y_debug_julia = np.array(time_result_julia['times_mean'][0])
    plt.plot(time_result_julia['num_grid_nodes'], y_debug_julia, label='CUDA64.jl')
    y_debug_julia = np.array(time_result_julia['times_mean'][1])
    plt.plot(time_result_julia['num_grid_nodes'], y_debug_julia, label='CUDA32.jl')
    """
    plt.legend()
    plt.ylim([0, 2])
    plt.xlim([2, 35])
    plt.grid()
    plt.ylabel('time mean')
    plt.xlabel('nodes')
    plt.show()
    asd = 1
