import json
import math
import time

import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':

    filenames = [  # 'time_result_python_nodes[2-70]_ts1e-5_control_stepwise.json',
        # 'time_result_python_nodes[2-24]_ts1e-5_RK45,RK23.json',
        # 'time_result_python_nodes[2-24]_ts1e-5_RK45,RK23_fixedstepsize.json',
        # 'time_result_python_nodes[2-40]_expm.json',
        # 'time_result_python_nodes[2-70]_ts1e-5_control.json'
        # 'time_result_python_nodes[2-70]_ts1e-5_controlCUDA_stepwise.json',
        #'time_result_python_nodes[2-70]_ts1e-5_env_agent.json'
    ]

    for filename in filenames:
        f = open(
            'saves/Jan_officePC/' + filename, "r")

        # Reading from file
        time_result_python = json.loads(f.read())
        # time_result_julia = json.loads(f.read())

        for cc in range(len(time_result_python['methode_args'])):
            # for cc in range(len(time_result_julia['methods'])):
            y_debug_python = np.array(time_result_python['times_mean'][cc])
            # y_debug_julia = np.array(time_result_julia['times_mean'][cc])

            # plt.plot(time_result_python['num_grid_nodes'], y_debug_python, label=str(num[cc]) + 'threads')
            # if cc == 200:
            # if time_result_python['methode_args'][cc] == 'cuda' and filename == 'time_result_python_nodes[2-70]_ts1e-5_control_stepwise.json':
            #    #plt.plot(time_result_python['num_grid_nodes'], y_debug_python, label='control')
            #    asd = 1
            # else:
            plt.plot(time_result_python['num_grid_nodes'], y_debug_python,
                     label=str(time_result_python['methode_args'][cc]) + '')
        # plt.plot(time_result_julia['num_grid_nodes'], y_debug_julia, label=time_result_julia['methods'][cc]+'.jl')

    plt.legend()
    plt.ylim([0, 0.7])
    # plt.xlim([2, 35])
    plt.grid()
    plt.ylabel('time mean')
    plt.xlabel('nodes')
    plt.title(filename)
    plt.show()
    asd = 1

    """
    filenames = [#'julia_control_2-70_stepwise_neu.json',
       #'julia_control_2-70_stepwise.json',
                    #'julia_BS3,DP5_2-24.json',
                    #'julia_BS3,DP5_2-24_fixedstepsize.json',
                    #'julia_expm_2-40.json',
                    'julia_control_2-70.json',
       'julia_controlCUDA_2-70.json',
       #'julia_env_agent_2-70.json',
               ]
    
    for filename in filenames:
       f = open(
           'saves/Jan_officePC/' + filename, "r")
    
       # Reading from file
       # time_result_python = json.loads(f.read())
       time_result_julia = json.loads(f.read())
    
       # for cc in range(len(time_result_python['methode_args'])):
       for cc in range(len(time_result_julia['methods'])):
           # y_debug_python = np.array(time_result_python['times_mean'][cc])
    
           if time_result_julia['methods'][cc] == 'controlCUDA' and filename == 'julia_control_2-70.json':
               asd = 1
           else:
               y_debug_julia = np.array(time_result_julia['times_mean'][cc])
    
           # plt.plot(time_result_python['num_grid_nodes'], y_debug_python, label=str(num[cc]) + 'threads')
    
           plt.plot(time_result_julia['num_grid_nodes'], y_debug_julia, label=time_result_julia['methods'][cc] + '.jl')
    
    plt.legend()
    plt.ylim([0, 0.3])
    plt.xlim([2, 35])
    plt.grid()
    plt.ylabel('time mean')
    plt.xlabel('nodes')
    plt.title(filename)
    plt.show()
    asd = 1
"""
