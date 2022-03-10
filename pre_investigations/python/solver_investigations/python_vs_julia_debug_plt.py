import json
import math
import time

import matplotlib.pyplot as plt
import numpy as np



if __name__ == '__main__':
    #f = open('saves/AAAAtime_result_all_python_debug_nodes8t_end0.03.json', "r")
    f = open('saves/time_result_all_python_debug_nodes8t_end0.03.json', "r")
    # Reading from file
    time_result_python = json.loads(f.read())

    #f = open('saves/julia_debug_lsoda,DP5,BS3_8_0.03_cm=1_loops=3_repeat=3.json', "r")
    f = open('saves/julia_debug_control_8_0.03_cm=3_loops=1_repeat=1.json', "r")
    # Reading from file
    time_result_julia = json.loads(f.read())

    plt_setting = dict({'c': [0,1,2]})


    for cc in range(len(plt_setting['c'])):

        y_debug_julia = np.array(time_result_julia['y_debug']['control_c'+str(cc+1)+'_nodes8_tend0.03'])
        y_debug_python = np.array(time_result_python['y_debug']['control.py_c'+str(cc)+'_nodes8_tend0.03'])

        diff = y_debug_julia - y_debug_python

        i, j = np.unravel_index(diff.argmax(), diff.shape)
        print(diff[i][j])

        #plt.figure(dpi=300, figsize=(15,9))
        plt.plot(y_debug_julia[i], 'r', label='julia')
        #plt.plot(y_debug_julia[i], 'xr', label='julia')

        plt.plot(y_debug_python[i], '--b', label='python')
        #plt.plot(y_debug_python[i], 'xb', label='python')

        plt.title('diff: '+str(diff[i][j])+ 'c: ' + str(cc))
        plt.legend()
        #plt.ylim([87,88])
        plt.grid()
        plt.ylabel('state (max_diff)')
        plt.show()
        asd = 1



