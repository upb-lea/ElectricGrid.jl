import json

import pandas as pd

df = pd.read_pickle('sim_time_python_repeats_2_loops_3.pkl.bz2')
data = json.load('sim_time_python_repeats_2_loops_3.json')
asd = 1
#plot_result(time_result, num_nodes_used, t_end_used)