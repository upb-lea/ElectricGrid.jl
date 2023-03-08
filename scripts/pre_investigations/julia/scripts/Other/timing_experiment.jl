using DrWatson
@quickactivate "MicroGridSimWithRL"

using JSON
include(srcdir("timing_evaluation.jl"))

#nodes = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ,15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130]
nodes = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ,15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 35]
nodes = [35, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2]
#nodes = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ,15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 35, 40]
nodes = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ,15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 35]

#r = timing_experiment_simulation(5, 5, [2,4,6,8,10], [0.001,0.01,0.03], 3, 1e-4, ["control", "env_without_agent", "env_with_agent", "lsoda"], Dict(), false)
#r = timing_experiment_simulation(3, 3, [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30], [0.03], 6, 1e-4, ["control", "env_without_agent"], Dict(), false)
#r = timing_experiment_simulation(3, 3, [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], [0.003], 1, 1e-4, ["env_without_agent"], Dict(), false)
#r = timing_experiment_simulation(3, 3, nodes, [0.03], 3, 1e-4, ["control", "control1", "control16", "controlCUDA"], Dict(), false)
#r = timing_experiment_simulation(3, 3, nodes, [0.03], 3, 1e-4, ["controlCUDA", "controlCUDA32"], Dict(), false)
#r = timing_experiment_simulation(3, 3, nodes, [0.03], 3, 1e-4, ["control32bit"], Dict(), false)
#r = timing_experiment_simulation(3, 3, nodes, [0.03], 1, 1e-4, ["BS3", "BS3_32", "BS3_CUDA", "BS3_CUDA32"], Dict(), false)
#r = timing_experiment_simulation(3, 3, nodes, [0.003], 1, 1e-5, ["control1", "control2", "control4", "control", "controlCUDA", "control32bit"], Dict(), false)
r = timing_experiment_simulation(1, 1, nodes, [0.003], 1, 1e-5, ["control1", "control"], Dict(), false, false, false)

#open(datadir(savename("julia_", r, "json", allowedtypes = (Vector,), ignores=["info","times_mean","times_std","num_samples"])),"w") do f 
open(datadir(r["name"] * ".json"),"w") do f 
    write(f, JSON.json(r)) 
end