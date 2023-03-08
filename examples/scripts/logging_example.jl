using JEG
using Logging

global_logger(JEGLogger(file_name = split(string(@__FILE__), "\\")[end], add_timestamp = false, log_file = true, log_level = Logging.Info, log_level_file = Logging.Debug))


Timestep = 100
t_final = 0.5
ts = Timestep*1e-6
t = 0:ts:t_final

fs = 1/ts

CM = [ 0. 0. 1.
        0. 0. 2
        -1. -2. 0.]

source1 = Dict("pwr" => 200e3,
                "vdc" => 800,
                "fltr" => "LC")

source2 = Dict("pwr" => 150e3,
                "vdc" => 800,
                "fltr" => "LCL")

parameters = Dict()

parameters["source"] = [source1, source2]
parameters["grid"] = Dict("fs" => fs, "phase" => 3, "v_rms" => 230)


num_sources = 2
num_loads = 1

env = ElectricGridEnv(ts = ts, use_gpu = false, CM = CM, num_sources = num_sources, num_loads = num_loads, 
parameters = parameters, maxsteps = length(t), action_delay = 1)


#example log
@info "env created with" maxsteps = env.maxsteps x0 = env.x0 A = env.A
