using Dare
using Logging

global_logger(DareLogger(file_name = split(string(@__FILE__), "\\")[end], add_date = true, log_file = true, log_level = Logging.Info))


Timestep = 100 #time step in μs ~ 100μs => 10kHz, 50μs => 20kHz, 20μs => 50kHz
t_final = 0.5 #time in seconds, total simulation run time
ts = Timestep*1e-6
t = 0:ts:t_final # time

fs = 1/ts # Hz, Sampling frequency of controller ~ 15 kHz < fs < 50kHz

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

#setup = create_setup(parameters)

# Define the environment

num_sources = 2
num_loads = 1

env = SimEnv(ts = ts, use_gpu = false, CM = CM, num_sources = num_sources, num_loads = num_loads, 
parameters = parameters, maxsteps = length(t), action_delay = 1)


#example log
@info "env created with" maxsteps = env.maxsteps x0 = env.x0
