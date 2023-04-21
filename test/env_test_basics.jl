using Test
using ElectricGrid

CM = [0. 1.
    -1. 0.]

parameters = Dict{Any, Any}(
        "source" => Any[
                        Dict{Any, Any}("pwr" => 200e3, "control_type" => "RL", "mode" => "user_def", "fltr" => "L", "L1" => 1e-4, "R1" => 1.1e-3, "i_limit"=>1),
                        ],
            "load"   => Any[
                        Dict{Any, Any}("impedance" => "RL", "R" => 2.64, "L" => 0.006, "v_limit"=>10000, "i_limit"=>10e8),
                        ],
        "cable"   => Any[
                        Dict{Any, Any}("R" => 1e-3, "L" => 1e-4, "C" => 1e-4, "i_limit" => 10e8,),
                        ],
        "grid" => Dict{Any, Any}("ramp_end" => 0.0)
        )

@testset "env_abort" begin
    """
    Checks if the env puts done = true if limit is exceeded
    """
    
    env = ElectricGridEnv(CM = CM, parameters = parameters, verbosity = 0, action_delay = 0)
    reset!(env)
    env([1,1,1])
    env([1,1,1])
    @test env.done == true
    
end

@testset "env_state_space" begin

    env = ElectricGridEnv(CM = CM, parameters = parameters, verbosity = 0, action_delay = 0)
    @test length(env.state_space) == 15
end

@testset "env_minimal_action_space" begin
    """
    Checks if it is possible to define an env least amount of arguments and checks the length of action space based on default values on case of 1 source
    """
    env = ElectricGridEnv(num_sources = 1, num_loads = 1)
    @test length(env.action_space) == 3
end