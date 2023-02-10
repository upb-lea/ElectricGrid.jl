using Test
using Dare

@testset "env_abort" begin
    """
    Checks if the env puts done = true if limit is exceeded
    """
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
end