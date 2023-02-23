using Test
using Dare
using JuMP

@testset "Optimizer status " begin
# Simple one source - one load model that optimiser can solve feasibly
    CM = [  
            0. 1.
            -1. 0.
        ]

    R_load, L_load, _, _ = Parallel_Load_Impedance(100e3, 0.95, 230)

    parameters = Dict{Any, Any}(
            "source" => Any[
                            Dict{Any, Any}("pwr" => 200e3, "mode" => 1),
                            # Dict{Any, Any}("pwr" => 100e3, "mode" => 4),
                            ],
            "load"   => Any[
                            Dict{Any, Any}("impedance" => "RL", "R" => R_load, "L" => L_load),
                            ],
            "grid" => Dict{Any, Any}("ramp_end" => 0.04)
        )

    env = SimEnv(ts = 0.0001, CM = CM, parameters = parameters, t_end = 0.1, verbosity = 0, action_delay = 1)
    Dare.optimizer_status

    @test Dare.optimizer_status["termination_status"] == LOCALLY_SOLVED
    @test Dare.optimizer_status["primal_status"] == FEASIBLE_POINT

end
   



    
    
