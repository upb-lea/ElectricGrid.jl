using Test
using Dare
using JuMP

include("utils.jl")

# Simple one source - one load model that optimiser can solve feasibly
CM = [  
    0. 1.
    -1. 0.
]


parameters = Dict{Any, Any}(
    "source" => Any[
                    Dict{Any, Any}("pwr" => 200e3, "mode" => 1),
                    # Dict{Any, Any}("pwr" => 100e3, "mode" => 4),
                    ],
    "grid" => Dict{Any, Any}("ramp_end" => 0.04)
)


populate_params(parameters, 100e3, 0.95, 230)


@testset "Optimizer status: FEASIBLE" begin
    env = SimEnv(ts = 0.0001, CM = CM, parameters = parameters, t_end = 0.1, verbosity = 2, action_delay = 1)
    # Dare.optimizer_status

    @test Dare.optimizer_status["termination_status"] == LOCALLY_SOLVED  #add other feasible enums from https://jump.dev/JuMP.jl/stable/moi/reference/models/#MathOptInterface.TerminationStatusCode
    @test Dare.optimizer_status["primal_status"] == FEASIBLE_POINT  # add other feasible status from above

end

@testset "Optimizer status: FEASIBLE" begin
    populate_params(parameters, 230e3, 0.95, 230)
    env = SimEnv(ts = 0.0001, CM = CM, parameters = parameters, t_end = 0.1, verbosity = 2, action_delay = 1)
    # Dare.optimizer_status

    @test Dare.optimizer_status["termination_status"] == LOCALLY_SOLVED  #add other feasible enums from https://jump.dev/JuMP.jl/stable/moi/reference/models/#MathOptInterface.TerminationStatusCode
    @test Dare.optimizer_status["primal_status"] == FEASIBLE_POINT  # add other feasible status from above

end
   



    
    
