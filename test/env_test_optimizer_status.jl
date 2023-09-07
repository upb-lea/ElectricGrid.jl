using Test
using ElectricGrid
using JuMP

import ElectricGrid

include("utils.jl")

# Simple one source - one load model that optimizer can solve feasibly
CM = [
        0 1
       -1 0
    ]


@testset "Optimizer status: FEASIBLE (Vrms = 230)" begin

    parameters = PopulateParams(100e3, 0.95)

    env = ElectricGridEnv(ts=0.0001, 
                            CM=CM, 
                            parameters=parameters, 
                            t_end=0.1, 
                            verbosity=2, 
                            action_delay=1)

    @test ElectricGrid.optimizer_status["termination_status"] == LOCALLY_SOLVED  #add other feasible enums from https://jump.dev/JuMP.jl/stable/moi/reference/models/#MathOptInterface.TerminationStatusCode
    @test ElectricGrid.optimizer_status["primal_status"] == FEASIBLE_POINT  # add other feasible status from above

end

@testset "Optimizer status: INFEASIBLE (Vrms = 230)" begin

    parameters = PopulateParams(230e3, 0.95)

    env = ElectricGridEnv(ts=0.0001, CM=CM, parameters=parameters, t_end=0.1, verbosity=2, action_delay=1)

    @test ElectricGrid.optimizer_status["termination_status"] == LOCALLY_INFEASIBLE  #add other feasible enums from https://jump.dev/JuMP.jl/stable/moi/reference/models/#MathOptInterface.TerminationStatusCode
    @test ElectricGrid.optimizer_status["primal_status"] == INFEASIBLE_POINT  # add other feasible status from above

end

# other volatges
@testset "Optimizer status: FEASIBLE (Vrms = 100)" begin

    parameters = PopulateParams(100e2, 0.95, 100)

    env = ElectricGridEnv(ts=0.0001, 
                            CM=CM, 
                            parameters=parameters, 
                            t_end=0.1, 
                            verbosity=2, 
                            action_delay=1)

    # @test 
    ElectricGrid.optimizer_status["termination_status"] == LOCALLY_SOLVED  #add other feasible enums from https://jump.dev/JuMP.jl/stable/moi/reference/models/#MathOptInterface.TerminationStatusCode
    @test ElectricGrid.optimizer_status["primal_status"] == FEASIBLE_POINT  # add other feasible status from above

end


@testset "Optimizer status: INFEASIBLE (Vrms = 100)" begin

    parameters = PopulateParams(100e3, 0.95, 100)

    env = ElectricGridEnv(ts=0.0001, 
                            CM=CM, 
                            parameters=parameters, 
                            t_end=0.1, 
                            verbosity=2, 
                            action_delay=1)

    # @test 
    ElectricGrid.optimizer_status["termination_status"] == LOCALLY_INFEASIBLE  #add other feasible enums from https://jump.dev/JuMP.jl/stable/moi/reference/models/#MathOptInterface.TerminationStatusCode
    @test ElectricGrid.optimizer_status["primal_status"] == INFEASIBLE_POINT  # add other feasible status from above

end

