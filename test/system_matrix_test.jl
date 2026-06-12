using Test
using ElectricGrid
using LinearAlgebra
using Printf

#=
Regression test: the system matrices a user extracts via NodeConstructor/GetSystem,
discretized with the exact-ZOH formula
    A_d = exp(A*ts),   B_d = A^{-1} (A_d - I) B,
must be identical to what an as-constructed ElectricGridEnv stores: env.A/env.B
(continuous) and env.sys.A/env.sys.B (discrete HeteroStateSpace).

This is the "oracle check" pattern from external downstream use (plant extraction
for differentiable rollouts): downstream code depends on the extracted matrices
being exactly the dynamics the env simulates. The env builds Ad/Bd with the same
formula (src/electric_grid_env.jl), so the comparison is expected to be exact;
the thresholds below (1e-9 continuous / 1e-12 discrete) allow only for benign
refactors of the same computation.

Explicit component values (and explicit cable "len") are used so the build
bypasses FilterDesign and the JuMP/Ipopt power-flow solver and stays
deterministic.
=#

@testset "system_matrix_extraction" begin

    CM = [0.0 1.0; -1.0 0.0]
    ts = 1e-4

    make_params(; Rload, vrms, fgrid) = Dict{Any,Any}(
        "source" => Any[Dict{Any,Any}("fltr"=>"L","L1"=>2.3e-3,"R1"=>0.4,"pwr"=>10e3,
                                       "vdc"=>800.0,"control_type"=>"classic","mode"=>"Droop")],
        "load"   => Any[Dict{Any,Any}("impedance"=>"RL","R"=>Rload,"L"=>10e-3)],
        "cable"  => Any[Dict{Any,Any}("R"=>0.1,"L"=>0.1e-3,"C"=>1e-6,"len"=>1.0)],
        "grid"   => Dict{Any,Any}("phase"=>3,"v_rms"=>vrms,"f_grid"=>fgrid))

    configs = (("50Hz/230V R=14", 14.0, 230.0, 50.0),
               ("60Hz/277V R=20", 20.0, 277.0, 60.0),
               ("60Hz/277V R=13", 13.0, 277.0, 60.0))

    for (label, R, vrms, fg) in configs
        # extraction path (fresh params dict each build; NodeConstructor mutates it)
        nc = NodeConstructor(num_sources=1, num_loads=1, CM=CM,
                             parameters=make_params(Rload=R, vrms=vrms, fgrid=fg), ts=ts)
        A, B, C, D = GetSystem(nc)
        A = Float64.(A); B = Float64.(B)
        Ad = exp(A * ts)
        Bd = factorize(A) \ ((Ad - I) * B)

        # as-constructed env path
        env = ElectricGridEnv(num_sources=1, num_loads=1, CM=CM,
                              parameters=make_params(Rload=R, vrms=vrms, fgrid=fg),
                              ts=ts, t_end=0.02, verbosity=0)

        dA  = maximum(abs.(A .- Matrix(env.A)))
        dB  = maximum(abs.(B .- Matrix(env.B)))
        dAd = maximum(abs.(Ad .- Matrix(env.sys.A)))
        dBd = maximum(abs.(Bd .- Matrix(env.sys.B)))
        rho = maximum(abs.(eigvals(Matrix(env.sys.A))))

        @printf("%s   n_states = %d   rho(A_d) = %.6f\n", label, size(A, 1), rho)
        @printf("  max|A  - env.A    | = %.3e   max|B  - env.B    | = %.3e\n", dA, dB)
        @printf("  max|Ad - env.sys.A| = %.3e   max|Bd - env.sys.B| = %.3e\n", dAd, dBd)

        @test dA  < 1e-9
        @test dB  < 1e-9
        @test dAd < 1e-12
        @test dBd < 1e-12
        @test rho < 1.0   # exact-ZOH discretization is stable for these plants
        @test env.sys.Ts == ts
    end
end
