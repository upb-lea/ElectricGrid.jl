# Use this file to include utility functions that you want to be available in tests
using Dare 
using JuMP

const DEFAULT_PARAMS = Dict{Any, Any}(
                    "source"    => Any[
                                    Dict{Any, Any}("pwr" => 200e3, "mode" => 1),
                                    ],
                    "grid"      => Dict{Any, Any}("ramp_end" => 0.04)
                )

function PopulateParams(power, pf, Vrms=230)

    R_load, L_load, _, _ = Parallel_Load_Impedance(power, pf, Vrms)

    # if haskey(params, "load")
    #     pop!(params, "load")
    # elseif haskey(params, "cable")
    #     pop!(params, "cable")
    # end
    params = copy(DEFAULT_PARAMS)

    params["load"] = Any[
                            Dict{Any, Any}("impedance" => "RL", "R" => R_load, "L" => L_load),
                        ]
    
    return params
end