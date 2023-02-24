# Use this file to include utility functions that you want to be available in tests
using Dare 
using JuMP


function populate_params(params, Power, pf, Vrms=230)

    R_load, L_load, _, _ = Parallel_Load_Impedance(Power, pf, Vrms)

    if haskey(params, "load")
        pop!(params, "load")
    elseif haskey(params, "cable")
        pop!(params, "cable")
    end

    params["load"] = Any[
                            Dict{Any, Any}("impedance" => "RL", "R" => R_load, "L" => L_load),
                        ]
    
    return params
end