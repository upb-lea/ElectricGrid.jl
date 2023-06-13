import Base: /, *, +, -, ^

/(a::Number, b::Function) = x -> a / b(x)
/(a::Function, b::Function) = x -> a(x) / b(x)
*(a::Function, b::Number) = x -> a(x) * b
*(b::Number, a::Function) = x -> a(x) * b
*(a::Function, b::Function) = x -> a(x) * b(x)
+(a::Function, b::Number) = x -> a(x) + b
+(b::Number, a::Function) = x -> a(x) + b
+(a::Function, b::Function) = x -> a(x) + b(x)
-(a::Function, b::Number) = x -> a(x) - b
-(b::Number, a::Function) = x -> b - a(x)
-(a::Function, b::Function) = x -> a(x) - b(x)
-(a::Function) = x -> -a(x)
^(a::Function, b::Number) = x -> a(x)^b

function funktion(v₁,v₂)
    return x -> v₁ * x^v₂
end

parameters = Dict(
    "L_value_1" => 1,
    "L_value_2" => 1,
    # "L" => funktion(parameters["L_value_1"],parameters["L_value_2"]),
    "L" => x->parameters["L_value_1"]*x^parameters["L_value_2"]
)


f = parameters["L"]


A = [1/f 4 -f
    1/f 5 -f]

function B(x)
    rows, columns = size(A)
    H = zeros(rows, columns)
    for row = 1:rows
        for column = 1:columns
            if isa(A[row, column], Function)
                H[row, column] = A[row, column](x[row])
            else
                H[row, column] = A[row, column]
            end
        end
    end
    return H
end
# statevariable:
x = [1,2]
@show B(x)