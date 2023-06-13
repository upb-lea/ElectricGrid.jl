import Base: /, *, +, -, ^

/(a::Number, b::Function) = x -> a / b(x)
/(a::Function, b::Function) = x -> a(x) / b(x)
/(a::Function, b::Number) = x -> a(x) / b
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
    "a" => x->x,
    "L_value_1" => 1,
    "L_value_2" => 3
)

parameters["L"] = funktion(parameters["L_value_1"],parameters["L_value_2"])

f = parameters["L"]

x = [1, 2, 3]

A = [2 1/f f; 
    2 1/f f;
    f f f]

b = Matrix{Function}(undef,size(A))

(rows,columns) = size(A)

# thats a little bit stupid. Because I have to make all numbers to a constant function
for row in 1:rows
    for column in 1:columns
        h = A[row,column]
        if isa(h,Number)
            b[row,column] = x->h
        else
            b[row,column] = h
        end
    end
end

# thats where the magic happens
C(x) = (|>).(x,b)
C(x)