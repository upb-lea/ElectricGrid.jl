


mutable struct NodeConstructor
    bar
    baz::Int
    qux::Float64
end

function NodeConstructor(;bar="hallo",baz=2,qux=3.0)
    NodeConstructor(bar,baz,qux)
end