using TimerOutputs

timer = TimerOutput()

mutable struct Abc
    matrix
end

abc = Abc(rand(2,3,4))

function roll1(abc)
    for i = 1:10_000
        abc.matrix[2, :, 1:end-1] = abc.matrix[2, :, 2:end]
    end
end

function roll2(matrix)
    for i = 1:10_000
        @views matrix[2, :, 1:end-1] = matrix[2, :, 2:end]
    end
end

roll1(abc)
roll2(abc.matrix)

abc = Abc(rand(2,3,4))
@timeit timer "Run1" begin
    roll1(abc)
end

abc = Abc(rand(2,3,4))
@timeit timer "Run2" begin
    roll2(abc.matrix)
end

show(timer)


