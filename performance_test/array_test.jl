using TimerOutputs

timer = TimerOutput()

mutable struct Test
    tmp
    list_iterator
end

test = Test(randn(10_000), 1)

function append_tmp(test, value)
    test.tmp[test.list_iterator] = value
    test.list_iterator += 1
end

function run1()
    for i = 1:10_000
        append_tmp(test, i)
    end
    test.list_iterator = 1
end

function run2()
    for i = 1:10_000
        test.tmp[test.list_iterator] = i
        test.list_iterator += 1
    end
    test.list_iterator = 1
end


run1()
run2()

@timeit timer "Run1" begin
    run1()
end

@timeit timer "Run2" begin
    run2()
end


function replace_values_loop!(array::Array, new_values::Array)
    for i in eachindex(array)
        array[i] = new_values[i]
    end
end

tmp2 = randn(10_000)

@timeit timer "Run3" begin
 replace_values_loop!(test.tmp, tmp2)
end


function run4(result)
    for i in 1:10_000
        push!(result, i)
    end
end

result = []

@timeit timer "Run4" begin
    run4(result)
end

show(timer)