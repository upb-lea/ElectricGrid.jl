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
        test.tmp[i] = i
    end
end


run1()
run2()

@timeit timer "Run1" begin
    run1()
end

@timeit timer "Run2" begin
    run2()
end

show(timer)