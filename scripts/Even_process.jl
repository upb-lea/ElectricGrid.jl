using DrWatson
@quickactivate "dare"

npast = 4 #Past Series size
nfuture = 5 #Future series size

N = 5 #Number of training samples

scale = 1 #bandwidth

#-------------------------------------------------------------

window_size = npast + nfuture

series_length = N + window_size - 1

series = zeros(series_length)

state = rand(0:1)

for i in 1:series_length

    global state

    if state == 1

        series[i] = 1
        state = 0

    else

        series[i] = rand(0:1)
        state = series[i]
    end
end

Gx, Gy, index_map = series_Gxy(series, scale, npast, nfuture)