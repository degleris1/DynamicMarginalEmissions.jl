using CarbonNetworks

n_threads = 10

Threads.@threads for i in 1:n_threads
    param = params[i]
    run_expansion_planning(param)
end