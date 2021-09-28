include("util.jl")

Threads.@threads for i in 1:length(STORAGE_CONFIG)
    println("Running storage model $i")
    run_rodm_dynamic(STORAGE_CONFIG[i], "storage$(i)")
end