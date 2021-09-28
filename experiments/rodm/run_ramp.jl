include("util.jl")

Threads.@threads for i in 1:length(RAMP_CONFIG)
    println("Running ramping model $i")
    run_rodm_dynamic(RAMP_CONFIG[i], "ramp$(i)")
end