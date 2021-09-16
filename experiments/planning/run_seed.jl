include("util.jl")

for (i, cfg) in enumerate(SEED_CONFIG)
    run_expansion_planning(cfg, "seed$(i)")
end