include("util.jl")

for (i, cfg) in enumerate(PARETO_CONFIG)
    run_expansion_planning(cfg, "pareto$(i)")
end