include("util.jl")

for (i, cfg) in enumerate(CARBON_TAX_CONFIG)
    run_expansion_planning(cfg, "tax$(i)")
end