function set_config(default, config)
    for (k, v) in default
        if !(k in keys(config))
            config[k] = v
        end
    end
end

DEFAULT_CONFIG = Dict(
    :net_seed => 1,
    :problem_seed => 3,
    :initialization_seed => 5,

    :renewable_penetration => 0.50,
    :demand_growth => 2.5,
    :renewable_archetypes => [0, 1, 3, 5],  # No wind!
    :horizon => 5*365.0,

    :storage_cost_per_mwh => 350.0,
    :charge_efficiency => 0.95,
    :discharge_efficiency => 0.95,
    :charge_rate => 0.25,

    :line_cost_per_mw_mile => 3.0,
    :line_length_min => 40.0,
    :line_length_max => 50.0,

    :emissions_rate => [0.35, 0.45, 1.1, 1.2, 1.3],
    :emissions_weight => 500.0,
)

CHEAP_STORAGE_CONFIG = set_config(DEFAULT_CONFIG, Dict{Symbol, Any}(
    :storage_cost_per_mwh => 100.0,
))