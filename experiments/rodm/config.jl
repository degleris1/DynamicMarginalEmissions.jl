function set_config(default, config)
    for (k, v) in default
        if !(k in keys(config))
            config[k] = v
        end
    end

    return config
end

DEFAULT_CONFIG = Dict(
    :datadir => "~/Data/carbon_networks/",
    :multithread => false,
    :storage_percentage => 0.001,
    :charge_rate => 0.25,
    :charge_efficiency => 0.95,
    :discharge_efficiency => 0.95,
    :carbon_tax => 0.00,
    :coal_ramping_rate => 1.5,
    :unit_commitment => false,
    :coal_min => 0.40,
)

RAMP_CONFIG = [
    set_config(DEFAULT_CONFIG, Dict{Symbol, Any}(
        :coal_ramping_rate => ρ,
    ))
    for ρ in [0.05, 0.10, 0.20, 0.40]
]

STORAGE_CONFIG = [
    set_config(DEFAULT_CONFIG, Dict{Symbol, Any}(
        :storage_percentage => s_rel,
        :carbon_tax => 0.04,
    ))
    for s_rel in [0.0, 0.02, 0.05, 0.10]
]

UNIT_CONFIG = set_config(DEFAULT_CONFIG, Dict{Symbol, Any}(
    :unit_commitment => true,
))