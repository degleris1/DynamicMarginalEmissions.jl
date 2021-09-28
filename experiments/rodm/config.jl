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
    :storage_percentage => 0.001,
    :charge_rate => 0.25,
    :charge_efficiency => 0.95,
    :discharge_efficiency => 0.95,
    :carbon_tax => 0.00,
    :coal_ramping_rate => 0.10,
)

STORAGE_CONFIG = set_config(DEFAULT_CONFIG, Dict{Symbol, Any}(
    :storage_percentage => 0.05,
    :carbon_tax => 0.04,
))