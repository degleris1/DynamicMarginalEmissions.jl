
function get_total_emissions(p, generator_devices, emissions_rates)
    total_emissions = 0.0

    for (g, e) in zip(generator_devices, emissions_rates)
        total_emissions += p[g, :]' * e
    end
    
    return total_emissions
end

function get_lmes(pmp, pmp_result, demand_devices, gen_devices, emissions_rates)
    T = get_time_horizon(pmp)
    (; p) = pmp_result

    # Compute the gradient `dE/dp` of total emissions with respect to `p`
    dz = zeros(get_dims(pmp))
    for (g, e) in zip(gen_devices, emissions_rates)
        dz[get_device_primary_var_inds(pmp, g)] = e
    end

    # Multiply by `dz/dD`, the Jacobian of optimal variables as a function of demand
    dD = pullback_pmin(pmp, pmp_result, demand_devices, dz)
    return reshape(dD, length(demand_devices), T)
end