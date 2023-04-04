using Revise
using CarbonNetworks
using Plots

using CarbonNetworks: sensitivity_demand_check, sensitivity_var_check


"""
    plot_sensitivity_check(dnet, d_dyn, node, t; varName="obj", unit=1, npoints=10, rel_inc=1e-1)   

Plot the result of the sensitivity check, i.e. actual optimal values (of objective or primal variables), along
with those estimated from implicit diff. 
"""
function plot_sensitivity_check(dnet, d_dyn, node, t; varName="obj", unit=1, npoints=10, rel_inc=1e-1)

    if varName=="obj"
        opt_vals, estimated_vals, rel_value = sensitivity_demand_check(dnet, d_dyn, node, t; npoints=npoints, rel_inc=rel_inc);
    else
        opt_vals, estimated_vals, rel_value = sensitivity_var_check(
            dnet, d_dyn, node, varName, unit, t; npoints=npoints, rel_inc=rel_inc
        )
    end
    _plot_check(opt_vals, estimated_vals, rel_value)
end


"""
    _plot_check(values, values_est, rel_value)

Helper function for plotting sensitivity check
"""
function _plot_check(values, values_est, rel_value)

    plt = plot(
        rel_value, values, marker=4, ylabel="variable value", xlabel="relative parameter value", label="Real",
        linewidth = 2, ls=:solid
    )
    plot!(rel_value, values_est, label="Estimates from diff", lw=2, ls=:dash, color="orange")
    plot!(legend=:best)

    return plt
end


