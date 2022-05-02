using Dates
using BSON
using StatsBase: mean
using LinearAlgebra

"""
"""
function load_results(p)
    return Dict(
	    DateTime(f[1:end-5], "yy-mm-dd-HH") => BSON.load(joinpath(p, f), @__MODULE__) 
	    for f in filter(f -> f[1:2] in ["04", "18"], readdir(p))
    )
    end


"""
    get_nodal_mefs(r, whichdates=d -> true; hybrid_mode=true, observed=false)

Get the nodal MEFs, filtering so that only MEFs from `whichdates` are returned.
"""
function get_nodal_mefs(r, whichdates=d -> true; hybrid_mode=true, observed=false)
    is_dynamic = (ndims(first(r)[2][:λ]) == 3)

    dates = sort(collect(keys(r)))
    is_valid = [r[d][:status] for d in dates] .== "OPTIMAL"

    # Get total mefs
    function get_total_mef(rd)
        if is_dynamic && hybrid_mode
            return reduce(hcat, rd[:λ_static])
        elseif is_dynamic
            return dropdims(sum(rd[:λ], dims=2), dims=2)
        else
            return rd[:λ]
        end
    end

    function get_observed_mef(rd)
        if is_dynamic && hybrid_mode
            return diag(rd[:λ_static])
        elseif is_dynamic
            return vcat([transpose(diag(rd[:λ][k, :, :])) for k in 1:size(rd[:λ])[1]]...)
        else
            return diag(rd[:λ])
        end
    end

    if observed
        mefs = [v ? get_observed_mef(r[d]) : missing for (d, v) in zip(dates, is_valid)]
    else
        mefs = [v ? get_total_mef(r[d]) : missing for (d, v) in zip(dates, is_valid)]
    end
    # Expand dates
    all_dates = [d .+ Hour.(0 : size(m, 2) - 1) for (d, m) in zip(dates, mefs) if !ismissing(m)]
    mefs = [m for m in mefs if !ismissing(m)]

    # Join lists of lists
    mefs = reduce(hcat, mefs)
    all_dates = reduce(vcat, all_dates)

    # Filter by hour
    mefs_hr = mefs[:, map(whichdates, all_dates)]

    return mefs_hr
end

"""
    get_average_nodal_mefs(r, whichdates=d -> true; hybrid_mode=true)

Get the nodal MEFs averaged across `whichdates`.
"""
function get_average_nodal_mefs(r, whichdates=d -> true; hybrid_mode=true)
	mefs = get_nodal_mefs(r, whichdates; hybrid_mode=hybrid_mode)
	return [mean(skipmissing(mefs[i, :])) for i in 1:size(mefs, 1)]
end


"""
    get_demand(r)

Returns a matrix of demand of size n_nodes x n_timesteps
"""
function get_demand(r)

	dates = sort(collect(keys(r)))
	is_valid = [r[d][:status] for d in dates] .== "OPTIMAL"

	demand = hcat([v ? hcat(r[d][:d]...) : missing for (d,v) in zip(dates, is_valid)]...)

	
	return demand
end

"""
"""
function get_total_emissions(r, co2_rates)
    dates = sort(collect(keys(r)))
	is_valid = [r[d][:status] for d in dates] .== "OPTIMAL"

	g_opt = hcat([ v ? hcat(r[d][:g]...) : missing for (d, v) in zip(dates, is_valid)]...)
    E = co2_rates * g_opt

	return E
end