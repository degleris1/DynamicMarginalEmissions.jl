using Dates
using BSON
using StatsBase: mean
using LinearAlgebra
using LightGraphs
using Statistics


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
function get_demand(r; n_nodes=243, n_timesteps=24)

	dates = sort(collect(keys(r)))
	is_valid = [r[d][:status] for d in dates] .== "OPTIMAL"

	demand = hcat([v ? hcat(r[d][:d]...) : Array{Union{Missing, String}}(missing, n_nodes, n_timesteps) for (d,v) in zip(dates, is_valid)]...)

	
	return demand
end

"""
"""
function get_total_emissions(r, co2_rates; n_gens=309, n_timesteps=24)
    dates = sort(collect(keys(r)))
	is_valid = [r[d][:status] for d in dates] .== "OPTIMAL"

	g_opt = reduce(hcat, [ v ? hcat(r[d][:g]...) : Array{Union{Missing, String}}(missing, n_gens, n_timesteps)  for (d, v) in zip(dates, is_valid)])
    E = co2_rates * g_opt

	return vec(E)
end

function get_deviations(metric, r1, r2)
    # Get all date-times in R1 and R2
    dts = intersect(keys(r1.data), keys(r2.data))

    # Compute deviations for all date-times in both
    mef = (ri, dt) -> get_nodal_mefs(ri.data, d -> d == dt, hybrid_mode=ri.hm)[:]
    @show size(mef(r1, dts[1]))
    @show size(mef(r2, dts[2]))

    return [metric(mef(r1, dt), mef(r2, dt)) for dt in dts]
end


"""
    function get_connected_clusters(x, cor_th)

Returns sets of features that are highly correlated with each other in a given
input matrix. 

`x` is a nxm matrix where n are observations and m are features. 
`cor_th` is the threshold for the correlation coefficient.
"""
function get_connected_clusters(x, cor_th)

    # XX is the correlation matrix
    n = size(x, 2) 
    XX = zeros(n, n)
    for i in 1:n, j in 1:n
        XX[i,j] = cor(x[:,i], x[:,j])
    end

    # construct a graph
    g = SimpleGraph(n);
	for i in 1:n, j in 1:n
        if XX[i,j]>=cor_th
            add_edge!(g, i, j);
        end
	end

    # extract the connected components
    return connected_components(g), XX
end