using Random

using CarbonNetworks

include("config.jl")




function create_network(config=CONFIG, seed=1)
    Random.seed!(seed)

    renewable_penetration = config[:renewable_penetration]
    demand_growth = config[:demand_growth]
    renewable_archetypes = config[:renewable_archetypes]

    # Load data
    net, d_peak, casefile = load_synthetic_network("case14.m")

    demand_data = load_demand_data("2021_07_01", normalize_rows=true)

    renew_data, _ = load_renewable_data("2021_07_01")
    renew_data ./= sum(renew_data, dims=1)

    # Parse dimensions
    n = length(d_peak)
    T, n_demand = size(demand_data)

    # Set demand / renewable profile for each node
    demand_profiles = rand(1:n_demand, n)
    renew_profiles = rand(renewable_archetypes, n)

    get_renew = (t, r) -> (r == 0 ? 0.0 : renew_data[t, r])

    # Modify demand
    d_dyn_no_renew = [
        demand_growth * d_peak .* demand_data[t, demand_profiles]
        for t in 1:T
    ]
    total_demands = [sum([d[i] for d in d_dyn_no_renew]) for i in 1:n]

    # Create renewable generation
    renew_gmax = [
        get_renew.(t, renew_profiles) .* total_demands * renewable_penetration 
        for t in 1:T
    ]

    # Create net demand
    d_dyn = d_dyn_no_renew .- renew_gmax

    # Modify generation
    net.gmax .*= (demand_growth * (1 - renewable_penetration))

    return net, d_dyn
end