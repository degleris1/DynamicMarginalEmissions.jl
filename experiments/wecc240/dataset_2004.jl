BRANCH_PATH = joinpath(DATA_DIR, "Branches-Table 1.csv")
DEMAND_PATH = joinpath(DATA_DIR, "Load & Gen Profiles-Table 1.csv")
DEMAND_PARTICIPATION_PATH = joinpath(DATA_DIR, "Load Participation Factors-Table 1.csv")
GEN_PATH = joinpath(DATA_DIR, "Generators-Table 1.csv")
HEAT_PATH = joinpath(DATA_DIR, "Generator Heat Rates-Table 1.csv")
STORAGE_PATH = joinpath(DATA_DIR, "Storage & DR-Table 1.csv")


"""
    get_storage_data(df_storage)

Return storage info.

Notes
-----
s_rate: constraint on the rate of charge
s_ramp: constraint on the change of the rate of charge
"""
function get_storage_data(df_storage)
    efficiency = df_storage.efficiency / 100
    s_capacity = df_storage.capacity * 1000
    # A multiplicative factor seems to be
    # needed to avoid Union{Missing, Float64}
    # which generates errors down the line
    s_rate = df_storage.max_power * 1 
    s_ramp = df_storage.ramp * 60

    return efficiency, s_capacity, s_rate, s_ramp
end

"""
    get_demand_map(hour, day, month, year, df_demand)

Return the demand for each WECC region and output for each renewable generator during the
specified time period.
"""
function get_demand_map(hour, day, month, year, df_demand)
    df = filter(r -> r.Year == year && r.Month == month && r.Day == day && r.Period == hour, df_demand)
    d = Dict()

    for (c, v) in zip(names(df), eachcol(df))
        if occursin("future", c)  # Skip these
            continue
        end
        d[replace(c, " base case" => "")] = v[1]
    end

    return d
end

"""
    get_generator_data(demand_map, node_ids, df_gen, df_heat)

Return the generator-node map `B`, the minimum outputs `gmin`, the maximum outputs `gmax`, 
the ramp rates `ramp`, the heat rates `heat`, and the fuel types `fuel`. Note that not
all parameters will be used in the model (e.g., ramp rates will not be used in a static 
case).
"""
function get_generator_data(demand_map, node_ids, df_gen, df_heat)
    k = nrow(df_gen)
    n = length(node_ids)

    B = spzeros(n, k)
    gmax = zeros(k) 
    fuel_types = [' ' for _ in 1:k]
    heat_rates = zeros(k)

    for (j, row) in enumerate(eachrow(df_gen))
        # Compute generator matrix
        gen_code = match(r"\{[0-9]{4} [A-Z]{1,2}\}", row.generator)
        gen_node = parse(Int, gen_code.match[2:5])

        i = findfirst(==(gen_node), node_ids)
        B[i, j] = 1

        # Compute max capacity
        try
            gmax[j] = demand_map[row.generator]
        catch
            try
                gmax[j] = parse(Float64, row.rating)
            catch
                gmax[j] = row.capacity
            end
        end

        # Compute fuel type
        fuel_types[j] = gen_code.match[end-1]

        # Compute heat rate
        hr = filter(r -> r["Aggregated Generator"] == row.generator, df_heat)
        if nrow(hr) > 0
            heat_rates[j] = mean(skipmissing(hr[1, 3:end]))
        elseif fuel_types[j] == 'C'
            heat_rates[j] = HEAT_RATE_COAL
        end
    end

    # Compute min output
    gmin = coalesce.(df_gen.min_output, 0.0)

    # Compute ramp rates
    ramp = coalesce.(df_gen.ramp*60, 2*gmax)

    return B, gmin, gmax, ramp, heat_rates, fuel_types
end

"""
    make_demand_vector(demand_map, nodes, df_participation)

Construct a vector containing demand at each node.
"""
function make_demand_vector(demand_map, nodes, df_participation)
    dfp = df_participation
    d = zeros(length(nodes))
    
    for (i, node) in enumerate(nodes)
        r = findfirst(r -> r.Node == node, eachrow(dfp))
        d[i] = dfp[r, "Load Participation Factor"] * demand_map[dfp[r, "Region"]]
    end

    return d
end

"""
    get_storage_map(df_storage, nodes)

Constructs matrix that maps storage data to vector-valued quantities
compatible with the codebase.
"""
function get_storage_map(df_storage, nodes_ids)

    storage_nodes_nms = df_storage.resource
    storage_nodes_ids = map(x -> parse(Int, match(r"[0-9]{4}", x).match), storage_nodes_nms)

    S = zeros(length(nodes_ids), length(storage_nodes_ids))
    for (i, s_ids) in enumerate(storage_nodes_ids)
        idx = findall(x->x==s_ids, nodes_ids)
        if length(idx)!=1
            # TODO: handle that in the proper way (e.g. exceptions)
            println("Duplicate nodes - breaking")
            break
        end
        S[idx, i] .= 1
    end
    return S
end

"""
    get_node_info(df_branch)

Return node names and IDs and nicknames
"""
function get_node_info(df_branch)
    nodes = sort(unique([df_branch.source; df_branch.sink]))
    node_ids = map(x -> parse(Int, match(r"\{....\}", x).match[2:end-1]), nodes)
    nicknames = map(x -> strip(rsplit(x, " ", limit=3)[1]), nodes)

    return nodes, node_ids, nicknames
end

"""
    get_network_structure(df_branch)

Return incidence matrix, line susceptances, line capacities, and line hurdle rates.
"""
function get_network_structure(df_branch)
    node_names, _, _ = get_node_info(df_branch)
    n, m = length(node_names), nrow(df_branch)

    β = abs.(1 ./ df_branch.reactance)
    fmax = Float64.(df_branch.capacity)
    cf = df_branch.hurdle_rate

    A = spzeros(n, m)
    for (j, row) in enumerate(eachrow(df_branch))
        src = findfirst(==(row.source), node_names)
        snk = findfirst(==(row.sink), node_names)

        A[src, j] = -1
        A[snk, j] = 1
    end

    return (A=A, β=β, fmax=fmax, cf=cf)
end

"""
    load_wecc_240_dataset()

Load raw tables from the 240-bus WECC test system.
"""
function load_wecc_240_dataset()
    # Branch data
    branch = select!(DataFrame(CSV.File(BRANCH_PATH)), 
        "From Node" => :source, "To Node" => :sink, 
        "Reactance (p.u.)" => :reactance, "Resistance (p.u.)" => :resistance,
        "Max Flow (MW)" => :capacity,
        "Hurdle Rate (Wheeling Charge) (\$/MWh)" => :hurdle_rate,
    )
    branch[!, :hurdle_rate] = coalesce.(branch[:, :hurdle_rate], 0)

    # Demand data
    demand = DataFrame(CSV.File(DEMAND_PATH))

    # Demand participation data
    participation = select!(DataFrame(CSV.File(DEMAND_PARTICIPATION_PATH)), 1:3)

    # Generation
    gen = select!(DataFrame(CSV.File(GEN_PATH)),
        "Generator" => :generator,
        "Category" => :category,
        "Max Capacity (MW)" => :capacity,
        "Ramp Rate (MW/min.)" => :ramp,
        "Rating (MW)" => :rating,
        "Min Stable Level (MW)" => :min_output
    )
    filter!(r -> !(occursin("New", r.generator)), gen)

    heat = filter!(r -> r.Property == "Heat Rate Incr", DataFrame(CSV.File(HEAT_PATH)))

    # Storage
    storage = select!(DataFrame(CSV.File(STORAGE_PATH)),
        "Resource" => :resource,
        "Gen Rating (MW)" => :max_power,
        "Ramp Rate (MW/min.)" => :ramp,
        "Storage Efficiency (%)" => :efficiency,
        "Storage Volume (GWh)" => :capacity,
    )
    filter!(r -> occursin("{", r.resource), storage)

    return (branch=branch, demand=demand, participation=participation, gen=gen, heat=heat, storage=storage)
end