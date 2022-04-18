using Pkg; Pkg.activate(joinpath(@__DIR__, "../../dev"))

using CSV
using Dates
using DataFrames
using SparseArrays: spzeros
using StatsBase: mean
using TOML
using XLSX
using BSON

using CarbonNetworks

config = TOML.parsefile(joinpath(@__DIR__, "../../config.toml"))
DATA_DIR = joinpath(config["data"]["DATA_DIR"], "wecc240")
SAVE_DIR = config["data"]["SAVE_DIR"]

BRANCH_PATH = joinpath(DATA_DIR, "Branches-Table 1.csv")
DEMAND_PATH = joinpath(DATA_DIR, "Load & Gen Profiles-Table 1.csv")
DEMAND_PARTICIPATION_PATH = joinpath(DATA_DIR, "Load Participation Factors-Table 1.csv")
GEN_PATH = joinpath(DATA_DIR, "Generators-Table 1.csv")
HEAT_PATH = joinpath(DATA_DIR, "Generator Heat Rates-Table 1.csv")
STORAGE_PATH = joinpath(DATA_DIR, "Storage & DR-Table 1.csv")

HEAT_RATE_COAL = 10.3

# https://www.eia.gov/totalenergy/data/annual/showtext.php?t=ptb0303
FUEL_COSTS = Dict('G' => 7.91, 'C' => 1.41, 'R' => 1.0)

# https://www.eia.gov/electricity/annual/html/epa_07_04.html
FUEL_COSTS_18 = Dict('G' => 3.55, 'C' => 2.06, 'R' => 1.0)

# https://www.epa.gov/sites/default/files/2015-07/documents/emission-factors_2014.pdf
FUEL_EMISSIONS = Dict('G' => 53.0, 'C' => 97.0, 'R' => 0.0)

include("nrel.jl")




function formulate_and_solve_dynamic(date, T; Z=1e3, line_max=100.0, line_weight=1.5)
    println("-------")
    @time case = make_dynamic_case(date, T)
    n, _ = size(case.A)

    # Construct flow matrix
    F = make_pfdf_matrix(case.A, case.β)

    # Get generator costs
    fl = case.fl
    fq = [zeros(length(fl[k])) for k in 1:length(fl)]

    # Get line capacities
    # (Specifically, increase them a little bit)
    pmax = [line_weight * min.(case.fmax / Z, line_max) for _ in 1:T]
    gmax = case.gmax / Z
    d = case.d / Z
    P = case.P / Z
    C = case.C / Z
    ρ = case.ramp[1] / Z

    # Formulate problem
    net = DynamicPowerNetwork(
        fq, fl, pmax, gmax, case.A, case.B, F, case.S,
        P, C, T; η_c=case.η_c, η_d=case.η_d, ρ=ρ,
    )
    pmp = DynamicPowerManagementProblem(net, d)

    # Solve
    @time solve!(pmp, CarbonNetworks.OPT)
    g = CarbonNetworks.evaluate.(pmp.g)
    p = CarbonNetworks.evaluate.(pmp.p)

    # Get generator emissions rates
    co2_rates = case.co2_rates

    # Compute MEFs
    mefs = zeros(n, T, T)
    @time λ = compute_mefs(pmp, net, d, co2_rates)
    for ind_t in 1:T
        mefs[:, :, ind_t] .= λ[ind_t];
    end

    # Solve again for equivalent static problems
    g_static = []
    p_static = []
    λ_static = []
    for t in 1:T
        @show t
        gmax_t = min.(gmax[t], g[t] .+ ρ)

        net_t = PowerNetwork(fq[t], fl[t], pmax[t], gmax_t, case.A, case.B, F)
        d_t = d[t] + case.S*(CarbonNetworks.evaluate(pmp.ch[t]) - CarbonNetworks.evaluate(pmp.dis[t]))
        pmp_t = PowerManagementProblem(net_t, d_t)

        # Solve
        @time solve!(pmp_t, CarbonNetworks.OPT)
        g_t = CarbonNetworks.evaluate(pmp_t.g)
        p_t = CarbonNetworks.evaluate(pmp_t.p)
        @time λ_t = compute_mefs(pmp_t, net_t, d_t, co2_rates)
        push!(g_static, g_t)
        push!(p_static, p_t)
        push!(λ_static, λ_t)
    end

    f_slack = [pmax[t] - abs.(p[t]) for t in 1:T]
    num_constr = mean(map(fs -> sum(fs .< 1e-2), f_slack))
    @show (date, pmp.problem.status, num_constr)

    return (
        g=g, p=p, λ=mefs, d=d, gmax=gmax, pmax=pmax, status=string(pmp.problem.status),
        g_static=g_static, p_static=p_static, λ_static=λ_static,
    )
end

function formulate_and_solve_static(date; Z=1e3, line_max=100.0, line_weight=1.5)
    case = make_static_case(date)

    # Construct flow matrix
    F = make_pfdf_matrix(case.A, case.β)

    # Get generator costs
    fl = case.fl
    fq = zeros(length(fl))

    # Get line capacities
    # (Specifically, increase them a little bit)
    pmax = line_weight * min.(case.fmax / Z, line_max)
    gmax = case.gmax / Z
    d = case.d / Z

    # Formulate problem
    net = PowerNetwork(fq, fl, pmax, gmax, case.A, case.B, F)
    pmp = PowerManagementProblem(net, d)

    # Solve
    solve!(pmp, CarbonNetworks.OPT)
    g = CarbonNetworks.evaluate(pmp.g)
    p = CarbonNetworks.evaluate(pmp.p)

    # Get generator emissions rates
    co2_rates = case.co2_rates

    # Compute MEFs
    λ = compute_mefs(pmp, net, d, co2_rates)

    f_slack = pmax - abs.(p)
    num_constr = sum(f_slack .< 1e-4)
    @show (date, pmp.problem.status, num_constr)

    return (g=g, p=p, λ=λ, d=d, gmax=gmax, pmax=pmax, status=string(pmp.problem.status))
end


"""
    get_costs(heat, fuel)

Return the cost of electricity for each generator.
"""
function get_costs(heat, fuel, fuel_costs)
    k = length(heat)

    fuel_costs = [get(fuel_costs, f, 0.0) for f in fuel]

    return heat .* fuel_costs
end

"""
    make_static_case(hour, day, month, year=2004)

Return data for specifying a static case (no storage or ramping).
"""
function make_static_case(date)
    @assert year(date) in [2004, 2018]

    if year(date) == 2004
        return _make_static_case2004(date)
    else  # year == 2018
        return _make_static_case2018(date)
    end
end

function _make_static_case2018(date)
    params = get_nrel_data(date)

    (; name, lat, lon) = params.node
    (; gmin, B, ramp, cost, fuel) = params.gen
    (; A, β, fmax) = params.line
    d = params.dt
    gmax = params.gt

    fl = cost

    # To get the CO2 rate, divide the cost by the fuel cost to get the heat rate
    # Then multiple by the emissions rate
    _fuel = [get(NREL_FUEL_MAP, f, 'R') for f in fuel]
    fuel_cost = [FUEL_COSTS_18[f] for f in _fuel]
    fuel_emissions_rates = [FUEL_EMISSIONS[f] for f in _fuel]
    co2_rates = (fl ./ fuel_cost) .* fuel_emissions_rates

    case = (
        A=A, β=β, fmax=fmax, 
        B=B, gmin=gmin, gmax=gmax, ramp=ramp, fuel=fuel, fl=fl, co2_rates=co2_rates,
        d=d,
        params=params,
    )
    return case
end

function _make_static_case2004(date)
    df = load_wecc_240_dataset()

    # Network structure
    node_names, node_ids = get_node_info(df.branch)
    A, β, fmax, cf = get_network_structure(df.branch)

    # Demand data
    demand_map = get_demand_map(hour(date)+1, day(date), month(date), year(date), df.demand)
    d = make_demand_vector(demand_map, node_names, df.participation)

    # Generator data
    B, gmin, gmax, ramp, heat, fuel = get_generator_data(demand_map, node_ids, df.gen, df.heat)

    fl = get_costs(heat, fuel, FUEL_COSTS)
    co2_rates = get_costs(heat, fuel, FUEL_EMISSIONS)

    case = (
        A=A, β=β, fmax=fmax, cf=cf, d=d, 
        B=B, gmin=gmin, gmax=gmax, ramp=ramp, fuel=fuel, fl=fl, co2_rates=co2_rates,
        node_names=node_names, node_ids=node_ids,
    )
    return case
end


"""
    make_dynamic_case(hour, day, month, duration, year=2004)

Return data for specifying a dynamic case, for a given duration T in number of timesteps.
"""
function make_dynamic_case(date, T)
    @assert year(date) in [2004, 2018]

    if year(date) == 2004
        return _make_dynamic_case2004(date, T)
    else  # year == 2018
        return _make_dynamic_case2018(date, T)
    end
end

function _make_dynamic_case2018(date, T, δ=1e-4)
    params = [get_nrel_data(date + Hour(t)) for t in 0:(T-1)]

    (; name, lat, lon) = params[1].node
    (; B, cost, fuel) = params[1].gen
    (; A, β, fmax) = params[1].line
    (; S, ηd, ηc, P, C) = params[1].storage

    d = [p.dt for p in params]
    gmin = [p.gen.gmin for p in params]
    gmax = [p.gt for p in params]
    ramp = [p.gen.ramp for p in params]
    
    fl = [cost for _ in 1:T]


    # To get the CO2 rate, divide the cost by the fuel cost to get the heat rate
    # Then multiple by the emissions rate
    _fuel = [get(NREL_FUEL_MAP, f, 'R') for f in fuel]
    fuel_cost = [FUEL_COSTS_18[f] for f in _fuel]
    fuel_emissions_rates = [FUEL_EMISSIONS[f] for f in _fuel]
    co2_rates = (cost ./ fuel_cost) .* fuel_emissions_rates

    case = (
        A=A, β=β, fmax=fmax, 
        B=B, gmin=gmin, gmax=gmax, ramp=ramp, fuel=fuel, fl=fl, co2_rates=co2_rates,
        S=S, P=P, C=C, η_c=mean(ηc), η_d=mean(ηd),
        d=d,
        params=params,
    )
    return case
end

function _make_dynamic_case2004(date, T, δ=1e-4)
    df = load_wecc_240_dataset()

    # Network structure
    node_names, node_ids = get_node_info(df.branch)
    A, β, fmax, cf = get_network_structure(df.branch) 

    # Demand and Generator data
    # Needs to be extracted at every time step
    # TODO: make sure this is the case
    d_dyn = [] 
    gmin_dyn = []
    gmax_dyn = []
    ramp_dyn = []

    # TODO: look for a better way to keep variables outside of for loops
    B = nothing
    heat = nothing
    fuel = nothing

    for dt = 0:T-1
        demand_map = get_demand_map(hour(date)+1+dt, day(date), month(date), year(date), df.demand)
        d = make_demand_vector(demand_map, node_names, df.participation)
        push!(d_dyn, d)

        B, gmin, gmax, ramp, heat, fuel = get_generator_data(demand_map, node_ids, df.gen, df.heat)
        push!(gmin_dyn, gmin)
        push!(gmax_dyn, gmax)
        push!(ramp_dyn, ramp)
    end

    fl = [get_costs(heat, fuel, FUEL_COSTS) for _ in 1:T]
    co2_rates = get_costs(heat, fuel, FUEL_EMISSIONS)

    # Storage data
    # TODO: clarify the meaning of each parameter
    # TODO: make vector-valued efficiencies compatible with the code
    efficiency, s_capacity, s_rate, _ = get_storage_data(df.storage)
    S = get_storage_map(df.storage, node_ids)

    case = (
        A=A, β=β, fmax=fmax, cf=cf, d=d_dyn, 
        B=B, gmin=gmin_dyn, gmax=gmax_dyn, ramp=ramp_dyn, fuel=fuel, fl=fl, co2_rates=co2_rates,
        η_c=sqrt(mean(efficiency)), η_d=sqrt(mean(efficiency)), C=s_capacity, P=s_rate, S=S,
        node_names=node_names, node_ids=node_ids,
    )
    return case
end

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
    get_generator_data(demand_map, df_gen)

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
    make_demand_vector(demand_map, df_participation)

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
    get_storage_map(storage_nodes, nodes)

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

Return node names and IDs.
"""
function get_node_info(df_branch)
    nodes = sort(unique([df_branch.source; df_branch.sink]))
    node_ids = map(x -> parse(Int, match(r"\{....\}", x).match[2:end-1]), nodes)

    return nodes, node_ids
end

"""
    get_network_structure(df_branch)

Return incidence matrix, line susceptances, line capacities, and line hurdle rates.
"""
function get_network_structure(df_branch)
    node_names, node_ids = get_node_info(df_branch)
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