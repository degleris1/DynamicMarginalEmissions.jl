using Pkg; Pkg.activate(joinpath(@__DIR__, "../../dev"))

using CSV
using Dates
using DataFrames
using SparseArrays: spzeros
using StatsBase: mean
using TOML
using XLSX
using BSON
using LinearAlgebra

using CarbonNetworks


Z = 1e3
line_max = 100.0
line_weight = 2.
use_NREL_network = true #only applies to 2004 dataset

config = TOML.parsefile(joinpath(@__DIR__, "../../config.toml"))
DATA_DIR = joinpath(config["data"]["DATA_DIR"], "wecc240")
SAVE_DIR = config["data"]["SAVE_DIR"]


HEAT_RATE_COAL = 10.3

# https://www.eia.gov/totalenergy/data/annual/showtext.php?t=ptb0303
FUEL_COSTS = Dict('G' => 7.91, 'C' => 1.41, 'R' => 1.0)

# https://www.eia.gov/electricity/annual/html/epa_07_04.html
FUEL_COSTS_18 = Dict('G' => 3.55, 'C' => 2.06, 'R' => 1.0)

# https://www.epa.gov/sites/default/files/2015-07/documents/emission-factors_2014.pdf
FUEL_EMISSIONS = Dict('G' => 53.0, 'C' => 97.0, 'R' => 0.0)

include("dataset_2004.jl")
include("nrel.jl")




function formulate_and_solve_dynamic(date, T; Z=Z, line_max=line_max, line_weight=line_weight, use_NREL_network=use_NREL_network)
    println("-------")
    case = make_dynamic_case(date, T; use_NREL_network=use_NREL_network)
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

    f_slack = [pmax[t] - abs.(p[t]) for t in 1:T]
    num_constr = mean(map(fs -> sum(fs .< 1e-2), f_slack))
    @show (date, pmp.problem.status, num_constr)

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


    return (
        g=g, p=p, λ=mefs, d=d, gmax=gmax, pmax=pmax, status=string(pmp.problem.status),
        g_static=g_static, p_static=p_static, λ_static=λ_static,
    )
end

function formulate_and_solve_static(date; Z=Z, line_max=line_max, line_weight=line_weight, use_NREL_network=use_NREL_network)

    case = make_static_case(date; use_NREL_network=use_NREL_network)

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

function get_F_and_pmax()
    case = make_static_case(DateTime(2018, 01, 01, 00))
    F = make_pfdf_matrix(case.A, case.β)
    pmax = case.fmax
    return F, pmax
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
function make_static_case(date; use_NREL_network=false)
    @assert year(date) in [2004, 2018]

    if year(date) == 2004
        return _make_static_case2004(date; use_NREL_network=use_NREL_network)
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

function _make_static_case2004(date; use_NREL_network=false)
    df = load_wecc_240_dataset()

    # you need to keep the 2004 formatting to match with all the
    # different files (e.g. df.participation)
    node_names, node_ids, _ = get_node_info(df.branch)

    # Demand data
    demand_map = get_demand_map(hour(date)+1, day(date), month(date), year(date), df.demand)
    d = make_demand_vector(demand_map, node_names, df.participation)

    # Generator data
    B, gmin, gmax, ramp, heat, fuel = get_generator_data(demand_map, node_ids, df.gen, df.heat)

    # Network structure
    if use_NREL_network

        # get the map between datasets
        node_nrel = get_node_params()
        node_ids_nrel = node_nrel.name
        map_nodes = get_map_nodes(node_ids, node_ids_nrel)

        # map parameters
        d, B = map_nodes*d, map_nodes*B
        
        (; A, β, fmax) = get_line_params(node_nrel.id_map)
    else
        A, β, fmax, _ = get_network_structure(df.branch)
    end

    fl = get_costs(heat, fuel, FUEL_COSTS)
    co2_rates = get_costs(heat, fuel, FUEL_EMISSIONS)

    case = (
        A=A, β=β, fmax=fmax, d=d, 
        B=B, gmin=gmin, gmax=gmax, ramp=ramp, fuel=fuel, fl=fl, co2_rates=co2_rates,
        node_names=node_names, node_ids=node_ids,
    )
    return case
end


"""
    make_dynamic_case(hour, day, month, duration, year=2004)

Return data for specifying a dynamic case, for a given duration T in number of timesteps.
"""
function make_dynamic_case(date, T; use_NREL_network=false)
    @assert year(date) in [2004, 2018]

    if year(date) == 2004
        return _make_dynamic_case2004(date, T; use_NREL_network=use_NREL_network)
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

function _make_dynamic_case2004(date, T, δ=1e-4; use_NREL_network=false)
    df = load_wecc_240_dataset()

    # Network structure
    node_names, node_ids, _ = get_node_info(df.branch)

    if use_NREL_network

        # get the map between datasets
        node_nrel = get_node_params()
        node_ids_nrel = node_nrel.name
        map_nodes = get_map_nodes(node_ids, node_ids_nrel)

        (; A, β, fmax) = get_line_params(node_nrel.id_map)

    else
        map_nodes = I(length(node_ids))
        A, β, fmax, _ = get_network_structure(df.branch) 
    end

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

        B, gmin, gmax, ramp, heat, fuel = get_generator_data(demand_map, node_ids, df.gen, df.heat)

        d, B = map_nodes*d, map_nodes*B


        push!(d_dyn, d)
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
    S = map_nodes*get_storage_map(df.storage, node_ids)

    case = (
        A=A, β=β, fmax=fmax, d=d_dyn, 
        B=B, gmin=gmin_dyn, gmax=gmax_dyn, ramp=ramp_dyn, fuel=fuel, fl=fl, co2_rates=co2_rates,
        η_c=sqrt(mean(efficiency)), η_d=sqrt(mean(efficiency)), C=s_capacity, P=s_rate, S=S,
        node_names=node_names, node_ids=node_ids,
    )
    return case
end


"""
    get_map_nodes(node_ids_04, node_ids_nrel, nicknames_04)

Returns a matrix that maps the nodes from the 2004 to the 2018 dataset. 
"""
function get_map_nodes(node_ids_04, node_ids_nrel)
	map_nodes = zeros(length(node_ids_nrel), length(node_ids_04))
	for (k, id_crt) in enumerate(node_ids_nrel)
		idx = findall(node_ids_04.==id_crt)
		if length(idx)==1
			map_nodes[k, idx[1]] = 1
		end
	end
	return map_nodes
end

