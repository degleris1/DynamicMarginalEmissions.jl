using Convex
using CSV
using DataFrames
using Dates
using ECOS
using Gurobi
using JLD

using CarbonNetworks

include("config.jl")

PATH_CARBON = "NoEVs_year2019_dpdf_20210816.csv"
PATH_GENERATOR = "NoEVs_year2019_bsdf_20210816.csv"

RESULTS_DIRECTORY = joinpath(@__DIR__, "../../results/")
FILE_PREFIX_STATIC = "rodm_static_"
FILE_PREFIX_DYNAMIC = "rodm_dynamic_"
FILE_SUFFIX = ".jld"

DATE_FORMAT = DateFormat("yyyy-mm-dd HH:MM:SS")

HOURS_PER_DAY = 24
DAYS_PER_WEEK = 7
HOURS_PER_WEEK = HOURS_PER_DAY*DAYS_PER_WEEK

WEEKS = 1:52
DAYS = 1 : (maximum(WEEKS)*DAYS_PER_WEEK)
TIMES = 1 : (maximum(DAYS)*HOURS_PER_DAY)

BLOCK_LENGTH = 1
HOURS_PER_BLOCK = BLOCK_LENGTH*HOURS_PER_DAY
BLOCKS = 1 : Int(maximum(DAYS) / BLOCK_LENGTH)

CHEMICALS = (:co2, :so2, :nox)

GUROBI_ENV = Gurobi.Env()
GUROBI_OPT = Convex.MOI.OptimizerWithAttributes(() -> Gurobi.Optimizer(GUROBI_ENV), "LogToConsole" => false)

solve_ecos!(x) = solve!(x, () -> ECOS.Optimizer(verbose=false))

solve_gurobi!(x) = solve!(x, GUROBI_OPT)

results_path_static(f) = joinpath(RESULTS_DIRECTORY, FILE_PREFIX_STATIC*f*FILE_SUFFIX)

results_path_dynamic(f) = joinpath(RESULTS_DIRECTORY, FILE_PREFIX_DYNAMIC*f*FILE_SUFFIX)

"""
    load_results_dynamic(savename)

Load results from a dynamic dispatch model on the WECC dataset.
"""
function load_results_dynamic(savename)
    r = JLD.load(results_path_dynamic(savename))

    if !(typeof(r["config"]) <: Dict)
        config = Dict(zip(r["config"].keys, r["config"].values))
    else
        config = r["config"]
    end

    config[:datadir] = DEFAULT_CONFIG[:datadir]
    data, df = create_dispatch_time_series(config)
    dynamic_data = make_dynamic_data(data)

    return config, df, data, dynamic_data, r["results"]
end

"""
    load_results_static(savename)

Load results from a static dispatch model on the WECC dataset.
"""
function load_results_static(savename)
    r = JLD.load(results_path_static(savename))

    if !(typeof(r["config"]) <: Dict)
        config = Dict(zip(r["config"].keys, r["config"].values))
    else
        config = r["config"]
    end

    data, df = create_dispatch_time_series(config)

    return r["config"], df, data, r["results"]
end

"""
    run_rodm_dynamic(config, savename)

Run dynamic model on the WECC dataset, with parameters specified in `config`.
"""
function run_rodm_dynamic(config, savename)
    println("Loading data...")
    data, df = create_dispatch_time_series(config)
    dynamic_data = make_dynamic_data(data)

    println("Solving problems...")
    results = Vector{Any}(undef, length(dynamic_data))

    if config[:multithread]
        Threads.@threads for i in 1:length(dynamic_data)
            results[i] = formulate_and_solve_problem_dynamic(dynamic_data[i], config, df.gen.is_coal) 
        end
    else
        for i in 1:length(dynamic_data)
            results[i] = formulate_and_solve_problem_dynamic(dynamic_data[i], config, df.gen.is_coal) 
        end
    end

    println("Saving results...")
    f = results_path_dynamic(savename)
    JLD.save(f, "config", config, "results", results)

    println("Done!")
    return config, df, data, dynamic_data, results
end

"""
    run_rodm_static(config, savename)

Run static model on the WECC dataset, with parameters specified in `config`.
"""
function run_rodm_static(config, savename)
    println("Loading data...")
    data, df = create_dispatch_time_series(config)

    println("Solving problems...")
    results = formulate_and_solve_problem_static.(data)

    println("Saving results...")
    f = results_path_static(savename)
    JLD.save(f, "config", config, "results", results)

    println("Done!")
    return config, df, data, results
end

"""
    make_dynamic_data(data)

Take a sequence of individual snapshots (in `data`) and convert them into block snapshots
that specify dynamic problems.
"""
function make_dynamic_data(data)
    dynamic_data = [
        data[((w-1)*HOURS_PER_BLOCK+1) : (w*HOURS_PER_BLOCK)]
        for w in BLOCKS
    ]

    return dynamic_data
end

"""
    formulate_and_solve_problem_dynamic(Ps, config, is_coal)

Solve the dynamic problem with data `Ps` and configurgation `config`.
"""
function formulate_and_solve_problem_dynamic(Ps, config, is_coal)
    # Unpack
    C_rel = config[:storage_percentage]
    C_rate = config[:charge_rate]
    η_c = config[:charge_efficiency]
    η_d = config[:discharge_efficiency]
    ρ_rel = config[:coal_ramping_rate]
    gmin_coal = config[:coal_min]

    T = length(Ps)

    # Formulate and solve problem
    # We assume the network parameters do not change during the day
    θ = first(Ps).θ
    q = first(Ps).q
    d = [[P.d] for P in Ps]
    C = C_rel * sum(d)
    P = C_rate * C
    ρ = (ρ_rel*is_coal + 2*(1 .- is_coal)) .* θ.gmax
    
    dnet = make_dynamic(θ, T, P, C; η_c=η_c, η_d=η_d, ρ=ρ)
    
    # Solve unit commitment problem
    if config[:unit_commitment]
        opf, dnet, d, uc_active = make_unit_commitment!(dnet, d, is_coal, gmin_coal)
    else
        uc_active = 1:length(dnet.gmax[1])
        opf = DynamicPowerManagementProblem(dnet, d)
    end

    solve_gurobi!(opf)

    if !(opf.problem.status in (Convex.MOI.OPTIMAL, Convex.MOI.ALMOST_OPTIMAL))
        @warn opf.problem.status
    end

    g = evaluate.(opf.g)

    # Reduce problem
    opfr, netr, rel = reduce_problem(opf, dnet)

    # Compute generations and marginal emissions rates
    mefs = [compute_mefs(opfr, netr, d, getfield(q, chem)[uc_active][rel]) for chem in CHEMICALS]

    return (g=g, dq=mefs)
end

"""
    make_unit_commitment!(net, d, is_uc, gmin_uc_pc)

Convert a dynamic problem with network `net` and demand `d` into a unit committment problem.
"""
function make_unit_commitment!(net, d, is_uc, gmin_uc_pc)
    # Make UC selector matrix
    T = length(net.gmax)
    l = size(net.B, 2)
    n_uc = sum(is_uc)

    uc_inds = findall(x -> x == 1, is_uc)
    free_inds = findall(x -> x == 0, is_uc)

    M_uc = zeros(n_uc, l)
    for (i, gen_id) in enumerate(uc_inds)
        M_uc[i, gen_id] = 1
    end

    # Introduce binary variables
    z = Variable(n_uc, BinVar)

    # New capacities
    gmax_uc = M_uc * net.gmax[1]
    gmin_uc = gmin_uc_pc * gmax_uc

    # # Update network
    # for t in 1:T
    #     uc_net.gmax[t][uc_inds] .-= gmin_uc
    # end

    # Update demand
    #uc_d = [d[t] - net.B*M_uc'*(z .* gmin_uc) for t in 1:T]

    # Solve unit commitment problem
    uc = DynamicPowerManagementProblem(net, d)
    for t in 1:T
        add_constraints!(uc.problem, [
            M_uc * uc.g[t] >= z .* gmin_uc,
            M_uc * uc.g[t] <= z .* gmax_uc,
        ])
    end
    solve_gurobi!(uc)
    @show uc.problem.status

    # Now reduce optimal power flow problem
    # Specifically, eliminate 'off' generators
    # Reduce capacities of 'on' generators by their minimum output
    # Reduce demand by the minimum outputs
    # This problem should have the same result as the unit commitment problem above

    # Update network --- eliminate off generators and lower capacities
    active_uc_gens = findall(x -> abs(x-1) < 1e-2, evaluate(z))
    active_inds = [free_inds; uc_inds[active_uc_gens]]

    new_net = deepcopy(net)
    new_net.B = net.B[:, active_inds]
    new_net.ρ = net.ρ[active_inds]
    new_net.fq = [fq[active_inds] for fq in net.fq]
    new_net.fl = [fl[active_inds] for fl in net.fl]
    for t in 1:T
        new_net.gmax[t] = net.gmax[t][active_inds]
        new_net.gmax[t][length(free_inds)+1:end] -= gmin_uc[active_uc_gens]
    end

    # Update demand --- decrease demand
    new_d = [
        d[t] - net.B * M_uc[active_uc_gens, :]' * gmin_uc[active_uc_gens]
        for t in 1:T
    ]

    opf = DynamicPowerManagementProblem(new_net, new_d)


    return opf, new_net, new_d, active_inds
end

"""
    reduce_problem(opf, net; tol=1e-4)

Given a solved dynamic OPF problem `opf`, eliminate inactive constraints and return the
(solved) reduced problem.
"""
function reduce_problem(opf, net; tol=1e-4)
    gmax = net.gmax

    T = length(net.gmax)
    l = length(net.gmax[1])

    # Find generators that are empty and at capacity at time t
    rel_gen = t -> evaluate(opf.g[t]) ./ gmax[t]
    empty_gens = [rel_gen(t) .< tol for t in 1:T]
    capped_gens = [rel_gen(t) .> (1-tol) for t in 1:T]

    # Find generators that are always active / always at capacity
    never_active = [prod(x -> x[i], empty_gens) for i in 1:l]
    always_capped = [prod(x -> x[i], capped_gens) for i in 1:l]

    # Keep relevant generators (that aren't always bound to the same constraint)
    relevant_gens = broadcast(!, (always_capped .| never_active))
	
    net_red = DynamicPowerNetwork(
        [net.fq[t][relevant_gens] for t in 1:T],
        [net.fl[t][relevant_gens] for t in 1:T],
        [net.pmax[t] for t in 1:T],
        [net.gmax[t][relevant_gens] for t in 1:T],
        net.A,
        net.B[:, relevant_gens],
        net.C,
        net.P,
        T,
        η_c=net.η_c,
        η_d=net.η_d,
        ρ=net.ρ[relevant_gens]
    )

    demand_red = [
        dt .- sum(net.gmax[t][always_capped]) 
        for (t, dt) in enumerate(opf.params.d)
    ]

    opf_red = DynamicPowerManagementProblem(net_red, demand_red)
    solve_ecos!(opf_red)

    return opf_red, net_red, relevant_gens
end

"""
    formulate_and_solve_problem_static(P)

Given a problem specification `P`, formulate and solve the static problem.
"""
function formulate_and_solve_problem_static(P)
    d, θ, q = P
    opf = PowerManagementProblem(θ, d)
    solve_ecos!(opf)

    g = evaluate(opf.g)
    mefs = compute_mefs(opf, θ, [d], hcat(q...))

    return (g=g, dq=mefs)
end

"""
    create_dispatch_time_series(config=DEFAULT_CONFIG)

Given a configuration file `config`, load the time series of demands, generation capacities,
costs, etc.
"""
function create_dispatch_time_series(config=DEFAULT_CONFIG)
    # Unpack
    carbon_tax = config[:carbon_tax]

    # Load datasets
    df_carbon, df_gen = load_dataset(config)

    # Extract problem data
    get_gmax = w -> df_gen[:, "mw"*w]
    get_emit = (w, chem) -> df_gen[:, chem*w]
    get_fl = w -> (
        df_gen[:, "fuel_price"*w] .* df_gen[:, "heat_rate"*w]
        + carbon_tax * get_emit(w, "co2")
        + df_gen.vom
    )

    n, m, l = 1, 1, nrow(df_gen)
    A = zeros(n, m)
    B = ones(n, l)
    pmax = ones(m)
    fq = zeros(l)

    # Get network, demand, and emissions for each time
    datetimes = DateTime.(df_carbon.datetime[TIMES], DATE_FORMAT)
    get_week = t -> string(((dayofyear(datetimes[t]) - 1) ÷ 7) + 1)

    demands = df_carbon.demand[TIMES]
    nets = [
        PowerNetwork(fq, get_fl(get_week(t)), pmax, get_gmax(get_week(t)), A, B)
        for t in TIMES
    ]
    emissions_rates = [
        NamedTuple(Dict(
            chem => get_emit(get_week(t), string(chem))
            for chem in CHEMICALS
        ))
        for t in TIMES
    ]

    return (
        data=[(d=d, θ=θ, q=q) for (d, θ, q) in zip(demands, nets, emissions_rates)], 
        df=(gen=df_gen, carbon=df_carbon),
    )
end

"""
    load_dataset(config=DEFAULT_CONFIG)

Load carbon and generator data
"""
function load_dataset(config=DEFAULT_CONFIG)
    # Functionals
    load_df = (dir, f) -> DataFrame(CSV.File(expanduser(joinpath(dir, f))))

    # Unpack
    datadir = config[:datadir]

    # Load carbon and generator data
    # First two generators are all zero, so we remove them
    # We also remove any rows with missing data
    df_carbon = load_df(datadir, PATH_CARBON)
    df_gen = dropmissing!(load_df(datadir, PATH_GENERATOR))[3:end, :]

    return df_carbon, df_gen
end
