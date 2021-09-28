using Convex
using CSV
using DataFrames
using Dates
using ECOS
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
WEEKS = 1:52
DAYS = 1:(maximum(WEEKS)*7)
TIMES = 1:(maximum(DAYS)*HOURS_PER_DAY)

CHEMICALS = (:co2,)


solve_ecos!(x) = solve!(x, () -> ECOS.Optimizer(verbose=false))

results_path_static(f) = joinpath(RESULTS_DIRECTORY, FILE_PREFIX_STATIC*f*FILE_SUFFIX)

results_path_dynamic(f) = joinpath(RESULTS_DIRECTORY, FILE_PREFIX_DYNAMIC*f*FILE_SUFFIX)

function load_results_dynamic(savename)
    r = JLD.load(results_path_dynamic(savename))
    data, df = create_dispatch_time_series(r["config"])
    dynamic_data = make_dynamic_data(data)

    return r["config"], df, data, dynamic_data, r["results"]
end

function load_results_static(savename)
    r = JLD.load(results_path_static(savename))
    data, df = create_dispatch_time_series(r["config"])

    return r["config"], df, data, r["results"]
end

function run_rodm_dynamic(config, savename)
    println("Loading data...")
    data, df = create_dispatch_time_series(config)
    dynamic_data = make_dynamic_data(data)

    println("Solving problems...")
    results = [
        formulate_and_solve_problem_dynamic(P, config, df.gen.is_coal) 
        for P in dynamic_data
    ]

    println("Saving results...")
    f = results_path_dynamic(savename)
    JLD.save(f, "config", config, "results", results)

    println("Done!")
    return config, df, data, dynamic_data, results
end

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

function make_dynamic_data(data)
    dynamic_data = [
        data[((d-1)*HOURS_PER_DAY+1) : (d*HOURS_PER_DAY)]
        for d in DAYS
    ]

    return dynamic_data
end

function formulate_and_solve_problem_dynamic(Ps, config, is_coal)
    println("Formulating and solving problem...")

    # Unpack
    C_rel = config[:storage_percentage]
    C_rate = config[:charge_rate]
    η_c = config[:charge_efficiency]
    η_d = config[:discharge_efficiency]
    ρ_rel = config[:coal_ramping_rate]

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
    opf = DynamicPowerManagementProblem(dnet, d)
    @time solve_ecos!(opf)
    @show opf.problem.status

    # TODO Reduce problem
    # @warn "Problem not reduced. Computing full Jacobian."

    # Compute generations and marginal emissions rates
    g = evaluate(opf.g)
    @time mefs = [compute_mefs(opf, dnet, d, getfield(q, chem)) for chem in CHEMICALS]

    return (g=g, dq=mefs)
end

function formulate_and_solve_problem_static(P)
    d, θ, q = P
    opf = PowerManagementProblem(θ, d)
    solve_ecos!(opf)

    g = evaluate(opf.g)
    mefs = compute_mefs(opf, θ, [d], hcat(q...))

    return (g=g, dq=mefs)
end

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