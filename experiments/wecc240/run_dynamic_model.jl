include("util.jl")

using Dates
using BSON
using CarbonNetworks
using TOML

config = TOML.parsefile("../../config.toml")
SAVE_DIR = config["data"]["SAVE_DIR"]

NUMBER_DAYS = 364
DATES = Date(2004, 01, 01) .+ Day.(0:NUMBER_DAYS)
HOURS = 1:24
DURATION = 24


function formulate_and_solve_dynamic(hour, day, month, T; Z=1e3, line_max=50_000.0, line_weight=1.3)
    case, _ = make_dynamic_case(hour, day, month, T)
    n, _ = size(case.A)
    # Construct flow matrix
    F = make_pfdf_matrix(case.A, case.β)

    # Get generator costs
    fl = [get_costs(case.heat, case.fuel, FUEL_COSTS) for _ in 1:T]
    fq = [zeros(length(fl[k])) for k in 1:length(fl)]

    # Get line capacities
    # (Specifically, increase them a little bit)
    pmax = [min.(line_weight * case.fmax, line_max) / Z for _ in 1:T]
    gmax = case.gmax / Z
    d = case.d / Z
    # TODO: probably need to divide C and P etc. by Z, too! Check. 

    # Formulate problem
    net = DynamicPowerNetwork(
        fq, fl, pmax, gmax, case.A, case.B, F, 
        case.P, case.C, T; η_c=case.η, η_d=case.η
        )
    pmp = DynamicPowerManagementProblem(net, d)

    # Solve
    solve!(pmp, CarbonNetworks.OPT)
    g = CarbonNetworks.evaluate(pmp.g)

    # Get generator emissions rates
    co2_rates = get_costs(case.heat, case.fuel, FUEL_EMISSIONS)

    # Compute MEFs
    mefs = zeros(n, T, T)
    λ = compute_mefs(pmp, net, d, co2_rates)
    for ind_t in 1:T
		mefs[:, :, ind_t] .= λ[ind_t];
	end

    @show (hour, day, month, pmp.problem.status)

    return (d=d, gmax=gmax, g=g, λ=mefs, status=pmp.problem.status)
end
case, meta = make_dynamic_case(1, 1, 1, DURATION)
results = [
    formulate_and_solve_dynamic(h, day(d), month(d), DURATION) for h in HOURS, d in DATES
    ]

bson(
    joinpath(SAVE_DIR, "wecc240_dynamic_results.bson")
    case=case, meta=meta, results=results
)