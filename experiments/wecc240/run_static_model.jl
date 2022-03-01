include("util.jl")

using Dates
using BSON
using CarbonNetworks

DATES = Date(2004, 01, 01) .+ Day.(0:365)
HOURS = 1:24


function formulate_and_solve_static(hour, day, month; Z=1e3, line_max=50_000.0, line_weight=1.3)
    case, meta = make_static_case(hour, day, month)

    # Construct flow matrix
    F = make_pfdf_matrix(case.A, case.β)

    # Get generator costs
    fl = get_costs(case.heat, case.fuel, FUEL_COSTS)
    fq = zeros(length(fl))

    # Get line capacities
    # (Specifically, increase them a little bit)
    pmax = min.(line_weight * case.fmax, line_max) / Z
    gmax = case.gmax / Z
    d = case.d / Z

    # Formulate problem
    net = PowerNetwork(fq, fl, pmax, gmax, case.A, case.B, F)
    pmp = PowerManagementProblem(net, d)

    # Solve
    solve!(pmp, CarbonNetworks.OPT)
    g = CarbonNetworks.evaluate(pmp.g)

    # Get generator emissions rates
    co2_rates = get_costs(case.heat, case.fuel, FUEL_EMISSIONS)

    # Compute MEFs
    λ = compute_mefs(pmp, net, d, co2_rates)

    @show (hour, day, month, pmp.problem.status)

    return (d=d, gmax=gmax, g=g, λ=λ, status=pmp.problem.status)
end

case, meta = make_static_case(1, 1, 1)
results = [formulate_and_solve_static(h, day(d), month(d)) for h in HOURS, d in DATES]

bson(
    "/Users/degleris/Data/carbon_networks/wecc240_static_results.bson", 
    case=case, meta=meta, results=results
)