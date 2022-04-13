include("util.jl")

using Dates
using BSON
using CarbonNetworks

DATES = Date(2004, 01, 01) .+ Day.(0:364)
HOURS = 1:24


function formulate_and_solve_static(hour, day, month; Z=1e3, line_max=100.0, line_weight=2.0)
    case, _ = make_static_case(hour, day, month)

    # Construct flow matrix
    F = make_pfdf_matrix(case.A, case.β)

    # Get generator costs
    fl = get_costs(case.heat, case.fuel, FUEL_COSTS)
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

    # Get generator emissions rates
    co2_rates = get_costs(case.heat, case.fuel, FUEL_EMISSIONS)

    # Compute MEFs
    λ = compute_mefs(pmp, net, d, co2_rates)

    f_slack = [pmax; pmax] - abs.(F*(case.B * g - d))
    num_constr = sum(f_slack .< 1e-4)
    @show (hour, day, month, pmp.problem.status, num_constr)

    return (d=d, gmax=gmax, g=g, λ=λ, status=pmp.problem.status)
end

case, meta = make_static_case(1, 1, 1)
results = [formulate_and_solve_static(h, day(d), month(d)) for h in HOURS, d in DATES]

bson(
    joinpath(SAVE_DIR, "wecc240_static_results.bson"),
    case=case, meta=meta, results=results
)