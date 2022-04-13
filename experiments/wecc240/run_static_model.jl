include("util.jl")

using BSON
using CarbonNetworks

yr = (length(ARGS) > 0) ? parse(Int, ARGS[1]) : 2018

NUM_HOURS = 24 * 3
DATES = DateTime(yr, 01, 01, 00) .+ Hour.(0:(NUM_HOURS-1))

@show unique(day.(DATES[month.(DATES) .== 2]))

function formulate_and_solve_static(date; Z=1e3, line_max=100.0, line_weight=1.5)
    case, _ = make_static_case(date)

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

    f_slack = pmax - abs.(F*(case.B * g - d))
    num_constr = sum(f_slack .< 1e-4)
    @show (date, pmp.problem.status, num_constr)

    return (g=g, λ=λ, d=d, gmax=gmax, status=string(pmp.problem.status))
end

case, meta = make_static_case(DATES[1])
results = [formulate_and_solve_static(d) for d in DATES]

bson(
    joinpath(SAVE_DIR, "wecc240_static_results_$(yr).bson"),
    meta=meta, case=case, results=results
)