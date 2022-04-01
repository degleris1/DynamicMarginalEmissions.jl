include("util.jl")

using Dates
using BSON
using CarbonNetworks

NUMBER_DAYS = 365
DATES = Date(2004, 01, 01) .+ Day.(0:(NUMBER_DAYS-1))
HOUR = 1
DURATION = 24

ECOS_OPT = CarbonNetworks.OPT

function formulate_and_solve_dynamic(
    hour, day, month, T; 
    Z=1e3, line_max=100e3, line_weight=2.0, OPT=ECOS_OPT
)
    println("-------")
    @time case, _ = make_dynamic_case(hour, day, month, T)
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
    P = case.P / Z
    C = case.C / Z
    ρ = case.ramp[1] / Z

    # Formulate problem
    net = DynamicPowerNetwork(
        fq, fl, pmax, gmax, case.A, case.B, F, case.S,
        P, C, T; η_c=case.η, η_d=case.η, ρ=ρ,
    )
    pmp = DynamicPowerManagementProblem(net, d)

    # Solve
    @time solve!(pmp, OPT)
    g = CarbonNetworks.evaluate(pmp.g)

    # Get generator emissions rates
    co2_rates = get_costs(case.heat, case.fuel, FUEL_EMISSIONS)

    # Compute MEFs
    mefs = zeros(n, T, T)
    @time λ = compute_mefs(pmp, net, d, co2_rates)
    for ind_t in 1:T
		mefs[:, :, ind_t] .= λ[ind_t];
	end

    @show (hour, day, month, pmp.problem.status)

    return (d=d, gmax=gmax, g=g, λ=mefs, co2_rates=co2_rates, status=pmp.problem.status)
end


case, meta = make_dynamic_case(1, 1, 1, DURATION)

results = []
for d in DATES
    try
        r = formulate_and_solve_dynamic(HOUR, day(d), month(d), DURATION)
        push!(results, r)
    catch
        @warn "No mefs on $d"
        push!(results, missing)
    end
end

bson(
    joinpath(SAVE_DIR, "wecc240_dynamic_results.bson"),
    case=case, meta=meta, results=results
)