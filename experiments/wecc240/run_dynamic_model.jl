include("util.jl")

using BSON
using CarbonNetworks


yr = (length(ARGS) > 0) ? parse(Int, ARGS[1]) : 2018

DURATION = 24
NUM_HOURS = 24*365
DATES = DateTime(yr, 01, 01, 00) .+ Hour.(0:DURATION:(NUM_HOURS-DURATION))

ECOS_OPT = CarbonNetworks.OPT

function formulate_and_solve_dynamic(
    date, T; 
    Z=1e3, line_max=100e3, line_weight=2.0, OPT=ECOS_OPT
)
    println("-------")
    @time case, _ = make_dynamic_case(date, T)
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
    @time solve!(pmp, OPT)
    g = CarbonNetworks.evaluate(pmp.g)

    # Get generator emissions rates
    co2_rates = case.co2_rates

    # Compute MEFs
    mefs = zeros(n, T, T)
    @time λ = compute_mefs(pmp, net, d, co2_rates)
    for ind_t in 1:T
		mefs[:, :, ind_t] .= λ[ind_t];
	end

    @show (date, pmp.problem.status)

    return (g=g, λ=mefs, status=pmp.problem.status, case=case)
end


case, meta = make_dynamic_case(DATES[1], DURATION)

results = []
for d in DATES
    # try
        r = formulate_and_solve_dynamic(d, DURATION)
        push!(results, r)
    # catch
    #     @warn "No mefs on $d"
    #     push!(results, missing)
    # end
end

bson(
    joinpath(SAVE_DIR, "wecc240_dynamic_results_$(yr).bson"),
    meta=meta, results=results
)