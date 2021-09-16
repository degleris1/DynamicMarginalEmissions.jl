using Convex, ECOS
using Random
using SparseArrays
using Zygote
using CarbonNetworks

using JLD

using Distributions: Uniform
using CarbonNetworks: compute_jacobian_kkt_dyn, _make_∇C

include("config.jl")

RESULTS_DIRECTORY = joinpath(@__DIR__, "../../results/")
FILE_PREFIX = "planning_"
FILE_SUFFIX = ".jld"

function load_results(savename)
    f = joinpath(RESULTS_DIRECTORY, FILE_PREFIX*savename*FILE_SUFFIX)
    r = JLD.load(f)

    θ, history = r["theta"], r["history"]
    config = Dict(zip(r["config"].keys, r["config"].values))
    P = create_planning_problem(config)

    return P, config, θ, history
end

function run_expansion_planning(config, savename)
    # Unpack
    init_seed = config[:initialization_seed]

    # Load problem
    P = create_planning_problem(config)
    
    # Initialize θ
    Random.seed!(init_seed)
    θ0 = P.θ_min .+ 3*rand(length(P.θ_min))

    # Run gradient descent
    θ, history = run_sgd(θ0, P; verbose=true)

    # Save results
    f = joinpath(RESULTS_DIRECTORY, FILE_PREFIX*savename*FILE_SUFFIX)
    JLD.save(f, "config", config, "theta", θ, "history", history)

    return config, θ, history
end


function create_planning_problem(config=DEFAULT_CONFIG)
    seed = config[:problem_seed]

    charge_efficiency = config[:charge_efficiency]
    discharge_efficiency = config[:discharge_efficiency]
    charge_rate = config[:charge_rate]
    emissions_rate = config[:emissions_rate]

    line_cost_per_mw_mile = config[:line_cost_per_mw_mile]
    line_length_min = config[:line_length_min]
    line_length_max = config[:line_length_max]
    storage_cost_per_mwh = config[:storage_cost_per_mwh]
    horizon = config[:horizon]
    emissions_weight = config[:emissions_weight]

    # Set up network
    net, d_dyn = create_network(config)
    η_c = charge_efficiency
    η_d = discharge_efficiency
    P0 = charge_rate
    c = emissions_rate

    # Set up planning parameters
    Random.seed!(seed)
    n, m = size(net.A)

    miles_per_line = rand(Uniform(line_length_min, line_length_max), m)
    θ_min = [copy(net.pmax); 1e-3 * ones(n)]
    Z = horizon / 1e3
    λ = emissions_weight
    ξ = [
        line_cost_per_mw_mile * miles_per_line;
        ones(n) * storage_cost_per_mwh
    ]

    return (
        net=net, d_dyn=d_dyn, η_c=η_c, η_d=η_d, P0=P0, c=c,
        θ_min=θ_min, ξ=ξ, λ=λ, Z=Z
    )
end


function create_network(config=DEFAULT_CONFIG)
    # Load config
    renewable_penetration = config[:renewable_penetration]
    demand_growth = config[:demand_growth]
    renewable_archetypes = config[:renewable_archetypes]
    emissions_tax = config[:emissions_tax]
    emissions_rate = config[:emissions_rate]
    seed = config[:net_seed]

    # Set seed
    Random.seed!(seed)

    # Load data
    net, d_peak, casefile = load_synthetic_network("case14.m")
    demand_data = load_demand_data("2021_07_01", normalize_rows=true)
    renew_data, _ = load_renewable_data("2021_07_01")
    renew_data ./= sum(renew_data, dims=1)

    # Parse dimensions
    n = length(d_peak)
    T, n_demand = size(demand_data)

    # Set demand / renewable profile for each node
    demand_profiles = rand(1:n_demand, n)
    renew_profiles = rand(renewable_archetypes, n)

    get_renew = (t, r) -> (r == 0 ? 0.0 : renew_data[t, r])

    # Set demand according to growth and profile
    d_dyn_no_renew = [
        demand_growth * d_peak .* demand_data[t, demand_profiles]
        for t in 1:T
    ]
    total_demands = [sum([d[i] for d in d_dyn_no_renew]) for i in 1:n]

    # Create renewable profiles
    renew_gmax = [
        get_renew.(t, renew_profiles) .* total_demands * renewable_penetration 
        for t in 1:T
    ]

    # Create net demand
    d_dyn = d_dyn_no_renew .- renew_gmax

    # Modify network
    net.gmax .= net.gmax * (sum(d_peak) * demand_growth / sum(net.gmax))
    net.fl ./= 100
    net.fq ./= 100

    # Add emissions tax
    net.fl .+= emissions_tax * emissions_rate

    return net, d_dyn
end


function run_sgd(θ0, P; num_iter=500, α=1e-5, verbose=false)
    # Unpack problem
    θ_min, ξ, c, λ, Z = P.θ_min, P.ξ, P.c, P.λ, P.Z

    # Initialize
    θ = deepcopy(θ0)
    loss_hist = []
    grad_hist = []
    θ_hist = []

    λ1 = 1
    λ2 = λ

    verbose && println("\n\n\n\n\nStarting planning problem...")
    for iter in 1:num_iter
        J, x, Dx, dnet, opf = formulate_and_solve_problem(θ, P)

        dJ = sum(_make_∇C(dnet, dnet.fl[1], dnet.fq[1], x), dims=2)[:, 1]
        dE = sum(_make_∇C(dnet, c), dims=2)[:, 1]
		
		if opf.problem.status in [Convex.MOI.ALMOST_OPTIMAL, Convex.MOI.OPTIMAL]
            # Compute objective
            E = sum([c'xt for xt in x])
    		O = ξ'θ + Z * (J + λ*E) - ξ'θ_min

            # Compute gradient
			∇J = -Dx.∂K_θT * (Dx.∂K_xT \ dJ)
			∇E = -Dx.∂K_θT * (Dx.∂K_xT \ dE)			
			dθ = λ1*ξ + Z * (λ1*∇J + λ2*∇E)

			# Update
			@. θ = θ - α*dθ
			θ .= max.(θ, θ_min)
	
			if iter > 1 && O < loss_hist[end]
				α = min(α * 1.05, Inf)
			elseif iter > 1
				α = 0.7 * α
			end
			
		else
			verbose && @warn (
                "Problem infeasible at iter $(iter)! Backtracking to last safe iterate."
            )
            
            # Compute objective
			O = Inf
            E = Inf
			
            # Compute gradient
			∇J = -2ξ
			∇E = zeros(length(θ))
			dθ = λ1*ξ + Z*(λ1*∇J + λ2*∇E)
			
			# Update
			@. θ = θ - α*dθ
			θ .= max.(θ, θ_min)
	
			α = α * 0.7
		end
		
		if verbose && iter % 100 == 0
			@show iter, O, J, E, α
			@show norm(dθ), norm(∇J), norm(∇E)
			println("---")
		end
			
		push!(loss_hist, O)
		push!(grad_hist, copy(dθ))
		push!(θ_hist, copy(θ))
	end

    return θ, (loss=loss_hist, grad=grad_hist, θ=θ_hist)
end


function formulate_and_solve_problem(θ, P)
    # Unpack
    net, d_dyn, η_c, η_d, P0 = P.net, P.d_dyn, P.η_c, P.η_d, P.P0
    n, m = size(net.A)
    pmax, C = θ[1:m], θ[m+1:end]


    T = length(d_dyn)
	
	# Construct problem
	dnet = make_dynamic(net, T, P0 * C, C, [net.gmax for _ in 1:T], η_c, η_d)
	dnet.pmax = [pmax for _ in 1:T]
	
	# Solve problem
	opf = DynamicPowerManagementProblem(dnet, d_dyn)
	solve!(opf, () -> ECOS.Optimizer(verbose=false))
	
	# Compute Jacobian
	x = flatten_variables_dyn(opf)
	
	# KKT operator wrt variables
	∂K_xT = adjoint(compute_jacobian_kkt_dyn(x, dnet, d_dyn))
	
	# KKT operator wrt params
	_, ∂K_pmaxT = Zygote.forward_jacobian(pmax -> kkt_dyn(x, dnet.fq, dnet.fl, d_dyn, [pmax for _ in 1:T], dnet.gmax, dnet.A, dnet.B, dnet.P, dnet.C, dnet.η_c, dnet.η_d), pmax)
	_, ∂K_CT = Zygote.forward_jacobian(C -> kkt_dyn(x, dnet.fq, dnet.fl, d_dyn, dnet.pmax, dnet.gmax, dnet.A, dnet.B, P0*C, C, dnet.η_c, dnet.η_d), C)
	
	∂K_θT = [∂K_pmaxT; ∂K_CT]
		
	return opf.problem.optval, evaluate.(opf.g), (∂K_xT=sparse(∂K_xT), ∂K_θT=∂K_θT), dnet, opf
end