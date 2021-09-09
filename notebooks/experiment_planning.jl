### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 9c1c6f4e-1004-11ec-3aea-818cf0adbf67
begin
	using Pkg; Pkg.activate("")
	
	using Distributions
	using Random	
	using Convex, ECOS
	using Zygote
end

# ╔═╡ 738ada66-448b-4bd0-b5c5-d7aea2cb842f
begin
	using Revise
	using CarbonNetworks
end

# ╔═╡ b14255d5-bdb0-4a5a-80ce-ab004eaa370e
using Plots

# ╔═╡ 82eec4ea-4780-43d1-81c4-d5a2b461149b
using SparseArrays

# ╔═╡ 53a5f459-87c3-41fe-9c4e-835cfb7d27a1
using CarbonNetworks: _make_∇C

# ╔═╡ c41f71b1-17fb-4a15-b2f4-885d96ab1bc4
using LinearAlgebra

# ╔═╡ 4be4fdac-474f-4d93-b984-335ae099bde7
theme(:default, label=nothing)

# ╔═╡ 809702b6-08fc-48aa-9605-248de37ecde8
md"""
## Load network
"""

# ╔═╡ 7bbad73e-f84d-4614-910b-a83112e6213b
demand_data = load_demand_data("2021_07_01", normalize_rows=true);

# ╔═╡ c54bfcd8-1beb-4e99-9ca7-bd48a8b6cf93
begin
	renew_data, renew_labels = load_renewable_data("2021_07_01"; normalize_rows=true)
	renew_data ./= sum(renew_data, dims=1)
end;

# ╔═╡ 423e9a6c-4e6d-4fd2-9703-9748736fa278
md"""
#### Map time series to network
"""

# ╔═╡ f2ce2fdd-b842-4e60-bac3-00bf80821d3c
renewable_penetration = 0.50

# ╔═╡ 12d0dcfd-80c6-429a-b2b8-852eff9481f9
demand_growth = 1.2

# ╔═╡ 482843af-3901-4a52-9b8f-6e85c3f94c28
begin
	net, d_peak, casefile = load_synthetic_network("case14.m")
	net.gmax .*= (1 + demand_growth/2)
end

# ╔═╡ 0db77332-1785-4326-a442-5c545587dfc2
casefile["bus"]["5"]

# ╔═╡ a3a173cc-3594-4379-8bd3-802e243ffa15
n = length(d_peak)

# ╔═╡ d61f3dbd-410a-4a5c-8176-62fed16f33c0
T, n_demand = size(demand_data)

# ╔═╡ 29c7cdf9-1494-46fa-a94e-e313c07cee4a
renewable_archetypes = [0, 1, 3, 5, 1, 3, 5, 2, 4]

# ╔═╡ 216b61ca-9b68-4cff-9786-75978e35bf24
begin
	Random.seed!(1)
	demand_profiles = rand(1:n_demand, n)
end;

# ╔═╡ 7a33afa0-96d0-47d1-a492-6c5f26e2fbd1
d_dyn_no_renew = [
	(1 + demand_growth) * d_peak .* demand_data[t, demand_profiles]
	for t in 1:T
];

# ╔═╡ 611d4989-3404-451c-8722-481dca7d8d90
total_demands = [sum([d[i] for d in d_dyn_no_renew]) for i in 1:n]

# ╔═╡ 4c8516f4-ceb6-4df7-8810-b095aacba324
begin
	Random.seed!(2)
	renew_profiles = rand(renewable_archetypes, n)
end;

# ╔═╡ 08ad9073-b2b0-4052-9a5d-fd60563f5f76
d_peak

# ╔═╡ 41a87c76-9e29-4284-ae69-91760db71c13
renew_profiles

# ╔═╡ 99c5f1f3-fe3c-48f0-b80d-c084da31888a
get_renew = 
	(t, renew_profile) -> renew_profile == 0 ? 0.0 : renew_data[t, renew_profile]

# ╔═╡ 14181899-b014-438b-beb9-b5c2aef64fc8
renew_gmax = [get_renew.(t, renew_profiles) .* total_demands * renewable_penetration for t in 1:T]

# ╔═╡ 1f4a0995-0f99-4460-8b3b-a5306df76aa5
d_dyn = d_dyn_no_renew .- renew_gmax

# ╔═╡ 22bb2f8d-77ac-4093-88c6-73dd17ba76b5
findall(x -> sum(x) < 0, d_dyn)

# ╔═╡ 140e655c-5c40-4ee0-9a6f-203c1a3ea284
let
	_k = 12
	
	bar(d_dyn_no_renew[_k])
	bar!(d_dyn[_k])
end

# ╔═╡ 8f0afe62-9e47-4f6c-b8bf-10ae80c46132
md"""
## Formulate planning problem

We seek to minimize the operation cost (scaled by 365) plus the investment cost of the new wires and storage.
"""

# ╔═╡ 9e721cd4-a588-45d6-a765-c40dd2c597e8
md"""
#### Wire parameters
"""

# ╔═╡ ed7b4437-4358-41df-a1b5-4f44adf0b996
pmax_min = copy(net.pmax) / 2

# ╔═╡ 25c3d688-d9b1-4f43-9c29-ced4bf252eab
m = length(pmax_min)

# ╔═╡ 683d7dfc-4392-4d84-aa72-c1f3bdc0a2d5
line_cost_per_mw_mile = 3  # in thousands of dollars

# ╔═╡ 502a8d13-5f99-466d-a8c4-a8be8cc2d9e4
miles_per_line = rand(Uniform(25, 35), m)

# ╔═╡ b0b2893e-74c1-4d11-830e-5b3ee9806d1c
η = line_cost_per_mw_mile * miles_per_line  # cost per mw transmission

# ╔═╡ 11a6aea1-0c20-45c8-9054-a7b1c0a3275d
md"""
#### Storage parameters
"""

# ╔═╡ e0250e6c-7836-475f-94d0-af4025e23570
C_min = ones(n) * 1e-3

# ╔═╡ b5e88f54-b7bd-4d0e-ace7-72fd3cb9a38a
charge_rate = 0.25

# ╔═╡ 014f3d8e-f4cc-4ed2-8c58-b3c389345ea7
storage_cost_per_mwh = 350  # in thousands

# ╔═╡ 5ee900ca-3c3f-4f0f-b776-2a5c5163303a
ξ = ones(n)*storage_cost_per_mwh

# ╔═╡ 2c70c0a4-c64a-454c-af77-77c70f2a988d
η_c = 0.8

# ╔═╡ d1c548f1-d939-443b-bab8-0065fa157583
η_d = 1.0

# ╔═╡ 041a1bd6-ffe2-43bd-901b-712856602d2d
md"""
#### Emissions parameters
"""

# ╔═╡ d1dcccdd-050a-4516-9193-5713015edd1c
c = [0.8e3, 0.7e3, 0.9e3, 1.4e3, 1.8e3]

# ╔═╡ 230ffe69-e7c9-4a5d-a8e5-3a77507f9c30
md"""
#### Operation problem scaling
"""

# ╔═╡ 7a37521d-ceea-4a8e-ada4-39bc02472d26
horizon = 5 * 365

# ╔═╡ b313fce9-d02d-474d-aea9-a3122c991dc6
production_cost_scaling = horizon / 1e3

# ╔═╡ 76b5b11b-e9a5-4e14-8726-a54230799767
md"""
## Solve expansion planning problem
"""

# ╔═╡ 29c4b385-d539-457a-8abe-f99b3e658d89
λ = 1.0  # emissions weight

# ╔═╡ 4782ec56-2a80-4ce6-89a9-7fafb8fc22df
α = 5e-4  # step size

# ╔═╡ 73840389-33f8-47b9-9098-798ddfc423b5
begin
	Random.seed!(1)
	C_init = 5*ones(n) + 2*rand(n)
end

# ╔═╡ 7122567d-0ceb-4681-bbcb-88b1eeefddb0
begin
	Random.seed!(2)
	pmax_init = 10*pmax_min + 2*rand(m)
end

# ╔═╡ 1829b639-2c04-4fc7-b13b-85702e8e9191
initialize(pmax_init, C_init) = deepcopy([pmax_init; C_init])

# ╔═╡ eb3bbd85-670b-4fad-a48d-d92e17d399d7
function formulate_and_solve_problem(θ)
	pmax, C = θ[1:m], θ[m+1:m+n]
	
	# Construct problem
	dnet = make_dynamic(
		net, 
		T, 
		charge_rate * C,
		C, 
		[net.gmax for _ in 1:T],
		η_c, 
		η_d
	)
	dnet.pmax = [pmax for _ in 1:T]
	
	# Solve problem
	opf = DynamicPowerManagementProblem(dnet, d_dyn)
	solve!(opf, () -> ECOS.Optimizer(verbose=false))
	
	@show opf.problem.status
	
	# Compute Jacobian
	x = flatten_variables_dyn(opf)
	_, ∂K_xT = Zygote.forward_jacobian(x -> kkt_dyn(x, dnet, d_dyn), x)
	
	_, ∂K_pmaxT = Zygote.forward_jacobian(pmax -> kkt_dyn(x, dnet.fq, dnet.fl, d_dyn, [pmax for _ in 1:T], dnet.gmax, dnet.A, dnet.B, dnet.P, dnet.C, dnet.η_c, dnet.η_d), pmax)
	_, ∂K_CT = Zygote.forward_jacobian(C -> kkt_dyn(x, dnet.fq, dnet.fl, d_dyn, dnet.pmax, dnet.gmax, dnet.A, dnet.B, charge_rate*C, C, dnet.η_c, dnet.η_d), C)
	
	∂K_θT = [∂K_pmaxT; ∂K_CT]
	
	
	# test_b = rand(size(∂K_xT, 2))
	# test_x = sparse(∂K_xT) \ test_b
	# @assert (norm(sparse(∂K_xT)*test_x - test_b) / length(test_x)) < 1e-4 
	
	return opf.problem.optval, evaluate.(opf.g), (∂K_xT=sparse(∂K_xT), ∂K_θT=∂K_θT), dnet, opf
end

# ╔═╡ 16513e3d-587d-4cdf-899b-2f876cfd4315
project(θ, pmax_min, C_min) = max.([pmax_min; C_min], θ)

# ╔═╡ 65082604-1a8a-48e8-9e87-4b9cf7ad8693
begin
	θ = initialize(pmax_init, C_init)
	println("\n\n\n\n\nStarting planning problem...")
	
	loss_hist = []
	grad_hist = []
	
	for iter in 1:3_000
		J, x, Dx, dnet, opf = formulate_and_solve_problem(θ)
		
		E = sum([c'xt for xt in x])
		
		O = [η; ξ]'θ + production_cost_scaling * (J + λ*E)
		@show iter, O
		
		if opf.problem.status in [Convex.MOI.ALMOST_OPTIMAL, Convex.MOI.OPTIMAL]
			∇J = sum(-Dx.∂K_θT * (Dx.∂K_xT \ _make_∇C(dnet, dnet.fl[1], dnet.fq[1], x)), dims=2)[:, 1]
			∇E = sum(-Dx.∂K_θT * (Dx.∂K_xT \ _make_∇C(dnet, c)), dims=2)[:, 1]
		else
			∇J = -[η; ξ] / production_cost_scaling - ones(length(θ))
			∇E = -[η; ξ] / production_cost_scaling - ones(length(θ))
		end
		
		dθ = [η; ξ] + production_cost_scaling * (∇J + λ*∇E)
		@show norm(dθ), norm(∇J), norm(∇E)
			

		θ .= θ - α*dθ
		θ .= project(θ, pmax_min, C_min)
		
		push!(loss_hist, O)
		push!(grad_hist, dθ)
		
		println()
	end
end

# ╔═╡ 1c3d0bc2-ac19-4d48-8273-513be139d9c1
let
	plt1 = plot(loss_hist, lw=4, title="cost", size=(600, 200))
	#plt2 = plot(norm.(grad_hist), lw=4, title="norm of gradient")
	
	#plot(plt1, plt2, layout=(2, 1))
	
	plot(map(x -> x[end], grad_hist)[1:100])
end

# ╔═╡ 3a58391a-3ec8-42a0-8bf7-c8b54e80a378
@show loss_hist[end]

# ╔═╡ ad5fc387-5b38-4d71-b9bb-2e50c25d62f6
[η; ξ]' * θ - [η; ξ]' * [pmax_min; C_min]

# ╔═╡ b7e89c1d-268c-4cfd-8ee8-3ad50e5266f9
@show θ[1:m] - pmax_min

# ╔═╡ 07f7e05b-e660-4b6f-95e0-1bcf6ab7d0c6
@show θ[m+1:end] - C_min

# ╔═╡ 3cf7c255-0e16-4b76-a151-ec44988dac34
let
	p1 = bar(θ[m+1:end], title="proposed storage capacity")
	p2 = bar(net.B*net.gmax, title="generation")
	p3 = bar(d_dyn[10], title="demand")
	
	plot(p2, p3, p1, layout=(3, 1))
end

# ╔═╡ 41cce822-ba68-4a8c-a0e8-8930ae53d3d8
let
	L = net.A * diagm(sqrt.(net.pmax)) * net.A'
	vals, vecs = eigen(Matrix(L))
	
	plot()
	for col in 1:size(net.A, 2)
		fr = findfirst(x -> x == -1, net.A[:, col])
		to = findfirst(x -> x == 1, net.A[:, col])
		plot!([vecs[fr, 2]; vecs[to, 2]], [vecs[fr, 3]; vecs[to, 3]], 
			c=:black, lw=3, alpha=0.5)
	end
	scatter!(vecs[:, 2], vecs[:, 3])
	annotate!(vecs[:, 2] .- 0.02, vecs[:, 3] .+ 0.03, 1:n)
	plot!()
	
	
	#vals
end

# ╔═╡ 285adcdc-4f6f-4ef8-ba4f-5c8cf994a9f8
net.pmax

# ╔═╡ d68c0ec3-cec4-4552-b5ea-4275fbdf9a8f
Matrix(net.A)

# ╔═╡ 5d4f8975-807d-4d9e-b734-badcb1db3597
diag(Matrix(net.A'net.A)

# ╔═╡ b4270990-e676-4591-8d0a-2dee91b8f793
md"""
### Parameter ideas
- Push it further---more renewables
- Look at cheap vs expensive storage options
- Have homogeneous storage options, eg a pumped hydro option
- * Change emissions rates to match gen properties
- * Increase the emissions weights!
- Consider a emissions focused vs moderate planning scheme


* Cheap vs expensive storage
* Fix storage efficiency parameter

### Plot ideas

- Show which nodes had storage improved
- Plot network, before and after (show flows before and after)
- Show how flows vary at different times of day
- 
"""

# ╔═╡ Cell order:
# ╠═9c1c6f4e-1004-11ec-3aea-818cf0adbf67
# ╠═738ada66-448b-4bd0-b5c5-d7aea2cb842f
# ╠═b14255d5-bdb0-4a5a-80ce-ab004eaa370e
# ╠═4be4fdac-474f-4d93-b984-335ae099bde7
# ╠═809702b6-08fc-48aa-9605-248de37ecde8
# ╠═0db77332-1785-4326-a442-5c545587dfc2
# ╠═482843af-3901-4a52-9b8f-6e85c3f94c28
# ╠═a3a173cc-3594-4379-8bd3-802e243ffa15
# ╠═7bbad73e-f84d-4614-910b-a83112e6213b
# ╠═c54bfcd8-1beb-4e99-9ca7-bd48a8b6cf93
# ╠═423e9a6c-4e6d-4fd2-9703-9748736fa278
# ╠═f2ce2fdd-b842-4e60-bac3-00bf80821d3c
# ╠═12d0dcfd-80c6-429a-b2b8-852eff9481f9
# ╠═d61f3dbd-410a-4a5c-8176-62fed16f33c0
# ╠═29c7cdf9-1494-46fa-a94e-e313c07cee4a
# ╠═216b61ca-9b68-4cff-9786-75978e35bf24
# ╠═7a33afa0-96d0-47d1-a492-6c5f26e2fbd1
# ╠═611d4989-3404-451c-8722-481dca7d8d90
# ╠═4c8516f4-ceb6-4df7-8810-b095aacba324
# ╠═08ad9073-b2b0-4052-9a5d-fd60563f5f76
# ╠═41a87c76-9e29-4284-ae69-91760db71c13
# ╠═99c5f1f3-fe3c-48f0-b80d-c084da31888a
# ╠═14181899-b014-438b-beb9-b5c2aef64fc8
# ╠═1f4a0995-0f99-4460-8b3b-a5306df76aa5
# ╠═22bb2f8d-77ac-4093-88c6-73dd17ba76b5
# ╠═140e655c-5c40-4ee0-9a6f-203c1a3ea284
# ╠═8f0afe62-9e47-4f6c-b8bf-10ae80c46132
# ╟─9e721cd4-a588-45d6-a765-c40dd2c597e8
# ╠═ed7b4437-4358-41df-a1b5-4f44adf0b996
# ╠═25c3d688-d9b1-4f43-9c29-ced4bf252eab
# ╠═683d7dfc-4392-4d84-aa72-c1f3bdc0a2d5
# ╠═502a8d13-5f99-466d-a8c4-a8be8cc2d9e4
# ╠═b0b2893e-74c1-4d11-830e-5b3ee9806d1c
# ╠═11a6aea1-0c20-45c8-9054-a7b1c0a3275d
# ╠═e0250e6c-7836-475f-94d0-af4025e23570
# ╠═b5e88f54-b7bd-4d0e-ace7-72fd3cb9a38a
# ╠═014f3d8e-f4cc-4ed2-8c58-b3c389345ea7
# ╠═5ee900ca-3c3f-4f0f-b776-2a5c5163303a
# ╠═2c70c0a4-c64a-454c-af77-77c70f2a988d
# ╠═d1c548f1-d939-443b-bab8-0065fa157583
# ╠═041a1bd6-ffe2-43bd-901b-712856602d2d
# ╠═d1dcccdd-050a-4516-9193-5713015edd1c
# ╠═230ffe69-e7c9-4a5d-a8e5-3a77507f9c30
# ╠═7a37521d-ceea-4a8e-ada4-39bc02472d26
# ╠═b313fce9-d02d-474d-aea9-a3122c991dc6
# ╠═76b5b11b-e9a5-4e14-8726-a54230799767
# ╠═29c4b385-d539-457a-8abe-f99b3e658d89
# ╠═4782ec56-2a80-4ce6-89a9-7fafb8fc22df
# ╠═73840389-33f8-47b9-9098-798ddfc423b5
# ╠═7122567d-0ceb-4681-bbcb-88b1eeefddb0
# ╟─1829b639-2c04-4fc7-b13b-85702e8e9191
# ╠═82eec4ea-4780-43d1-81c4-d5a2b461149b
# ╠═53a5f459-87c3-41fe-9c4e-835cfb7d27a1
# ╟─eb3bbd85-670b-4fad-a48d-d92e17d399d7
# ╟─16513e3d-587d-4cdf-899b-2f876cfd4315
# ╠═65082604-1a8a-48e8-9e87-4b9cf7ad8693
# ╠═1c3d0bc2-ac19-4d48-8273-513be139d9c1
# ╠═3a58391a-3ec8-42a0-8bf7-c8b54e80a378
# ╠═ad5fc387-5b38-4d71-b9bb-2e50c25d62f6
# ╠═b7e89c1d-268c-4cfd-8ee8-3ad50e5266f9
# ╠═07f7e05b-e660-4b6f-95e0-1bcf6ab7d0c6
# ╠═3cf7c255-0e16-4b76-a151-ec44988dac34
# ╠═c41f71b1-17fb-4a15-b2f4-885d96ab1bc4
# ╠═41cce822-ba68-4a8c-a0e8-8930ae53d3d8
# ╠═285adcdc-4f6f-4ef8-ba4f-5c8cf994a9f8
# ╠═d68c0ec3-cec4-4552-b5ea-4275fbdf9a8f
# ╠═5d4f8975-807d-4d9e-b734-badcb1db3597
# ╠═b4270990-e676-4591-8d0a-2dee91b8f793
