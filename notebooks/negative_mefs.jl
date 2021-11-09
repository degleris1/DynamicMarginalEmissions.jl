### A Pluto.jl notebook ###
# v0.16.3

using Markdown
using InteractiveUtils

# ╔═╡ 67394e67-6efc-46ea-92b7-2abb7be21b62
begin
	import Pkg
	Pkg.activate();
	using Random, Distributions
	using Convex, ECOS, Gurobi
	using Plots
	using PlutoUI
	using JLD
	using LinearAlgebra
end;

# ╔═╡ 3230e7ee-88ca-4825-9b74-e5e0877a6971
begin
	using Revise
	using CarbonNetworks
end

# ╔═╡ 822006ee-3f48-11ec-0651-b14492d06fd4
using LightGraphs, SimpleWeightedGraphs, GraphPlot

# ╔═╡ baef8089-356e-4c54-a522-eaa6d544563b
begin
	ECOS_OPT = () -> ECOS.Optimizer(verbose=false)
	OPT = ECOS_OPT
	δ = 1e-5;
end;

# ╔═╡ c7c88fb0-5b29-4eda-a0eb-26f2edbd287d
theme(:default, label=nothing, 
		tickfont=(:Times, 8), guidefont=(:Times, 8), legendfont=(:Times, 8), titlefont=(:Times,8), framestyle=:box)

# ╔═╡ 3a45e70f-e0ef-4945-a882-f91b365b758f
begin
n = 3
l = 2
end;

# ╔═╡ 1dd4f028-4986-4f3e-a4a4-4f2559a5ed34
begin
	Random.seed!(2)
	A = [-1 1 0; 1 0 -1; 0 -1 1]
	B = [1 0; 0 0; 0 1]
	_, m = size(A)
	β = rand(Uniform(0, 1), m)
	F = make_pfdf_matrix(A, β)
	fq = [0, 0]
	fl = [1, 0]
	pmax = [1, 1e4, 1e4]
	gmax = [1e3, 1e3]


	net = PowerNetwork(fq, fl, pmax, gmax, A, B, F)

	
	emissions_rates = [1, 0]
	
	
end;

# ╔═╡ 9fea2ea2-17f5-4213-ba62-6cf618630503
size(B)

# ╔═╡ f7f63ea9-be3d-4652-a194-05aea3a5a4f2
net.A

# ╔═╡ cb3db331-5845-45e8-8eda-a2d9ee3ce97f
net.B

# ╔═╡ be0f0e8b-7d15-471d-96b5-09c4ed04be3d
net.F

# ╔═╡ 3815a528-37cb-42a0-b696-fcb5764265e6
T = 1

# ╔═╡ 1ecf9020-cc94-4743-a2c5-fb66cbfeca69
d = [[100.; 1.; 0.]]

# ╔═╡ eb0fbeea-4dd5-45cc-8017-6f553faa9ed8
mef_times = [1]

# ╔═╡ eed8f462-edc4-469a-81a6-d036d0e471b8
spen = 0

# ╔═╡ 9fcd00ea-f444-4143-9e32-cc07afe25bc8
begin
	Random.seed!(2)
	mefs = zeros(
			n, T, length(mef_times)
		)
	# println("Solving problem with parameters:")
	# println("s = $spen, η = $η, α=$α, renewable_pc = $(100*percent_renewable)%")
	
	# Construct dynamic network
	C = zeros(n) .+ δ
	P = 100*C 
	η = 1

	net_dyn = make_dynamic(net, T, P, C, η);


	# Construct and solve OPF problem
	opf_dyn = DynamicPowerManagementProblem(net_dyn, d)
	solve!(opf_dyn, OPT, verbose=false)

	if opf_dyn.problem.status != Convex.MOI.OPTIMAL
		@show opf_dyn.problem.status
	end

	# Compute MEFs
	mefs_ = compute_mefs(opf_dyn, net_dyn, d, emissions_rates)

	for ind_t in 1:T
		mefs[:, :, ind_t] .= mefs_[ind_t];
	end
	
	opf_dyn.problem.status
end

# ╔═╡ aa950896-4edc-4fa6-8076-f38322e154ea
evaluate(opf_dyn.g[1])

# ╔═╡ 293a54a9-b381-4361-a583-23b72f3ab7e9
d[1] - B*evaluate(opf_dyn.g[1])

# ╔═╡ 8c331431-0376-4776-be20-e783cae37f0f
F*(d[1] - B*evaluate(opf_dyn.g[1])) # equals to p

# ╔═╡ 89faceec-4484-4727-9408-b174c94f0582
mefs

# ╔═╡ 17e0bdf9-14f1-46e6-b6bc-107dedb61395
md"""## sensitivity analysis"""

# ╔═╡ 73d6193d-9db5-44a9-8482-a4040b057dbd
begin
	s_vals = zeros(n, T+1)
	g_vals = zeros(l, T)
	p_vals = zeros(m, T)
	d_vals = zeros(n, T)
	for t in 1:T
    	s_vals[:, t+1] = evaluate(opf_dyn.s[t])./net_dyn.C
		g_vals[:, t] = evaluate(opf_dyn.g[t])./net_dyn.gmax[t]
		p_vals[:, t] = evaluate(opf_dyn.p[t])./net_dyn.pmax[t]
		d_vals[:, t] = d[t]
	end
end

# ╔═╡ 03a77bce-9e11-49cf-9e4d-7bb1fbb01e0b
begin 
	ε = 1e-3
	npoints = 10
	cons_time = 1
	node = 2
end

# ╔═╡ 7a97bc79-449e-4d7d-aed4-9a372d7eb826
node

# ╔═╡ 2afe8d77-b77b-4a41-bf66-01bd82d2682e
d_crt = deepcopy(d)

# ╔═╡ dd855679-5a2e-4061-8709-dc9159d657b9
begin
	println("Running sensitivity analysis")
	# size of the matrices are
	# 2npoints+1: number of different values of demand for which we solve the problem
	# n: number of nodes in the graph
	# l: number of generators (= length(emissions_rates))
	# T: the time horizon

	ref_val = deepcopy(d[cons_time][node])
	if ref_val > 0 
			perturb_vals = [ref_val * (1+i*ε) for i in -npoints:npoints]
			x_axis_vals = [1+i*ε for i in -npoints:npoints]
			idx_ref = npoints+1
	else
			println("""Demand is zero!""")
			perturb_vals = [i*ε for i in 0:npoints]
			x_axis_vals = perturb_vals
			idx_ref = 1
	end
	
	L = length(perturb_vals)
	E_sensitivity = zeros(L, length(emissions_rates), T);
	s_sensitivity = zeros(L, n, T)
	g_sensitivity = zeros(L, l, T);
	
	for i in -npoints:npoints
		d_crt = deepcopy(d)
		d_crt[cons_time][node] = ref_val * (1+i*ε)

		opf_ = DynamicPowerManagementProblem(net_dyn, d_crt)
		solve!(opf_, OPT, verbose=false)
		if opf_.problem.status != Convex.MOI.OPTIMAL
			@show opf_.problem.status
		end

		for t in 1:T
			s_sensitivity[i+npoints+1, :, t] = evaluate(opf_.s[t])
			g_sensitivity[i+npoints+1, :, t] = evaluate(opf_.g[t])
			E_sensitivity[i+npoints+1, :, t] = evaluate(opf_.g[t]).*emissions_rates
		end

	end
end

# ╔═╡ cbcee548-7b9a-4898-8c78-7a966ae1203f
t_display = 1

# ╔═╡ 5b0761ec-7836-4998-8042-f07a80e2536c
length(x_axis_vals)

# ╔═╡ 191e2bb2-96a2-49be-8c8f-e252bb18f6ee
mefs

# ╔═╡ 4c6bdc82-c6d7-4d97-b643-a94118ba35d5
let
	println("----------------")
	println(" Making plots")
	γ = 1e-4
	Δ = .002
	ylims = (1-Δ, 1+Δ)
	plt_s = plot(
		x_axis_vals, 
		[s_sensitivity[:, k, t_display]/(s_sensitivity[idx_ref, k, t_display]+γ) for k in 1:n], ylim=ylims
	)
	title!("Storage at time $t_display")
	xlabel!("Change in demand at node $node at time $cons_time")
	ylabel!("Change in storage at all nodes at time $t_display")
	
	plt_E = plot(
		x_axis_vals, 
		[E_sensitivity[:, k, t_display]./(E_sensitivity[idx_ref, k, t_display]+γ) for k in 1:length(emissions_rates)], ylim=ylims
		)
	title!("Emissions at time $t_display")
	xlabel!("Change in demand at node $node at time $cons_time")
	ylabel!("Change in emissions at all generators at time $t_display")
	
	plt_g = plot(
		x_axis_vals, 
		[g_sensitivity[:, k, t_display]./(g_sensitivity[idx_ref, k, t_display]+γ) for k in 1:length(emissions_rates)], ylim=ylims
		)
	title!("Generators at time $t_display")
	xlabel!("Change in demand at node $node at time $cons_time")
	ylabel!("Change in generation at all generators at time $t_display")
	
	plt_E_tot = plot(
		x_axis_vals, 
		sum(E_sensitivity[:, :, t_display], dims=2)./sum(E_sensitivity[idx_ref, :, t_display]), ylim=ylims
		)
	xlabel!("Change in demand at node $node at time $cons_time")
	ylabel!("Change in total emissions")
	
	#adding the theoretical curve for the sensitivity
	E_th = (
		sum(E_sensitivity[idx_ref, :, t_display]) .+ (perturb_vals.-ref_val) .* mefs[node, t_display, cons_time]
		)./sum(E_sensitivity[idx_ref, :, t_display])
	plot!(x_axis_vals, E_th, ls=:dash)
	title!("Total emissions at time $t_display")

	
	plot([plt_s, plt_E, plt_g, plt_E_tot]..., size = (650, 650), lw = 3)
	
end

# ╔═╡ Cell order:
# ╠═67394e67-6efc-46ea-92b7-2abb7be21b62
# ╠═822006ee-3f48-11ec-0651-b14492d06fd4
# ╠═3230e7ee-88ca-4825-9b74-e5e0877a6971
# ╠═baef8089-356e-4c54-a522-eaa6d544563b
# ╠═c7c88fb0-5b29-4eda-a0eb-26f2edbd287d
# ╠═3a45e70f-e0ef-4945-a882-f91b365b758f
# ╠═9fea2ea2-17f5-4213-ba62-6cf618630503
# ╠═1dd4f028-4986-4f3e-a4a4-4f2559a5ed34
# ╠═f7f63ea9-be3d-4652-a194-05aea3a5a4f2
# ╠═cb3db331-5845-45e8-8eda-a2d9ee3ce97f
# ╠═be0f0e8b-7d15-471d-96b5-09c4ed04be3d
# ╠═aa950896-4edc-4fa6-8076-f38322e154ea
# ╠═293a54a9-b381-4361-a583-23b72f3ab7e9
# ╠═8c331431-0376-4776-be20-e783cae37f0f
# ╠═3815a528-37cb-42a0-b696-fcb5764265e6
# ╠═1ecf9020-cc94-4743-a2c5-fb66cbfeca69
# ╠═eb0fbeea-4dd5-45cc-8017-6f553faa9ed8
# ╠═eed8f462-edc4-469a-81a6-d036d0e471b8
# ╠═9fcd00ea-f444-4143-9e32-cc07afe25bc8
# ╠═89faceec-4484-4727-9408-b174c94f0582
# ╠═17e0bdf9-14f1-46e6-b6bc-107dedb61395
# ╠═73d6193d-9db5-44a9-8482-a4040b057dbd
# ╠═03a77bce-9e11-49cf-9e4d-7bb1fbb01e0b
# ╠═7a97bc79-449e-4d7d-aed4-9a372d7eb826
# ╠═2afe8d77-b77b-4a41-bf66-01bd82d2682e
# ╠═dd855679-5a2e-4061-8709-dc9159d657b9
# ╠═cbcee548-7b9a-4898-8c78-7a966ae1203f
# ╠═5b0761ec-7836-4998-8042-f07a80e2536c
# ╠═191e2bb2-96a2-49be-8c8f-e252bb18f6ee
# ╠═4c6bdc82-c6d7-4d97-b643-a94118ba35d5
