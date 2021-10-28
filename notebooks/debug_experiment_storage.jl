### A Pluto.jl notebook ###
# v0.16.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 2498bfac-3108-11ec-2b8b-7fb26f96afbb
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

# ╔═╡ 6a260a4f-9f96-464b-b101-1127e6ec48fe
begin
	using Revise
	using CarbonNetworks
end

# ╔═╡ 58ebe2b7-ea23-41fd-9cca-a3e02fdb4012
md"""
## Questions
- Why do the results barely change when I change the selection for `node`? Probably because the power flow constraints are not tight - can we make some of them tight? 
- Why do batteries SOC barely change when I increase `spen` even to 100%?

- What is the impact of having `cq = 0`, really? Is that important? How does that change the condition number of the jacobian? can it be reliable inverted?
"""

# ╔═╡ f94d2b5b-779a-4de0-9753-c077bc925fa1
begin
	ECOS_OPT = () -> ECOS.Optimizer(verbose=false)
	ECOS_OPT_2 = Convex.MOI.OptimizerWithAttributes(
		ECOS.Optimizer, "maxit"=> 100, "reltol"=>1e-6, "LogToConsole"=>false
		)
	OPT = ECOS_OPT
	δ = 1e-5;
end;

# ╔═╡ 1be301f4-31fe-44e9-895a-49bb1eec512f
theme(:default, label=nothing, 
		tickfont=(:Times, 8), guidefont=(:Times, 8), legendfont=(:Times, 8), titlefont=(:Times,8), framestyle=:box)

# ╔═╡ 5aef540d-719e-4212-9e2d-40c75cb685a7
LOAD_RANDOM_GRAPH = true

# ╔═╡ 06070be9-fae4-41b1-bdd3-066f7e785439
QUAD_COSTS = true

# ╔═╡ fb3e2e73-8b0b-43ba-a102-39ad9599941f
begin
	if LOAD_RANDOM_GRAPH
		n_ = 4
		l_ = 4
		A, B, cq, cl, _, gmax, pmax, _, _ = generate_random_data(n_, l_, 1)
		m_, _ = size(B)
		β = rand(Uniform(0, 1), m_)
		F = make_pfdf_matrix(A, β)
		
		if ~ QUAD_COSTS
			cq = [zeros(l_)]
		end
		net = PowerNetwork(cq[1], cl[1], pmax[1], gmax[1], A, B, F)
	else
		net, _, _ = load_synthetic_network("case14.m");
	end
	
end;

# ╔═╡ 4aaf967a-2c29-49c4-80f0-4d3311b8c1ff
net.fq

# ╔═╡ 0b806e4a-73a9-48f1-993c-875d0c9a01c3
pmax_ref = net.pmax;

# ╔═╡ f8072762-643c-485d-86d7-caf5a54c7d1e
n, m, l = get_problem_dims(net);

# ╔═╡ 058f7c12-d237-4130-aeb2-99279953d3f8
begin
	using LightGraphs, SimpleWeightedGraphs, GraphPlot
	Random.seed!(5)
	
	G = SimpleWeightedGraph(n)
	for j in 1:size(net.A, 2)
		inds = findall(x -> x != 0, net.A[:, j])
		add_edge!(G, inds[1], inds[2], net.pmax[j])
	end
	
	plt = gplot(G, nodelabel=1:n)
	plt
end

# ╔═╡ 73dac19a-3663-4f25-ac96-24e3a9997d3e
@bind T Slider(1:1:24)

# ╔═╡ 6b706efa-3343-4b9a-bd2f-1bf263707836
xticks_hr = 1:3:T;

# ╔═╡ 7d4a1f79-2a2f-4179-bdb6-ea394b7ca5fb
mef_times = 1:1:T;

# ╔═╡ 65b961b0-cc87-450e-bedb-43c5f38aafc5
md"""T = $T"""

# ╔═╡ b77196c7-987d-4175-9508-88c11cedbc3c
@bind spen Slider(0:0.1:3)

# ╔═╡ 2085073b-4770-4815-9e38-c36b850ab8d4
md"""
spen = $spen
"""

# ╔═╡ 01ec04b7-e317-43b9-86b8-8a7e9ecbcdea
@bind η Slider(0.8:0.05:1.0)

# ╔═╡ 69e3a691-b3a7-49c8-b178-01c82d9fa30f
md"""η = $η"""

# ╔═╡ 971d508d-c862-44ac-a0a4-a2fd9aa9156b
@bind node Slider(1:1:n)

# ╔═╡ c418c302-2bd5-40c5-aca6-e0302372903a
md"""node = $node"""

# ╔═╡ c9a5eeec-389b-4dbf-9c57-3ddd7ebff264
begin
α_max = 20
step = .1
αs = 1:step:α_max;
end;

# ╔═╡ c35434eb-c850-4650-870e-b4b719ff1a9b
@bind αi Slider(1:1:length(αs))

# ╔═╡ 6198b321-6573-4100-8166-6048c5fa2980
α = αs[αi];

# ╔═╡ 5c07c744-639b-4f56-af9a-a4b2366864bb
md"""α = $α"""

# ╔═╡ 377a4e9e-8d8c-4a76-9965-8df0a21dcf9c
begin
	pmax_mat = zeros(m, length(αs));
	for i in 1:length(αs)
		pmax_mat[:, i] = αs[i]*ones(m);
	end
	d_dyn = [rand(Bernoulli(0.8), n) .* rand(Exponential(2), n) for _ in 1:T];
	for i in 1:T
		d_dyn[i] = d_dyn[i]/sum(d_dyn[i]) * (0.75 * sum(net.gmax)) * rand(Uniform(0, 1));
	end
	emissions_rates = rand(Exponential(1), l)
	
end;

# ╔═╡ 141ec2da-8ab0-4fa4-91ea-01d88dbf5c42
@bind node_storage Slider(0:n)

# ╔═╡ c058136b-e536-41be-8c7a-28a8c51f3b22
md"""
Put the storage in a single location: $node_storage. 

If 0, then all the nodes are taken. 
"""

# ╔═╡ 6726672d-ab82-4364-bcea-bec33034bdac
begin
	
	mefs = zeros(
			n, T, length(mef_times)
		)
	println("Solving problem with parameters:")
	println("s = $spen, η = $η, α=$α")
	
	# Construct dynamic network
	C = sum(d_dyn) .* (spen + δ)/T
	if node_storage >= 1
		nodes_on = zeros(n)
		nodes_on[node_storage] = 1
		C = C.* nodes_on
	end
	P = 100*C 
	net.pmax = pmax_mat[:, αi];

	net_dyn = make_dynamic(net, T, P, C, η);

	# Construct and solve OPF problem
	opf_dyn = DynamicPowerManagementProblem(net_dyn, d_dyn)
	solve!(opf_dyn, OPT, verbose=false)

	if opf_dyn.problem.status != Convex.MOI.OPTIMAL
		@show opf_dyn.problem.status
	end

	# Compute MEFs
	mefs_ = compute_mefs(opf_dyn, net_dyn, d_dyn, emissions_rates)

	for ind_t in 1:T
		mefs[:, :, ind_t] .= mefs_[ind_t];
	end
	
	opf_dyn.problem.status
end

# ╔═╡ f59c420b-a1fe-4847-8115-a33166be54aa
begin
	s_vals = zeros(n, T+1)
	g_vals = zeros(l, T)
	p_vals = zeros(m, T)
	d_vals = zeros(n, T)
	for t in 1:T
    	s_vals[:, t+1] = evaluate(opf_dyn.s[t])./net_dyn.C
		g_vals[:, t] = evaluate(opf_dyn.g[t])./net_dyn.gmax[t]
		p_vals[:, t] = evaluate(opf_dyn.p[t])./net_dyn.pmax[t]
		d_vals[:, t] = d_dyn[t]
	end
end

# ╔═╡ 488db7c4-a048-47bd-b2bf-430d7f5664ae
let
	
	subplts = []
		
	
	
	plot_p = plot(abs.(p_vals)', xlabel="t", ylabel="p(t)/pmax", ylim=(-0.1, 1.1));
	push!(subplts, plot_p)
	
	t_axis = [t for t in 0:T]
	s_subplt = plot(t_axis, s_vals', ylim=(-.1, 1.1))
	title!("Relative SOC \n η = $η")
	xlabel!("Hour")
	ylabel!("SOC")
	
	push!(subplts, s_subplt)
	g_plt = plot(g_vals', xlabel="t", ylabel="g(t)/gmax", ylim=(-.1, 1.1))
	# plot!([sum(d_dyn[i]) for i in 1:T], ls=:dash)
	# plot!(sum(g_vals', 1), ls=:dash)
	push!(subplts, g_plt)
	
	g_s_plt = plot(
		[sum(d_vals[:, i]) for i=1:T], label="demand", lw=2, legend=:outertopright,
		xlabel="t", ylabel="E(t)"
	)
	# p_tot = [sum(p_vals[:, i]) for i=1:T]
	# plot(p_tot)
	g_tot = [sum(g_vals[:, i] .* net_dyn.gmax[i]) for i=1:T]
	plot!(g_tot, ls=:dash, label="g", lw = 4)
	# ds = [sum(evaluate(opf_dyn.s[1]))] 
	ds = vcat(
			sum(evaluate(opf_dyn.s[1])),
			[sum((evaluate(opf_dyn.s[i]) - evaluate(opf_dyn.s[i-1]))) for i=2:T]
			)

	plot!(ds, ls=:dash, label="ds", lw=4)
	plot!(g_tot - ds, ls=:dash, label="g-ds", lw=4)
	
	push!(subplts, g_s_plt)
	
	
	crt_mefs = mefs[node, :, :]
	lim = max(abs(maximum(crt_mefs)), abs(minimum(crt_mefs)));
	clims = (-lim, lim)
	subplt = heatmap(crt_mefs, 
		c=:balance, colorbar=true,
		xlabel="Consumption Time",
		title="$(100*spen)% storage, η=$η", 
		xticks=xticks_hr, yticks=xticks_hr, 
		clim=clims
	)
		
	plot!(ylabel="Emissions Time")
	# plot!(colorbar=true)
	push!(subplts, subplt)
	
	lay = @layout [Plots.grid(1, 3); Plots.grid(1, 2)]
	
	plt_emissions_heatmap = plot(subplts..., 
		layout=lay
		# size=(650, 200), 
		# bottom_margin=8Plots.pt
	)

	
	plt_emissions_heatmap
	
	
	
	
end

# ╔═╡ 2b8ccbdc-bbbe-4a0a-a71b-52fad873bc70
md""" node = $node"""

# ╔═╡ 79eb4f29-1457-41b4-96c3-9bfbe7a54208
begin
	
	subplts_nodes = []
	lim = max(abs(minimum(mefs)), abs(maximum(mefs)));
	clims = (-lim, lim)
	for node_ in 1:n
		crt_map = mefs[node_, :, :]
		subplt = heatmap(crt_map, 
		c=:bluesreds, colorbar=false,
		# xlabel="Consumption Time",
		title="node=$node_",
		# xticks=nothing, yticks=nothing,
		clim=clims,
	)
		# plot!(ylabel="Emissions Time")
	plot!(colorbar=true)
	push!(subplts_nodes, subplt)
	end
	
	
	
	plt_emissions_heatmap_nodes = plot(subplts_nodes..., 
		# layout=Plots.grid(5, 3, widths=[.29, 0.29, 0.42]), 
		# size=(650, 200), 
		# bottom_margin=8Plots.pt
	)

	
	plt_emissions_heatmap_nodes
	
	
	
	
end

# ╔═╡ bdd10c3a-5b1f-4f63-860c-268b9c4f035d
begin
	npoints = 20
	ε = 1e-3
	em_times = 1:1:T
end;

# ╔═╡ 124bc7c3-bb3b-4b36-bb7e-4f3f52b7de58
@bind cons_time Slider(1:1:T)

# ╔═╡ c6d3c265-d3e9-4f4e-b5e7-23a3762beb49
mefs[node, :, cons_time]

# ╔═╡ c1f12b9d-60fe-4c10-b74b-a63dfd965577
md"""Cons time = $cons_time"""

# ╔═╡ a133a850-22a8-4f3f-a2df-07c56735fbe3
begin
	println("Running sensitivity analysis")
	# size of the matrices are
	# 2npoints+1: number of different values of demand for which we solve the problem
	# n: number of nodes in the graph
	# l: number of generators (= length(emissions_rates))
	# T: the time horizon
	E_sensitivity = zeros(2npoints+1, length(emissions_rates), T);
	s_sensitivity = zeros(2npoints+1, n, T)
	g_sensitivity = zeros(2npoints+1, l, T);
	# println("initial value of the demand:")
	# println(d_dyn[cons_time][node])
	ref_val = deepcopy(d_dyn[cons_time][node])
	for i in -npoints:npoints
		d_crt = deepcopy(d_dyn)
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
		# println(d_dyn[cons_time][node])
		# println(d_crt[cons_time][node])
		# println(ref_val)
	end
end

# ╔═╡ c62e4fe2-2971-4479-af09-11b5467b2fee
@bind t_display Slider(1:1:T)

# ╔═╡ 5bd2052c-8f4b-4969-8864-8906ece79df5
md"""
Emissions time: $t_display
"""

# ╔═╡ 2aa296d4-7f47-4585-8a01-c70e9dddd2ac
let
	γ = 1e-4
	Δ = .02
	ylims = (1-Δ, 1+Δ)
	plt_s = plot(
		[1+i*ε for i in -npoints:npoints], 
		[s_sensitivity[:, k, t_display]/(s_sensitivity[npoints+1, k, t_display]+γ) for k in 1:n], ylim=ylims
	)
	title!("Storage at time $t_display")
	xlabel!("Change in demand at node $node at time $cons_time")
	ylabel!("Change in storage at all nodes at time $t_display")
	
	plt_E = plot(
		[1+i*ε for i in -npoints:npoints], 
		[E_sensitivity[:, k, t_display]./(E_sensitivity[npoints+1, k, t_display]+γ) for k in 1:length(emissions_rates)], ylim=ylims
		)
	title!("Emissions at time $t_display")
	xlabel!("Change in demand at node $node at time $cons_time")
	ylabel!("Change in emissions at all generators at time $t_display")
	
	plt_g = plot(
		[1+i*ε for i in -npoints:npoints], 
		[g_sensitivity[:, k, t_display]./(g_sensitivity[npoints+1, k, t_display]+γ) for k in 1:length(emissions_rates)], ylim=ylims
		)
	title!("Generators at time $t_display")
	xlabel!("Change in demand at node $node at time $cons_time")
	ylabel!("Change in generation at all generators at time $t_display")
	
	plt_E_tot = plot(
		[1+i*ε for i in -npoints:npoints], 
		sum(E_sensitivity[:, :, t_display], dims=2)./sum(E_sensitivity[npoints+1, :, t_display]), ylim=ylims
		)
	xlabel!("Change in demand at node $node at time $cons_time")
	ylabel!("Change in total emissions")
	
	#adding the theoretical curve for the sensitivity
	E_th = (
		sum(E_sensitivity[npoints+1, :, t_display]) .+ [ref_val * i * ε for i in -npoints:npoints] .* sum(mefs[node, t_display, cons_time])
		)./sum(E_sensitivity[npoints+1, :, t_display])
	plot!([1+i*ε for i in -npoints:npoints], E_th, ls=:dash)
	title!("Total emissions at time $t_display")
	
	plot([plt_s, plt_E, plt_g, plt_E_tot]..., size = (650, 650), lw = 3)

end

# ╔═╡ 9f8ef7ba-8e3d-4ac8-bb7f-44be3dee822d
md"""
---
"""

# ╔═╡ 2892cc44-a2a6-4e36-b9fc-a2d5e4e564bc
md"""
Fun GIFs with Anthony: 
"""

# ╔═╡ 0ed86b74-66e6-428e-84cc-5afe96e49060
@bind c_t Slider(1:1:T)

# ╔═╡ 1132387c-80b1-416d-b167-909be1e3e3ab
let
	foo = :Paired_12
	plt1 = plot()
	for node in 1:n
		plot!([sum(mefs[node, _t, :]) for _t in 1:T], palette=foo)
	end
	plot!()
	#bar([crt_map = results[node, _t, _t]' for node in 1:n])
	
	plt2 = plot()
	for node in 1:n
			plot!(mefs[node, c_t, :])
	end

	plot([plt1, plt2]...)
end

# ╔═╡ Cell order:
# ╟─58ebe2b7-ea23-41fd-9cca-a3e02fdb4012
# ╟─2498bfac-3108-11ec-2b8b-7fb26f96afbb
# ╟─6a260a4f-9f96-464b-b101-1127e6ec48fe
# ╟─f94d2b5b-779a-4de0-9753-c077bc925fa1
# ╟─1be301f4-31fe-44e9-895a-49bb1eec512f
# ╠═6b706efa-3343-4b9a-bd2f-1bf263707836
# ╠═7d4a1f79-2a2f-4179-bdb6-ea394b7ca5fb
# ╟─5aef540d-719e-4212-9e2d-40c75cb685a7
# ╠═06070be9-fae4-41b1-bdd3-066f7e785439
# ╠═fb3e2e73-8b0b-43ba-a102-39ad9599941f
# ╠═4aaf967a-2c29-49c4-80f0-4d3311b8c1ff
# ╟─058f7c12-d237-4130-aeb2-99279953d3f8
# ╟─0b806e4a-73a9-48f1-993c-875d0c9a01c3
# ╟─f8072762-643c-485d-86d7-caf5a54c7d1e
# ╟─65b961b0-cc87-450e-bedb-43c5f38aafc5
# ╟─73dac19a-3663-4f25-ac96-24e3a9997d3e
# ╟─2085073b-4770-4815-9e38-c36b850ab8d4
# ╟─b77196c7-987d-4175-9508-88c11cedbc3c
# ╟─69e3a691-b3a7-49c8-b178-01c82d9fa30f
# ╟─01ec04b7-e317-43b9-86b8-8a7e9ecbcdea
# ╟─c418c302-2bd5-40c5-aca6-e0302372903a
# ╟─971d508d-c862-44ac-a0a4-a2fd9aa9156b
# ╟─c9a5eeec-389b-4dbf-9c57-3ddd7ebff264
# ╟─5c07c744-639b-4f56-af9a-a4b2366864bb
# ╟─c35434eb-c850-4650-870e-b4b719ff1a9b
# ╟─6198b321-6573-4100-8166-6048c5fa2980
# ╟─377a4e9e-8d8c-4a76-9965-8df0a21dcf9c
# ╟─c058136b-e536-41be-8c7a-28a8c51f3b22
# ╟─141ec2da-8ab0-4fa4-91ea-01d88dbf5c42
# ╟─6726672d-ab82-4364-bcea-bec33034bdac
# ╟─488db7c4-a048-47bd-b2bf-430d7f5664ae
# ╟─f59c420b-a1fe-4847-8115-a33166be54aa
# ╟─2b8ccbdc-bbbe-4a0a-a71b-52fad873bc70
# ╠═c6d3c265-d3e9-4f4e-b5e7-23a3762beb49
# ╟─79eb4f29-1457-41b4-96c3-9bfbe7a54208
# ╟─bdd10c3a-5b1f-4f63-860c-268b9c4f035d
# ╟─c1f12b9d-60fe-4c10-b74b-a63dfd965577
# ╟─124bc7c3-bb3b-4b36-bb7e-4f3f52b7de58
# ╟─a133a850-22a8-4f3f-a2df-07c56735fbe3
# ╟─5bd2052c-8f4b-4969-8864-8906ece79df5
# ╟─c62e4fe2-2971-4479-af09-11b5467b2fee
# ╟─2aa296d4-7f47-4585-8a01-c70e9dddd2ac
# ╟─9f8ef7ba-8e3d-4ac8-bb7f-44be3dee822d
# ╟─2892cc44-a2a6-4e36-b9fc-a2d5e4e564bc
# ╟─0ed86b74-66e6-428e-84cc-5afe96e49060
# ╠═1132387c-80b1-416d-b167-909be1e3e3ab
