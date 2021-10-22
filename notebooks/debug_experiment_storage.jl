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
"""

# ╔═╡ f94d2b5b-779a-4de0-9753-c077bc925fa1
begin
	ECOS_OPT = () -> ECOS.Optimizer(verbose=false)
	ECOS_OPT_2 = Convex.MOI.OptimizerWithAttributes(
		ECOS.Optimizer, "maxit"=> 100, "reltol"=>1e-6, "LogToConsole"=>false
		)
	OPT = ECOS_OPT
	δ = 1e-5
end

# ╔═╡ 1be301f4-31fe-44e9-895a-49bb1eec512f
theme(:default, label=nothing, 
		tickfont=(:Times, 8), guidefont=(:Times, 8), legendfont=(:Times, 8), titlefont=(:Times,8), framestyle=:box)

# ╔═╡ fb3e2e73-8b0b-43ba-a102-39ad9599941f
net, d_peak, _ = load_synthetic_network("case14.m");

# ╔═╡ 0b806e4a-73a9-48f1-993c-875d0c9a01c3
pmax_ref = net.pmax;

# ╔═╡ 4edf0252-38d2-415a-b7bf-a67d4448d188
pmax_ref

# ╔═╡ f8072762-643c-485d-86d7-caf5a54c7d1e
n, m, l = get_problem_dims(net);

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
α_max = 5
step = .1
αs = 1:step:α_max;
end

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

# ╔═╡ 6726672d-ab82-4364-bcea-bec33034bdac
begin
	
	results = zeros(
			n, T, length(mef_times)
		)
	println("s = $spen, η = $η, α=$α")
	
	# Construct dynamic network
	C = sum(d_dyn) .* (spen + δ)/T
	P = 100*C 
	net.pmax = pmax_mat[:, αi];

	net_dyn = make_dynamic(net, T, P, C, η);

	# Construct and solve OPF problem|
	opf_dyn = DynamicPowerManagementProblem(net_dyn, d_dyn)
	solve!(opf_dyn, OPT, verbose=false)

	@show opf_dyn.problem.status

	# Compute MEFs
	mefs = compute_mefs(opf_dyn, net_dyn, d_dyn, emissions_rates)

	for ind_t in 1:T
		results[:, :, ind_t] .= mefs[ind_t];
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
begin
	
	subplts = []
		
	crt_results = results[node, :, :]'
	
	subplt = heatmap(max.(crt_results, δ)', 
		c=:Blues_9, colorbar=false,
		xlabel="Consumption Time",
		title="$(100*spen)% storage, η=$η", xticks=xticks_hr, yticks=xticks_hr
	)
		
	plot!(ylabel="Emissions Time")
	# plot!(colorbar=true)
	push!(subplts, subplt)
	
	t_axis = [t for t in 0:T]
	s_subplt = plot(t_axis, s_vals', ylim=(-.1, 1.1))
	title!("Relative SOC \n η = $η")
	xlabel!("Hour")
	ylabel!("SOC")
	
	push!(subplts, s_subplt)
	g_plt = plot(g_vals', xlabel="t", ylabel="g(t)/gmax")
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
	
	plot_p = plot(p_vals', xlabel="t", ylabel="p(t)/pmax");
	push!(subplts, plot_p)
	
	lay = @layout [grid(1, 3); grid(1, 2)]
	
	plt_emissions_heatmap = plot(subplts..., 
		layout=lay
		# size=(650, 200), 
		# bottom_margin=8Plots.pt
	)

	
	plt_emissions_heatmap
	
	
	
	
end

# ╔═╡ 79eb4f29-1457-41b4-96c3-9bfbe7a54208
begin
	
	subplts_nodes = []
	for node_ in 1:n
		crt_map = results[node_, :, :]'
		subplt = heatmap(crt_map, 
		c=:Blues_9, colorbar=false,
		# xlabel="Consumption Time",
		title="node=$node_",
		xticks=nothing, yticks=nothing,
		#clim=(1.0, 1.5),
	)
		# plot!(ylabel="Emissions Time")
	# plot!(colorbar=true)
	push!(subplts_nodes, subplt)
	end
	
	
	
	plt_emissions_heatmap_nodes = plot(subplts_nodes..., 
		# layout=Plots.grid(5, 3, widths=[.29, 0.29, 0.42]), 
		# size=(650, 200), 
		# bottom_margin=8Plots.pt
	)

	
	plt_emissions_heatmap_nodes
	
	
	
	
end

# ╔═╡ 0ed86b74-66e6-428e-84cc-5afe96e49060
@bind c_t Slider(1:1:T)

# ╔═╡ 1132387c-80b1-416d-b167-909be1e3e3ab
let
	foo = :Paired_12
	plt1 = plot()
	for node in 1:n
		plot!([sum(results[node, _t, :]) for _t in 1:T], palette=foo)
	end
	plot!()
	#bar([crt_map = results[node, _t, _t]' for node in 1:n])
	
	plt2 = plot()
	for node in 1:n
			plot!(results[node, c_t, :])
	end

	plot([plt1, plt2]...)
end

# ╔═╡ Cell order:
# ╠═58ebe2b7-ea23-41fd-9cca-a3e02fdb4012
# ╠═2498bfac-3108-11ec-2b8b-7fb26f96afbb
# ╠═6a260a4f-9f96-464b-b101-1127e6ec48fe
# ╠═f94d2b5b-779a-4de0-9753-c077bc925fa1
# ╠═1be301f4-31fe-44e9-895a-49bb1eec512f
# ╠═6b706efa-3343-4b9a-bd2f-1bf263707836
# ╠═7d4a1f79-2a2f-4179-bdb6-ea394b7ca5fb
# ╠═fb3e2e73-8b0b-43ba-a102-39ad9599941f
# ╠═0b806e4a-73a9-48f1-993c-875d0c9a01c3
# ╠═4edf0252-38d2-415a-b7bf-a67d4448d188
# ╠═f8072762-643c-485d-86d7-caf5a54c7d1e
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
# ╟─6726672d-ab82-4364-bcea-bec33034bdac
# ╠═488db7c4-a048-47bd-b2bf-430d7f5664ae
# ╟─f59c420b-a1fe-4847-8115-a33166be54aa
# ╟─79eb4f29-1457-41b4-96c3-9bfbe7a54208
# ╟─0ed86b74-66e6-428e-84cc-5afe96e49060
# ╟─1132387c-80b1-416d-b167-909be1e3e3ab
