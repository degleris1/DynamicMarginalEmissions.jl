### A Pluto.jl notebook ###
# v0.15.1

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

# ╔═╡ db59921e-e998-11eb-0307-e396d43191b5
begin
	using Pkg; Pkg.activate()
	using Random
	using Convex, ECOS
	using Plots
	using PlutoUI
	using JLD
	using LinearAlgebra
	
	OPT = () -> ECOS.Optimizer(verbose=false)
	δ = 1e-4
end;

# ╔═╡ 0aac9a3f-a477-4095-9be1-f4babe1e2803
begin
	using Revise
	using CarbonNetworks
end

# ╔═╡ 257a6f74-d3c3-42eb-8076-80d26cf164ca
theme(:default, label=nothing, 
		tickfont=(:Times, 8), guidefont=(:Times, 8), legendfont=(:Times, 8))

# ╔═╡ 9bd515d4-c7aa-4a3d-a4fb-28686290a134
md"""
## Generate data
"""

# ╔═╡ 75dfaefd-abec-47e2-acc3-c0ff3a01048e
md"""
### Network
"""

# ╔═╡ f999d732-14b3-4ac5-b803-3df7a96ef898
begin
	net, d_peak = load_synthetic_network("case30.m")
	n = length(d_peak)
	
	# Limit line flows and generation capacities
	net.pmax *= 0.75
	net.gmax *= 0.7
	
	# Remove lines
	net.pmax[[2, 3, 5, 11, 12, 20, 22, 26, 27, 30, 39]] .*= δ
	
	# Add generators for each renewable source
	net.B = [net.B I(n)]
	net.gmax = [net.gmax; zeros(n)]
	net.fq = [net.fq; zeros(n)]
	net.fl = [net.fl; zeros(n)]
	
	"Network loaded."
end

# ╔═╡ 23690382-3d30-46e3-b26a-a30875be78ec
begin
	using LightGraphs, SimpleWeightedGraphs, GraphPlot
	Random.seed!(5)
	
	G = SimpleWeightedGraph(n)
	for j in 1:size(net.A, 2)
		inds = findall(x -> x != 0, net.A[:, j])
		add_edge!(G, inds[1], inds[2], net.pmax[j])
	end
	
	plt = gplot(G)
	plt
end

# ╔═╡ 1bd72281-4a7f-44f4-974d-632e9d0aaf28
md"""
### Demand and renewable time series
"""

# ╔═╡ 0c786da1-7f44-40af-b6d6-e0d6db2242b2
demand_data = load_demand_data("2021_07_01", normalize_rows=true);

# ╔═╡ 5b80ca83-0719-437f-9e51-38f2bed02fb4
begin
	renew_data, renew_labels = load_renewable_data("2021_07_01"; normalize_rows=true)
	renew_data ./= sum(renew_data, dims=1)
end;

# ╔═╡ 34d4bd62-6be2-4089-8caa-1a8715bee433
md"""
### Map time series to network data
"""

# ╔═╡ b7476391-30b9-4817-babf-7c9078531ee7
T, n_demand = size(demand_data)

# ╔═╡ cfcba5ad-e516-4223-860e-b1f18a6449ba
begin
	plt1 = plot(demand_data[:, rand(1:n_demand, 5)], lw=2, palette=:Blues)
	plt2 = plot(renew_data, lw=2, palette=:Reds)
	
	plt_time_series = plot(plt1, plt2, layout=(2, 1), size=(600, 300))
	
	savefig(plt_time_series, "../img/storage_penetration_time_series.png")
	plt_time_series
end

# ╔═╡ c82ef027-740a-49b1-93d2-1554c411a896
renewable_penetration = 0.25

# ╔═╡ 0239e1da-caf5-4593-af1b-5d1e8d2f2b3e
begin
	Random.seed!(1)
	demand_profiles = rand(1:n_demand, n)
	
	d_dyn = [
		d_peak .* demand_data[t, demand_profiles]
		for t in 1:T
	]
	
	total_demands = [sum([d[i] for d in d_dyn]) for i in 1:n] 
end;

# ╔═╡ c75d996b-0847-4cc7-bc93-3d42db39f1bd
n_renew = size(renew_data, 2)

# ╔═╡ 522e95d5-15a8-47ec-a79c-c4cc17cf86fd
begin
	Random.seed!(2)
	renew_profiles = rand(0:n_renew, n)
	
	get_renew = (t, renew_profile) -> 
		renew_profile == 0 ? 0.0 : renew_data[t, renew_profile]
	
	renew_gmax = [get_renew.(t, renew_profiles) .* total_demands * renewable_penetration for t in 1:T]
	
	dyn_gmax = [deepcopy(net.gmax) for t in 1:T]
	for t in 1:T
		dyn_gmax[t][end-n+1:end] .= max.(δ, renew_gmax[t])
	end
end;

# ╔═╡ bfa4a8f9-cfcd-4e22-b2dc-751226f3a73c
begin
	net_d_dyn = hcat(d_dyn...)' - hcat(renew_gmax...)'
	Random.seed!(7)
	plot(size=(600, 200))
	plot!(net_d_dyn[:, rand(1:n, 10)], lw=2, c=repeat(renew_profiles, outer=T)')
end

# ╔═╡ a8ccbc8e-24e6-4214-a179-4edf3cf26dad
md"""
### Carbon emissions data
"""

# ╔═╡ e8ee5cbb-4afc-4737-b006-90071f6138cd
begin	
	# [:NG, :NG_NEW, :NUC, :PEAK, :COL, :COL], roughly
	emissions_rates = [0.8e3, 0.7e3, 0.9e3, 1.4e3, 1.8e3, 1.8e3]
	emissions_rates = [emissions_rates; zeros(n)]
end;

# ╔═╡ c8f71644-9371-443b-b976-1734cc7ae583
md"""
## What is the impact of charging/discharging efficiency on storage dynamics?
"""

# ╔═╡ 57a37eb0-a547-428c-9f8c-e5f3e30f5260
c_rate = 0.25

# ╔═╡ 491da4e6-03e6-49d3-907d-43ddcffdfb42
s_rel = .1

# ╔═╡ 9a45c280-d8a4-4849-8977-2dc763ed633b
@bind η_c Slider(0.0:0.1:1.0)

# ╔═╡ e55e0566-7bd6-4126-8f60-0d940f6d8111
@bind η_d Slider(0.0:0.1:1.0)

# ╔═╡ 9deb6bdf-54f4-42ee-855d-7e82cef6f4bb
begin
	# Construct dynamic network
	C = total_demands * (s_rel)
	P = C * c_rate
	net_dyn = make_dynamic(net, T, P, C, dyn_gmax, η_c, η_d)

	# Construct and solve OPF problem
	opf_dyn = DynamicPowerManagementProblem(net_dyn, d_dyn)
	solve!(opf_dyn, OPT, verbose=false)
end

# ╔═╡ af6d4ef0-f59f-42be-af36-7cf447478e4c
begin
	s_vals = zeros(n, T)
	for t in 1:T
    	s_vals[:, t] = opf_dyn.s[t].value./C
	end
end

# ╔═╡ dd9efed9-f6ed-43fb-93f2-4d21d1360091
begin
	plot(s_vals', ylim=(0, 1))
	title!("Relative SOC \n η_c = $η_c, η_d = $η_d")
end

# ╔═╡ 3f9eb091-059c-44a5-9b50-ae3cabe24060
md"""
## How does charging efficiency affect MEFs
"""

# ╔═╡ 4ee8d0aa-9a68-44e5-8b72-f3158c4ba7f8
md"""
TODO LUCAS: basically implement the same loop as below with storage penetration fixed and with varying η_c and η_d
"""

# ╔═╡ 4a8c7752-6de8-4ea7-aafe-702b17507185
storage_pen = .05

# ╔═╡ 6e6b15b1-7685-4a20-9d94-dd703caa2fe9
η_vals = [0.1, 0.5, 0.9]

# ╔═╡ dfa2a03b-6925-4be0-aeac-076c4cf25969
interesting_nodes = 2 : 2 : 30 

# ╔═╡ 30ec492c-8d21-43f6-bb09-32810494f21e
md"""
## How does storage penetration affect MEFs?
"""

# ╔═╡ 98a0d7c5-b1a8-4ebe-bb73-7ca88b475592
storage_penetrations = [0.0, 0.05, 0.10]

# ╔═╡ 8fc06205-0227-4b46-a2e9-72bdf9d57926
mef_times = 1:24

# ╔═╡ 008ed573-67ec-4908-a51a-c5d2a01e5b0e
refresh = true

# ╔═╡ 47e2e177-3701-471f-ae3c-38276ec69370
begin
	# Recompute results
	if refresh || !isfile("../results/storage_eta.jld")
		println("Recomputing results")
		options_η = (
			c_rate=c_rate,
			renewable_penetration=renewable_penetration,
			storage_penetrations=storage_pen,
			mef_times=mef_times,
			emissions_rates=emissions_rates,
			d_dyn=d_dyn,
			η_vals=η_vals
		)
		
		meta_η = Dict()
		
		results_η = zeros(n, T, length(mef_times), length(η_vals))
		for (ind_η, η) in enumerate(η_vals)

			# Construct dynamic network
			C = total_demands * (storage_pen + δ)
			P = C * c_rate
			net_dyn = make_dynamic(net, T, P, C, dyn_gmax, η, η)

			# Construct and solve OPF problem
			opf_dyn = DynamicPowerManagementProblem(net_dyn, d_dyn)
			solve!(opf_dyn, OPT, verbose=false)
			@show opf_dyn.problem.optval / T
			
			# Compute MEFs
			@time mefs = compute_mefs(opf_dyn, net_dyn, d_dyn, emissions_rates)
			for ind_t in 1:length(mef_times)
				results_η[:, :, ind_t, ind_η] .= mefs[ind_t]
			end
			
			meta_η[ind_η] = (opf=opf_dyn, net=net_dyn)
		end
		
		println("")

		JLD.save("../results/storage_eta.jld", 
			"results", results_η, "options", options_η)
		
	else
		
		jld_file_η = JLD.load("../results/storage_eta.jld")
		results_η, options_η = jld_file_η["results"], jld_file_η["options"]
		
	end
	
	"Results loaded."
end

# ╔═╡ 457c1959-94fa-4267-8645-3ed1409cd0a0
total_mefs_η = sum(results_η, dims=2)[:, 1, :, :];

# ╔═╡ a7e75e49-5031-4cc4-b96e-6227277ec3ba
begin
	subplots_η = []
	curves_η = 1:length(η_vals)
	
	for (ind_plt, i) in enumerate(interesting_nodes)
		plt = plot(xticks=[6, 12, 18, 24], xlim=(1, 24))
		plot!(legend=nothing)
		plot!(mef_times, total_mefs_η[i, :, curves_η], lw=4, alpha=0.8,
			labels=η_vals[curves_η]')
		
		ind_plt in [1, 6, 11] && plot!(ylabel="Δco2 (lbs) / Δmwh")
		ind_plt in 11:15 && plot!(xlabel="hour")
		# ind_plt in [2] && plot!(legend=:topleft)
		
		push!(subplots_η, plt)
	end
	
	plot(subplots_η..., layout=(3, 5), leftmargin=4Plots.mm, size=(800, 400))
end

# ╔═╡ d9617524-76c3-447d-9b94-0a690f83a7b9
begin
	highlighted_node_ = 5
	plt_dyn_mef_eta = plot(subplots_η[highlighted_node_])
	plot!(size=(600, 200), legend=:outertopright)
	plot!(title="node $(interesting_nodes[highlighted_node_])", bottommargin=3Plots.mm)
	plot!(ylabel="Δco2 (lbs) / Δmwh", xlabel="hour")
	
	# savefig(plt_dynamic_mef, 
	# 	"../img/storage_penetration_dynamic_mef_$(highlighted_node).png")
	plt_dyn_mef_eta
end

# ╔═╡ 3dce8c04-8a5c-41a7-b18a-be06caa628d3
begin
	#fix values to avoid rerunning the whole cell whenever playing with the sliders
	η_c_ = 1.0
	η_d_ = 1.0
end

# ╔═╡ 6f08828b-4c4b-4f50-bd40-35805a37aae0
begin
	# Recompute results
	if refresh || !isfile("../results/storage.jld")
		println("Recomputing results")
		options = (
			c_rate=c_rate,
			renewable_penetration=renewable_penetration,
			storage_penetrations=storage_penetrations,
			mef_times=mef_times,
			emissions_rates=emissions_rates,
			d_dyn=d_dyn,
			η_c=η_c_,
			η_d=η_d_
		)
		
		meta = Dict()
		
		results = zeros(n, T, length(mef_times), length(storage_penetrations))
		for (ind_s, s_rel) in enumerate(storage_penetrations)

			# Construct dynamic network
			C = total_demands * (s_rel + δ)
			P = C * c_rate
			net_dyn = make_dynamic(net, T, P, C, dyn_gmax, η_c_, η_d_)

			# Construct and solve OPF problem
			opf_dyn = DynamicPowerManagementProblem(net_dyn, d_dyn)
			solve!(opf_dyn, OPT, verbose=false)
			@show opf_dyn.problem.optval / T
			
			# Compute MEFs
			@time mefs = compute_mefs(opf_dyn, net_dyn, d_dyn, emissions_rates)
			for ind_t in 1:length(mef_times)
				results[:, :, ind_t, ind_s] .= mefs[ind_t]
			end
			
			meta[ind_s] = (opf=opf_dyn, net=net_dyn)
		end
		
		println("")

		JLD.save("../results/storage.jld", 
			"results", results, "options", options)
		
	else
		
		jld_file = JLD.load("../results/storage.jld")
		results, options = jld_file["results"], jld_file["options"]
		
	end
	
	"Results loaded."
end

# ╔═╡ ec676d81-d56e-452b-b739-3c3040bf6c8d
begin
	test_t = 18
	test_s = 1
	bar(evaluate(meta[test_s].opf.g[test_t]) ./ meta[test_s].net.gmax[test_t])
	plot!(size=(700, 150))
end

# ╔═╡ cd5fe410-6264-4381-b19f-31d050bc3930
begin
	total_emissions = []
	for s in 1:length(storage_penetrations)
		E = 0
		for t in 1:T
			E += evaluate(meta[s].opf.g[t])' * emissions_rates
		end
		push!(total_emissions, E)
	end
	
	plt_total_emissions = plot()
	plot!(storage_penetrations, total_emissions, lw=4)
	plot!(size=(650, 150), xlabel="storage capacity (% total demand)")
	plot!(ylabel="co2 emissions", xlim=(0, Inf))
	plot!(bottom_margin=5Plots.mm, left_margin=5Plots.mm)
	
	savefig(plt_total_emissions, "../img/storage_penetration_total_emissions.png")
	plt_total_emissions
end

# ╔═╡ f5b6479f-3bc8-4b86-987a-a14968d60e25
"PC Renewable: $(sum(evaluate(meta[test_s].opf.g[test_t])[7:end]) / sum(d_dyn[test_t]))"

# ╔═╡ 0740dc70-a532-4818-b09d-b3b8d60fa6ba
total_mefs = sum(results, dims=2)[:, 1, :, :];

# ╔═╡ 6186798f-6711-4222-94bb-f53b2d0fad7d
begin
	subplots = []
	curves = 1:length(storage_penetrations)
	
	for (ind_plt, i) in enumerate(interesting_nodes)
		plt = plot(xticks=[6, 12, 18, 24], xlim=(1, 24))
		plot!(legend=nothing)
		plot!(mef_times, total_mefs[i, :, curves], lw=4, alpha=0.8,
			labels=storage_penetrations[curves]')
		
		ind_plt in [1, 6, 11] && plot!(ylabel="Δco2 (lbs) / Δmwh")
		ind_plt in 11:15 && plot!(xlabel="hour")
		# ind_plt in [2] && plot!(legend=:topleft)
		
		push!(subplots, plt)
	end
	
	plot(subplots..., layout=(3, 5), leftmargin=4Plots.mm, size=(800, 400))
end

# ╔═╡ d27ef0d8-70b2-4897-9000-8fa70b1862fc
begin
	highlighted_node = 5
	plt_dynamic_mef = plot(subplots[highlighted_node])
	plot!(size=(600, 200), legend=:outertopright)
	plot!(title="node $(interesting_nodes[highlighted_node])", bottommargin=3Plots.mm)
	plot!(ylabel="Δco2 (lbs) / Δmwh", xlabel="hour")
	
	savefig(plt_dynamic_mef, 
		"../img/storage_penetration_dynamic_mef_$(highlighted_node).png")
	plt_dynamic_mef
end

# ╔═╡ f7e0d09c-40bf-4936-987a-a3bcadae5487
begin
	node = interesting_nodes[highlighted_node]
	
	heatmap_subplts = []
	for s_idx in 1:length(storage_penetrations)
		
		subplt = heatmap(results[node, :, :, s_idx]', 
			c=:grayC, clim=(0, 2000), colorbar=false,
			xlabel="consumption time",
			title="$(100*storage_penetrations[s_idx])% storage"
		)
		
		s_idx == 1 && plot!(ylabel="emissions time")
		s_idx == 3 && plot!(colorbar=true)
		
		push!(heatmap_subplts, subplt)
	end
	
	plt_emissions_heatmap = plot(heatmap_subplts..., 
		layout=Plots.grid(1, 3, widths=[0.3, 0.3, 0.4]), 
		size=(650, 200), 
		bottom_margin=8Plots.pt
	)
	savefig(plt_emissions_heatmap, 
		"../img/storage_penetration_emissions_heatmap.png")
	plt_emissions_heatmap
end

# ╔═╡ 418861e0-a35e-47db-9f00-a6a7fcf733fe
md"""
## Appendix
"""

# ╔═╡ 0d07df73-89e1-4dbf-8f8c-82c202ad84c7
md"""
### Solve static problem

Solve the static optimal power flow problem and display the LMPs.
"""

# ╔═╡ 43ab37f7-bf1f-44fb-8858-e3bf7d4e8880
begin
	for t in 1:T
		opf_t = PowerManagementProblem(net, d_dyn[t])
		solve!(opf_t, OPT, verbose=true)
		@assert opf_t.problem.status == Convex.MOI.OPTIMAL
	end
		
	opf = PowerManagementProblem(net, d_dyn[18])
	solve!(opf, OPT, verbose=true)
	
	opf.problem.status
end

# ╔═╡ e1ef0db3-3130-45b0-9f07-c5776d72c31a
begin
	λ = get_lmps(opf)
	bar(λ, size=(600, 200), ylim=(100, Inf), title="locational marginal prices")
end

# ╔═╡ 08cce787-8118-4792-829a-153d2b637a78
bar(abs.(evaluate(opf.p)) ./ (net.pmax), size=(600, 100), title="line flows")

# ╔═╡ 67ef9083-dbfe-48e8-a741-2a5fb035b8d7
bar(evaluate(opf.g)[1:6], size=(600, 100))

# ╔═╡ Cell order:
# ╠═db59921e-e998-11eb-0307-e396d43191b5
# ╠═0aac9a3f-a477-4095-9be1-f4babe1e2803
# ╠═257a6f74-d3c3-42eb-8076-80d26cf164ca
# ╟─9bd515d4-c7aa-4a3d-a4fb-28686290a134
# ╟─75dfaefd-abec-47e2-acc3-c0ff3a01048e
# ╠═f999d732-14b3-4ac5-b803-3df7a96ef898
# ╠═23690382-3d30-46e3-b26a-a30875be78ec
# ╟─1bd72281-4a7f-44f4-974d-632e9d0aaf28
# ╠═0c786da1-7f44-40af-b6d6-e0d6db2242b2
# ╠═5b80ca83-0719-437f-9e51-38f2bed02fb4
# ╠═cfcba5ad-e516-4223-860e-b1f18a6449ba
# ╟─34d4bd62-6be2-4089-8caa-1a8715bee433
# ╠═b7476391-30b9-4817-babf-7c9078531ee7
# ╠═c82ef027-740a-49b1-93d2-1554c411a896
# ╠═0239e1da-caf5-4593-af1b-5d1e8d2f2b3e
# ╠═c75d996b-0847-4cc7-bc93-3d42db39f1bd
# ╠═522e95d5-15a8-47ec-a79c-c4cc17cf86fd
# ╠═bfa4a8f9-cfcd-4e22-b2dc-751226f3a73c
# ╟─a8ccbc8e-24e6-4214-a179-4edf3cf26dad
# ╠═e8ee5cbb-4afc-4737-b006-90071f6138cd
# ╟─c8f71644-9371-443b-b976-1734cc7ae583
# ╠═57a37eb0-a547-428c-9f8c-e5f3e30f5260
# ╠═491da4e6-03e6-49d3-907d-43ddcffdfb42
# ╠═9a45c280-d8a4-4849-8977-2dc763ed633b
# ╠═e55e0566-7bd6-4126-8f60-0d940f6d8111
# ╟─9deb6bdf-54f4-42ee-855d-7e82cef6f4bb
# ╠═af6d4ef0-f59f-42be-af36-7cf447478e4c
# ╠═dd9efed9-f6ed-43fb-93f2-4d21d1360091
# ╟─3f9eb091-059c-44a5-9b50-ae3cabe24060
# ╟─4ee8d0aa-9a68-44e5-8b72-f3158c4ba7f8
# ╠═4a8c7752-6de8-4ea7-aafe-702b17507185
# ╠═6e6b15b1-7685-4a20-9d94-dd703caa2fe9
# ╠═47e2e177-3701-471f-ae3c-38276ec69370
# ╠═457c1959-94fa-4267-8645-3ed1409cd0a0
# ╠═dfa2a03b-6925-4be0-aeac-076c4cf25969
# ╠═a7e75e49-5031-4cc4-b96e-6227277ec3ba
# ╠═d9617524-76c3-447d-9b94-0a690f83a7b9
# ╟─30ec492c-8d21-43f6-bb09-32810494f21e
# ╠═98a0d7c5-b1a8-4ebe-bb73-7ca88b475592
# ╠═8fc06205-0227-4b46-a2e9-72bdf9d57926
# ╠═008ed573-67ec-4908-a51a-c5d2a01e5b0e
# ╠═3dce8c04-8a5c-41a7-b18a-be06caa628d3
# ╠═6f08828b-4c4b-4f50-bd40-35805a37aae0
# ╠═ec676d81-d56e-452b-b739-3c3040bf6c8d
# ╠═cd5fe410-6264-4381-b19f-31d050bc3930
# ╟─f5b6479f-3bc8-4b86-987a-a14968d60e25
# ╠═0740dc70-a532-4818-b09d-b3b8d60fa6ba
# ╠═6186798f-6711-4222-94bb-f53b2d0fad7d
# ╠═d27ef0d8-70b2-4897-9000-8fa70b1862fc
# ╠═f7e0d09c-40bf-4936-987a-a3bcadae5487
# ╟─418861e0-a35e-47db-9f00-a6a7fcf733fe
# ╟─0d07df73-89e1-4dbf-8f8c-82c202ad84c7
# ╠═43ab37f7-bf1f-44fb-8858-e3bf7d4e8880
# ╟─e1ef0db3-3130-45b0-9f07-c5776d72c31a
# ╟─08cce787-8118-4792-829a-153d2b637a78
# ╟─67ef9083-dbfe-48e8-a741-2a5fb035b8d7
