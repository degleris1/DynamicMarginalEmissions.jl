### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

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
	δ = 0.01
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

# ╔═╡ 30ec492c-8d21-43f6-bb09-32810494f21e
md"""
## How does storage penetration affect MEFs?
"""

# ╔═╡ 64b06f05-b1ab-44c4-80e2-faeeb1466421
c_rate = 0.25

# ╔═╡ 98a0d7c5-b1a8-4ebe-bb73-7ca88b475592
storage_penetrations = [0.0, 0.01, 0.02, 0.03, 0.04]

# ╔═╡ 8fc06205-0227-4b46-a2e9-72bdf9d57926
mef_times = 1:24

# ╔═╡ 008ed573-67ec-4908-a51a-c5d2a01e5b0e
refresh = false

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

# ╔═╡ 7e1bec31-9cdf-466c-83b3-dc792fd5cc53
md"""
### Utilities
"""

# ╔═╡ a3105ef4-67e7-41ad-bd2e-5d458c853d80
function make_dynamic(net, T, P, C, dyn_gmax)
	fqs = [net.fq for t in 1:T]
	fls = [net.fl for t in 1:T]
	pmaxs = [net.pmax for t in 1:T]
	return DynamicPowerNetwork(fqs, fls, pmaxs, dyn_gmax, net.A, net.B, P, C)
end

# ╔═╡ 6f08828b-4c4b-4f50-bd40-35805a37aae0
begin
	# Recompute results
	if refresh || !isfile("../results/storage.jld")
		options = (
			c_rate=c_rate,
			renewable_penetration=renewable_penetration,
			storage_penetrations=storage_penetrations,
			mef_times=mef_times,
			emissions_rates=emissions_rates,
			d_dyn=d_dyn,
		)
		
		meta = Dict()
		
		results = zeros(n, length(mef_times), length(storage_penetrations))
		for (ind_s, s_rel) in enumerate(storage_penetrations)

			# Construct dynamic network
			C = total_demands * s_rel .+ δ
			P = C * c_rate
			net_dyn = make_dynamic(net, T, P, C, dyn_gmax)

			# Construct and solve OPF problem
			opf_dyn = DynamicPowerManagementProblem(net_dyn, d_dyn)
			solve!(opf_dyn, OPT, verbose=false)
			@show opf_dyn.problem.optval / T

			# Compute MEFs
			@time mefs = compute_mefs(opf_dyn, net_dyn, d_dyn, emissions_rates)
			for ind_t in 1:length(mef_times)
				results[:, ind_t, ind_s] .= mefs[ind_t]
			end
			
			meta[ind_s] = (opf=opf_dyn, net=net_dyn)
		end
		
		println("")
		println("\007")
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

# ╔═╡ 6186798f-6711-4222-94bb-f53b2d0fad7d
begin
	subplots = []
	interesting_nodes = 1 : 2 : 29 
	curves = [1, 3, 5]
	
	for (ind_plt, i) in enumerate(interesting_nodes)
		plt = plot(xticks=[6, 12, 18, 24], xlim=(1, 24))
		plot!(ylim=(400, 1.1*maximum(emissions_rates)))
		plot!(legend=nothing)
		plot!(mef_times, results[i, :, curves], lw=3, alpha=0.8,
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
# ╟─30ec492c-8d21-43f6-bb09-32810494f21e
# ╠═64b06f05-b1ab-44c4-80e2-faeeb1466421
# ╠═98a0d7c5-b1a8-4ebe-bb73-7ca88b475592
# ╠═8fc06205-0227-4b46-a2e9-72bdf9d57926
# ╠═008ed573-67ec-4908-a51a-c5d2a01e5b0e
# ╠═6f08828b-4c4b-4f50-bd40-35805a37aae0
# ╠═ec676d81-d56e-452b-b739-3c3040bf6c8d
# ╠═cd5fe410-6264-4381-b19f-31d050bc3930
# ╟─f5b6479f-3bc8-4b86-987a-a14968d60e25
# ╠═d27ef0d8-70b2-4897-9000-8fa70b1862fc
# ╠═6186798f-6711-4222-94bb-f53b2d0fad7d
# ╟─418861e0-a35e-47db-9f00-a6a7fcf733fe
# ╟─0d07df73-89e1-4dbf-8f8c-82c202ad84c7
# ╠═43ab37f7-bf1f-44fb-8858-e3bf7d4e8880
# ╟─e1ef0db3-3130-45b0-9f07-c5776d72c31a
# ╟─08cce787-8118-4792-829a-153d2b637a78
# ╟─67ef9083-dbfe-48e8-a741-2a5fb035b8d7
# ╟─7e1bec31-9cdf-466c-83b3-dc792fd5cc53
# ╠═a3105ef4-67e7-41ad-bd2e-5d458c853d80
