### A Pluto.jl notebook ###
# v0.16.1

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
	import Pkg
	Pkg.activate();
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

# ╔═╡ 1a7af4e0-2608-422c-bafe-d200f30bc4f3
using Statistics

# ╔═╡ c39005df-61e0-4c08-8321-49cc5fe71ef3
md"""
## Description
"""

# ╔═╡ 5867a4eb-470a-4a8a-84e1-6f150de1dcde
md"""
## ToDo
- Integrate ramping constraints
- Add a curve with 0% storage in the emissions vs time for different values of η plot
- Discuss validity of results of delayed and EARLY (i.e. impact before) emissions
- Discuss why there is this cross pattern
"""

# ╔═╡ 44275f74-7e7c-48d5-80a0-0f24609ef327
md"""
## Loading
"""

# ╔═╡ 257a6f74-d3c3-42eb-8076-80d26cf164ca
theme(:default, label=nothing, 
		tickfont=(:Times, 8), guidefont=(:Times, 8), legendfont=(:Times, 8), titlefont=(:Times,8), framestyle=:box)

# ╔═╡ 113e61a9-3b21-48d0-9854-a2fcce904e8a
xticks_hr = [0, 6, 12, 18, 24]

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
	net, d_peak, _ = load_synthetic_network("case30.m")
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
	plt1 = plot(
		demand_data[:, rand(1:n_demand, 5)], lw=2, palette=:Blues, title="Demand Data", xlabel="Hour", ylabel="Energy [?]", xticks=xticks_hr)
	plt2 = plot(renew_data, lw=2, palette=:Reds, title="Renewable Data", xlabel="Hour", ylabel="Energy [?]", xticks=xticks_hr)
	
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

# ╔═╡ 9e5a1672-d452-42f5-ba5a-a2fa0b1eaada
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
c_rate = 0.25 #charging rate

# ╔═╡ 491da4e6-03e6-49d3-907d-43ddcffdfb42
s_rel = .1 #storage penetration

# ╔═╡ 147a8ae3-a910-4d76-9689-4e77f7012914
md"""
This cursor allows for the definition/choice of the charge efficiency
"""

# ╔═╡ 9a45c280-d8a4-4849-8977-2dc763ed633b
@bind η_c Slider(0.8:0.01:1.0)

# ╔═╡ 5bd157b8-3dac-464e-85b0-c63830a7e817
md"""
This cursor allows for the definition/choice of the discharge efficiency
"""

# ╔═╡ e55e0566-7bd6-4126-8f60-0d940f6d8111
@bind η_d Slider(0.8:0.01:1.0)

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
	s_vals = zeros(n, T+1)
	for t in 1:T
    	s_vals[:, t+1] = evaluate(opf_dyn.s[t])./C
	end
end

# ╔═╡ dd9efed9-f6ed-43fb-93f2-4d21d1360091
begin
	t_axis = [t for t in 0:T]
	plot(t_axis, s_vals', ylim=(0, 1))
	title!("Relative SOC \n η_c = $η_c, η_d = $η_d")
	xlabel!("Hour")
	ylabel!("SOC")
end

# ╔═╡ 3f9eb091-059c-44a5-9b50-ae3cabe24060
md"""
## How does charging efficiency affect MEFs
"""

# ╔═╡ 355ebed2-42e6-41bb-b1db-a72d1aaae56f
mef_times = 1:24

# ╔═╡ 7bc3de9b-45e8-4a5b-a6b9-d816ee695bd9
refresh = true

# ╔═╡ 4a8c7752-6de8-4ea7-aafe-702b17507185
storage_pen = s_rel

# ╔═╡ 6e6b15b1-7685-4a20-9d94-dd703caa2fe9
η_vals = [0.9, 0.95, 0.99, 1.] #charge-discharge efficiency values

# ╔═╡ 47e2e177-3701-471f-ae3c-38276ec69370
begin

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

		println("...computing for η=$η")
		# Construct dynamic network
		C = total_demands * (storage_pen + δ) #δ is to avoid inversion errors
		P = C * c_rate
		net_dyn = make_dynamic(net, T, P, C, dyn_gmax, η, η)

		# Construct and solve OPF problem
		opf_dyn = DynamicPowerManagementProblem(net_dyn, d_dyn)
		solve!(opf_dyn, OPT, verbose=false)
		println("...The objective value is:")
		@show opf_dyn.problem.optval / T

		# Compute MEFs
		mefs = compute_mefs(opf_dyn, net_dyn, d_dyn, emissions_rates)
		println("...MEFs computed")
		for ind_t in 1:length(mef_times)
			results_η[:, :, ind_t, ind_η] .= mefs[ind_t]
		end

		meta_η[ind_η] = (opf=opf_dyn, net=net_dyn)
	end
		

end

# ╔═╡ 457c1959-94fa-4267-8645-3ed1409cd0a0
total_mefs_η = sum(results_η, dims=2)[:, 1, :, :];

# ╔═╡ dfa2a03b-6925-4be0-aeac-076c4cf25969
interesting_nodes = 2 : 2 : 30 

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
	plt_dyn_mef_eta = plot(subplots_η[highlighted_node_], xticks=xticks_hr)
	plot!(size=(600, 200), legend=:outertopright)
	plot!(title="node $(interesting_nodes[highlighted_node_])", bottommargin=3Plots.mm)
	plot!(ylabel="Δco2 (lbs) / Δmwh", xlabel="hour")
	
	# savefig(plt_dynamic_mef, 
	# 	"../img/storage_penetration_dynamic_mef_$(highlighted_node).png")
	plt_dyn_mef_eta
end

# ╔═╡ fb43471f-6aed-4fd6-a9e0-b165f6c77003
md"""
we also have to look into what *emissions* are, not only marginal emissions
"""

# ╔═╡ 9aded04b-e55f-4ebd-97c4-90c3adf62547
begin
	total_emissions_η = []
	for s in 1:length(η_vals)
		E = zeros(T)
		for t in 1:T
			E[t] = evaluate(meta_η[s].opf.g[t])' * emissions_rates
		end
		push!(total_emissions_η, E)
	end
	
	plt_total_emissions_η = plot()
	for (ind_η, η) in enumerate(η_vals)
		plot!(total_emissions_η[ind_η], label="η=$η", legend=:topleft, xticks=xticks_hr)
	end
	xlabel!("Hour")
	ylabel!("Total emissions")
	# legend!(:topleft)
	tot_emissions_vs_η = plot(η_vals, [sum(total_emissions_η[i]) for i in 1:length(η_vals)], xlabel="η", ylabel="Total Emissions", xticks=xticks_hr)
	
	# plt_tot_emissions_vs_η = plot(plt_total_emissions_η, tot_emissions_vs_η, layout=(2, 1), size=(600, 300))
	
	plt_tot_emissions_vs_η = plt_total_emissions_η
	
	plt_tot_emissions_vs_η
end

# ╔═╡ 30ec492c-8d21-43f6-bb09-32810494f21e
md"""
## How does storage penetration affect MEFs?
"""

# ╔═╡ 98a0d7c5-b1a8-4ebe-bb73-7ca88b475592
storage_penetrations = [0.0, 0.05, 0.10]

# ╔═╡ bbbb358c-e645-4989-bed3-73d9217f7447
md"""
The below cell is where computation of MEFs for different storage penetrations/different η vals happens
"""

# ╔═╡ 6f08828b-4c4b-4f50-bd40-35805a37aae0
begin
	println("Recomputing results")
	options = (
		c_rate=c_rate,
		renewable_penetration=renewable_penetration,
		storage_penetrations=storage_penetrations,
		mef_times=mef_times,
		emissions_rates=emissions_rates,
		d_dyn=d_dyn,
		η_vals=η_vals
	)

	meta = Dict()

	results = Dict()

	for (ind_η, η) in enumerate(η_vals)
		meta[ind_η]	= Dict()
		results[ind_η] = zeros(
			n, T, length(mef_times), length(storage_penetrations)
		)
		for (ind_s, s_rel) in enumerate(storage_penetrations)

			# Construct dynamic network
			C = total_demands * (s_rel + δ)
			P = C * c_rate
			net_dyn = make_dynamic(net, T, P, C, dyn_gmax, η, η)

			# Construct and solve OPF problem
			opf_dyn = DynamicPowerManagementProblem(net_dyn, d_dyn)
			solve!(opf_dyn, OPT, verbose=false)
			
			@show opf_dyn.problem.optval
			@show norm([evaluate(opf_dyn.g[t]) for t in 1:T])

			# Compute MEFs
			mefs = compute_mefs(opf_dyn, net_dyn, d_dyn, emissions_rates)
			for ind_t in 1:length(mef_times)
				results[ind_η][:, :, ind_t, ind_s] .= mefs[ind_t]
			end

			meta[ind_η][ind_s] = (opf=opf_dyn, net=net_dyn)
		end
	end
		
end

# ╔═╡ cd5fe410-6264-4381-b19f-31d050bc3930
begin
	
	plt_total_emissions = plot()
	for i in 1:length(η_vals)
		total_emissions = []
	for s in 1:length(storage_penetrations)
		E = 0
		for t in 1:T
			E += evaluate(meta[i][s].opf.g[t])' * emissions_rates
		end
		push!(total_emissions, E)
	end
		crt_η = η_vals[i]
		plot!(storage_penetrations, total_emissions, lw=4, label="η=$crt_η")
	end
	
	
	
	plot!(size=(650, 150), xlabel="storage capacity (% total demand)")
	plot!(ylabel="co2 emissions", xlim=(0, Inf))
	plot!(bottom_margin=5Plots.mm, left_margin=5Plots.mm)
	
	savefig(plt_total_emissions, "../img/storage_penetration_total_emissions.png")
	plt_total_emissions
end

# ╔═╡ 0740dc70-a532-4818-b09d-b3b8d60fa6ba
total_mefs = [sum(results[i], dims=2)[:, 1, :, :] for i in 1:length(η_vals)];

# ╔═╡ f26187fb-d4b2-4f0d-8a80-5d831e0de6c3
md"""
The cursor below allows choosing for a specific value of charge/discharge efficiency
"""

# ╔═╡ 19f4e0cc-0c93-42dc-8ee4-17f52d4e5e90

@bind idx_η Slider(1:1:length(η_vals))


# ╔═╡ 4d9a4b36-6b3d-4836-8501-7f46cd7ab5cc
md"""
The chosen value is:
"""

# ╔═╡ 75d956fc-bcf6-40fe-acd5-b5eef0fc7902
crt_η_ = η_vals[idx_η]

# ╔═╡ c6ee857d-8019-4c4f-bb07-a370a88ea3cf
md"""
MEFs as a function of time, for different charge/discharge efficiencies
"""

# ╔═╡ 6186798f-6711-4222-94bb-f53b2d0fad7d
begin
	subplots = []
	curves = 1:length(storage_penetrations)
	
	for (ind_plt, i) in enumerate(interesting_nodes)
		plt = plot(xticks=[6, 12, 18, 24], xlim=(1, 24))
		plot!(legend=nothing)
		plot!(mef_times, total_mefs[idx_η][i, :, curves], lw=4, alpha=0.8,
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
	plt_dynamic_mef = plot(subplots[highlighted_node], xticks=xticks_hr)
	plot!(size=(600, 200), legend=:outertopright)
	plot!(title="node $(interesting_nodes[highlighted_node]), η=$crt_η_", bottommargin=3Plots.mm)
	plot!(ylabel="Δco2 (lbs) / Δmwh", xlabel="hour")
	
	savefig(plt_dynamic_mef, 
		"../img/storage_penetration_dynamic_mef_$(highlighted_node).png")
	plt_dynamic_mef
end

# ╔═╡ 6fc320b1-b60d-4f49-89ab-bf029ead6b55
md"""
MEF as a function of consumption and emission times
"""

# ╔═╡ f7e0d09c-40bf-4936-987a-a3bcadae5487
begin
	node = interesting_nodes[highlighted_node] # we focus on a single node
	
	heatmap_subplts = []
	for s_idx in 1:length(storage_penetrations)
		
		crt_results = results[idx_η][node, :, :, s_idx]'
		# @show maximum(crt_results)
		subplt = heatmap(log10.(max.(crt_results, δ))', 
			c=:Blues_9, clim=(0, 4), colorbar=false,
			xlabel="Consumption Time",
			title="$(100*storage_penetrations[s_idx])% storage, η=$crt_η_", 					xticks=xticks_hr, yticks=xticks=xticks_hr
		)
		
		s_idx == 1 && plot!(ylabel="Emissions Time")
		s_idx == 3 && plot!(colorbar=true)
		
		push!(heatmap_subplts, subplt)
	end
	
	plt_emissions_heatmap = plot(heatmap_subplts..., 
		layout=Plots.grid(1, 3, widths=[.29, 0.29, 0.42]), 
		size=(650, 200), 
		bottom_margin=8Plots.pt
	)
	savefig(plt_emissions_heatmap, 
		"../img/storage_penetration_emissions_heatmap.png")
	plt_emissions_heatmap
	
	
	
	
end

# ╔═╡ edabacdd-8d25-4d64-9d4a-ecf1263ac02e
md"""
## Trying to understand where the cross patterns come from, and what they mean, etc.
"""

# ╔═╡ ee233685-474b-4162-bfef-5f3bdbe03c73
md"""
### Looking at charge/discharge patterns
"""

# ╔═╡ b5ac767e-ccde-4ee0-b50c-3584c5320e74
@show η_vals

# ╔═╡ a6df5496-a98a-49a8-9115-d73b86c9a7bd
@show storage_penetrations

# ╔═╡ ac6479c1-0cd9-4164-98e2-7081033f83e6
@bind s_idx_ Slider(1:1:length(storage_penetrations))

# ╔═╡ b5d3320e-0033-46b5-9aec-02b120a96cdc
md"""
The corresponding value of storage penetrations is
"""

# ╔═╡ f9ef362a-91c5-46ee-b4cb-1567b88f4dd7
storage_penetrations[s_idx_]

# ╔═╡ 4fb2a0e8-db61-4602-a01c-cfa7aeac6afc
meta_crt = meta[idx_η];

# ╔═╡ 62a605bc-26f8-43ac-8b66-ea9bf75b75c3
s_crt = meta_crt[s_idx_].opf.s;

# ╔═╡ 9c9d0fbc-58b8-440b-87ba-8cbf4dfaf71c
begin
	socs = zeros(n, T+1);
	for t in 1:T
		socs[:, t+1] = evaluate(s_crt[t])./C;
	end
end

# ╔═╡ 68f34570-6774-4232-b24c-bb634ec7d3a3
socs

# ╔═╡ 39be43d5-c14a-4340-8941-762cf873bb1b
md"""
In the above we see that some nodes have zero demand, and therefore zero capacity. The reason lies with `d_peak`, that has some components that are zero (i.e. some nodes don't consume any energy). 
"""

# ╔═╡ 88bf19ff-1635-4453-a611-a45775724c70
begin
	plot(t_axis, socs', ylim=(0, 1))
	title!("Relative SOC \n η = $crt_η_")
	xlabel!("Hour")
	ylabel!("SOC")
end

# ╔═╡ 0b9523aa-665f-4cb3-b985-b0d122a02a12
md"""
### Actually checking sensitivity at some points
"""

# ╔═╡ 3817ef32-586d-48fb-8f9b-ab4a9be7dbef
md"""
In order to check the sensitivity, we want to compute the emissions of the network at each time step for a given value of the demand and then re-compute the emissions for a slightly different value of the emissions at a given node and a given timestep. 

We can then visualize the dynamics of the system (e.g. charging patterns) to see if anything interesting is happening.
"""

# ╔═╡ 732cf0bc-24e7-44cd-a208-524a83706fdf
md"""
Some inspiration can be drawn from `sensitivity_demand_check()` function
"""

# ╔═╡ 3c5edbc5-8fc7-4d09-98a9-85f2efb699a8
node

# ╔═╡ e4c203b2-70b4-46f3-b025-3be46efa5e19
begin
cons_time = 22
em_times = [21, 22, 23];
end

# ╔═╡ 5365b74f-595f-4ade-a7af-e8dba53b84f7
md"""
Reference (as in computed) values
"""

# ╔═╡ a9b770e0-b137-40f7-b59a-35ad355b98bd
ref_mefs = results[idx_η][node, :, :, s_idx_]';

# ╔═╡ 8a5aa2a7-3d00-42c2-a26d-117380a46a7c
md"""
We see that there is a patch forming, i.e. the sensitivity is non zero even for nondiagonal elements.

Therefore, if we change the demand at `cons_time`, the emissions at those `em_times` should vary proportionately to what is below:
"""

# ╔═╡ b4c36ce0-7f31-453d-b3d5-ee54fcab7398
mefs_crt = ref_mefs[em_times, cons_time]

# ╔═╡ 4aed3df5-441b-445b-9277-a38690eb8603
npoints = 10

# ╔═╡ 0113ed6f-a39f-4167-a3c3-429bedbed7b0
ε = 1e-2

# ╔═╡ 0eac1f8d-2671-431c-84cc-a97225bb3120
evaluate(opf_dyn.g[1]).*emissions_rates

# ╔═╡ 87aa4097-62c5-439a-b5af-88a4cb5eb697
evaluate(opf_dyn.g[1])

# ╔═╡ 91f7d63c-9e30-4fd4-ab39-9fbf58d101dc
begin
	E_sensitivity = zeros(2npoints+1, length(emissions_rates), T);
	s_sensitivity = zeros(2npoints+1, n, T)
	g_sensitivity = zeros(2npoints+1, length(emissions_rates), T);
	crt_net = meta[idx_η][s_idx_].net;
	for i in -npoints:npoints
		d_crt = copy(d_dyn)
		d_crt[cons_time][node] = d_dyn[cons_time][node] * (1+i*ε)
		opf_ = DynamicPowerManagementProblem(crt_net, d_dyn)
		solve!(opf_, OPT, verbose=false)

		for t in 1:T
			s_sensitivity[i+npoints+1, :, t] = evaluate(opf_.s[t])
			g_sensitivity[i+npoints+1, :, t] = evaluate(opf_.g[t])
			E_sensitivity[i+npoints+1, :, t] = evaluate(opf_.g[t]).*emissions_rates
		end
	end
end

# ╔═╡ 77943ac8-36fe-4a13-a36d-db957780d869
begin
	E_ref = zeros(T)
	
	for t in 1:T
		E_ref[t] = evaluate(meta[idx_η][s_idx_].opf.g[t])' * emissions_rates
	end

end

# ╔═╡ 6fcd6e19-58c3-462d-964f-8cd3127b47a4
sum(E_sensitivity[npoints+1, :, em_times], dims=1)

# ╔═╡ 5c3117b3-7555-4366-9b43-bdba592a0877
md"""
it still does not add up... I don't understand what is going on...
Basically, the result in `E_sensitivity` at Δ = 0 should be the exact same as the result in `E_ref`
"""

# ╔═╡ af7f4492-4e97-4879-b151-2d8e38bf04d9
md"""
I have the feeling that the code does not result always in the same solution. One test I did was to run computation cells several times in a row and realized that I did not get the same emissions. 
"""

# ╔═╡ 2973af52-0bd0-4ba8-855d-297427627e22
E_ref[em_times]

# ╔═╡ 0c066c4c-2e9f-4bd6-b1de-5c61becf9ddd
md"""
I have currently ignored the mismatch between `E_sensitivity` and `E_ref` because I thought it's not that big of a deal. 
"""

# ╔═╡ 921cb430-ff44-4206-b5f0-60e2ef2a573a
begin
plot(
	[1-i*ε for i in -npoints:npoints], 
	E_sensitivity./repeat(E_sensitivity[npoints+1, :]', 2npoints+1, 1), 
	label=em_times', linestyle=[:solid :solid :dash], linewidth=3
)
	xlabel!("Relative demand change")
	ylabel!("Relative total emissions")
	title!("Consumption time=$cons_time")
end

# ╔═╡ b85b85d0-e1bc-4fc9-81cf-3792b55e3684
@bind t_display Slider(1:1:T)

# ╔═╡ 506e9360-2c25-4ea7-830b-68b4a6bf9026
t_display

# ╔═╡ d5db1b45-0c94-4c15-81ce-abb4ede4ad87
begin
	plot(
		[1-i*ε for i in -npoints:npoints], 
		[s_sensitivity[:, k, t_display]/(s_sensitivity[npoints+1, k, t_display]+1e-8) for k in 1:n], ylim=(.95, 1.05)
	)
	title!("Storage at time $t_display")
	xlabel!("Change in demand at node $node at time $cons_time")
	ylabel!("Change in storage at all nodes at time $t_display")

end

# ╔═╡ db1ae6bb-3923-4f7f-896e-0e42985ef380
begin
	plot(
		[1-i*ε for i in -npoints:npoints], 
		[E_sensitivity[:, k, t_display]./(E_sensitivity[npoints+1, k, t_display]+1e-6) for k in 1:length(emissions_rates)], ylim=(.95, 1.05)
		)
	title!("Emissions at time $t_display")
	xlabel!("Change in demand at node $node at time $cons_time")
	ylabel!("Change in emissions at all generators at time $t_display")
end

# ╔═╡ 2fbb3f97-42ae-4f58-887e-4b9a76c9cb49
begin
	plot(
		[1-i*ε for i in -npoints:npoints], 
		[g_sensitivity[:, k, t_display]./(g_sensitivity[npoints+1, k, t_display]+1e-8) for k in 1:length(emissions_rates)], ylim=(.99, 1.01)
		)
	title!("Generators at time $t_display")
	xlabel!("Change in demand at node $node at time $cons_time")
	ylabel!("Change in generation at all generators at time $t_display")
end

# ╔═╡ 9374abf1-78e4-4e60-875f-115cae7e7144
@show mefs_crt

# ╔═╡ 366d6ff5-4759-4f5d-8bd3-9a93a8f4eb9e
E_sensitivity

# ╔═╡ b5d3a530-2296-4e4a-83f7-98c4b6f116a6
md"""
showing the different charging patterns -- as a function of demand that changes
essentially, we need all the charging at the different amounts of time *around* the consumption time I guess. 

Basically I have to show: 
- charge at time t, node n (which node), as a function of the demand?
- Imagine you do all nodes, you have charge at time t (for multiple t's), node n (for all nodes), as a function of the demand

"""

# ╔═╡ 4d1aae3f-7976-41da-aa14-fff4c5574741
md"""
can you do the same with emissions?
"""

# ╔═╡ 7bdf5b4e-ea3b-4f79-9e16-ef6ea63d354a
md"""
### A plot "explaining" what we see?
"""

# ╔═╡ d8d1fb74-0018-4685-a283-e768ae877fe4
md"""
## Complete figure
"""

# ╔═╡ 5f73f4e6-4eff-41b9-b68d-3baa5e77e924
begin
l_ = @layout [
		a [b; c]
		d{.3h}
		]
Fig = plot(
		plt_time_series, plt_dynamic_mef, plt_tot_emissions_vs_η, plt_emissions_heatmap,  
		layout = l_, size=(800, 600), lw=2, legend=:outertopright, title = ["($i)" for j in 1:1, i in 1:7], titleloc = :right
	)
	
#save
savefig(Fig, 
	"../img/Fig_storage.png")
Fig
end

# ╔═╡ Cell order:
# ╟─c39005df-61e0-4c08-8321-49cc5fe71ef3
# ╠═5867a4eb-470a-4a8a-84e1-6f150de1dcde
# ╠═44275f74-7e7c-48d5-80a0-0f24609ef327
# ╠═db59921e-e998-11eb-0307-e396d43191b5
# ╠═0aac9a3f-a477-4095-9be1-f4babe1e2803
# ╠═257a6f74-d3c3-42eb-8076-80d26cf164ca
# ╠═113e61a9-3b21-48d0-9854-a2fcce904e8a
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
# ╠═9e5a1672-d452-42f5-ba5a-a2fa0b1eaada
# ╠═522e95d5-15a8-47ec-a79c-c4cc17cf86fd
# ╠═bfa4a8f9-cfcd-4e22-b2dc-751226f3a73c
# ╟─a8ccbc8e-24e6-4214-a179-4edf3cf26dad
# ╠═e8ee5cbb-4afc-4737-b006-90071f6138cd
# ╟─c8f71644-9371-443b-b976-1734cc7ae583
# ╠═57a37eb0-a547-428c-9f8c-e5f3e30f5260
# ╠═491da4e6-03e6-49d3-907d-43ddcffdfb42
# ╟─147a8ae3-a910-4d76-9689-4e77f7012914
# ╠═9a45c280-d8a4-4849-8977-2dc763ed633b
# ╟─5bd157b8-3dac-464e-85b0-c63830a7e817
# ╠═e55e0566-7bd6-4126-8f60-0d940f6d8111
# ╠═9deb6bdf-54f4-42ee-855d-7e82cef6f4bb
# ╠═af6d4ef0-f59f-42be-af36-7cf447478e4c
# ╠═dd9efed9-f6ed-43fb-93f2-4d21d1360091
# ╟─3f9eb091-059c-44a5-9b50-ae3cabe24060
# ╠═355ebed2-42e6-41bb-b1db-a72d1aaae56f
# ╠═7bc3de9b-45e8-4a5b-a6b9-d816ee695bd9
# ╟─4a8c7752-6de8-4ea7-aafe-702b17507185
# ╠═6e6b15b1-7685-4a20-9d94-dd703caa2fe9
# ╠═47e2e177-3701-471f-ae3c-38276ec69370
# ╠═457c1959-94fa-4267-8645-3ed1409cd0a0
# ╟─dfa2a03b-6925-4be0-aeac-076c4cf25969
# ╟─a7e75e49-5031-4cc4-b96e-6227277ec3ba
# ╠═d9617524-76c3-447d-9b94-0a690f83a7b9
# ╟─fb43471f-6aed-4fd6-a9e0-b165f6c77003
# ╠═9aded04b-e55f-4ebd-97c4-90c3adf62547
# ╟─30ec492c-8d21-43f6-bb09-32810494f21e
# ╠═98a0d7c5-b1a8-4ebe-bb73-7ca88b475592
# ╟─bbbb358c-e645-4989-bed3-73d9217f7447
# ╠═6f08828b-4c4b-4f50-bd40-35805a37aae0
# ╠═cd5fe410-6264-4381-b19f-31d050bc3930
# ╠═0740dc70-a532-4818-b09d-b3b8d60fa6ba
# ╟─f26187fb-d4b2-4f0d-8a80-5d831e0de6c3
# ╟─19f4e0cc-0c93-42dc-8ee4-17f52d4e5e90
# ╟─4d9a4b36-6b3d-4836-8501-7f46cd7ab5cc
# ╟─75d956fc-bcf6-40fe-acd5-b5eef0fc7902
# ╟─c6ee857d-8019-4c4f-bb07-a370a88ea3cf
# ╟─6186798f-6711-4222-94bb-f53b2d0fad7d
# ╟─d27ef0d8-70b2-4897-9000-8fa70b1862fc
# ╟─6fc320b1-b60d-4f49-89ab-bf029ead6b55
# ╠═f7e0d09c-40bf-4936-987a-a3bcadae5487
# ╟─edabacdd-8d25-4d64-9d4a-ecf1263ac02e
# ╟─ee233685-474b-4162-bfef-5f3bdbe03c73
# ╠═b5ac767e-ccde-4ee0-b50c-3584c5320e74
# ╠═a6df5496-a98a-49a8-9115-d73b86c9a7bd
# ╠═ac6479c1-0cd9-4164-98e2-7081033f83e6
# ╟─b5d3320e-0033-46b5-9aec-02b120a96cdc
# ╠═f9ef362a-91c5-46ee-b4cb-1567b88f4dd7
# ╠═4fb2a0e8-db61-4602-a01c-cfa7aeac6afc
# ╠═62a605bc-26f8-43ac-8b66-ea9bf75b75c3
# ╠═9c9d0fbc-58b8-440b-87ba-8cbf4dfaf71c
# ╠═68f34570-6774-4232-b24c-bb634ec7d3a3
# ╟─39be43d5-c14a-4340-8941-762cf873bb1b
# ╠═88bf19ff-1635-4453-a611-a45775724c70
# ╟─0b9523aa-665f-4cb3-b985-b0d122a02a12
# ╟─3817ef32-586d-48fb-8f9b-ab4a9be7dbef
# ╟─732cf0bc-24e7-44cd-a208-524a83706fdf
# ╠═3c5edbc5-8fc7-4d09-98a9-85f2efb699a8
# ╠═e4c203b2-70b4-46f3-b025-3be46efa5e19
# ╟─5365b74f-595f-4ade-a7af-e8dba53b84f7
# ╠═a9b770e0-b137-40f7-b59a-35ad355b98bd
# ╟─8a5aa2a7-3d00-42c2-a26d-117380a46a7c
# ╠═b4c36ce0-7f31-453d-b3d5-ee54fcab7398
# ╠═4aed3df5-441b-445b-9277-a38690eb8603
# ╠═0113ed6f-a39f-4167-a3c3-429bedbed7b0
# ╠═0eac1f8d-2671-431c-84cc-a97225bb3120
# ╠═87aa4097-62c5-439a-b5af-88a4cb5eb697
# ╠═91f7d63c-9e30-4fd4-ab39-9fbf58d101dc
# ╠═77943ac8-36fe-4a13-a36d-db957780d869
# ╠═6fcd6e19-58c3-462d-964f-8cd3127b47a4
# ╟─5c3117b3-7555-4366-9b43-bdba592a0877
# ╟─af7f4492-4e97-4879-b151-2d8e38bf04d9
# ╠═2973af52-0bd0-4ba8-855d-297427627e22
# ╟─0c066c4c-2e9f-4bd6-b1de-5c61becf9ddd
# ╠═921cb430-ff44-4206-b5f0-60e2ef2a573a
# ╟─b85b85d0-e1bc-4fc9-81cf-3792b55e3684
# ╟─506e9360-2c25-4ea7-830b-68b4a6bf9026
# ╠═d5db1b45-0c94-4c15-81ce-abb4ede4ad87
# ╠═db1ae6bb-3923-4f7f-896e-0e42985ef380
# ╠═2fbb3f97-42ae-4f58-887e-4b9a76c9cb49
# ╠═1a7af4e0-2608-422c-bafe-d200f30bc4f3
# ╠═9374abf1-78e4-4e60-875f-115cae7e7144
# ╠═366d6ff5-4759-4f5d-8bd3-9a93a8f4eb9e
# ╠═b5d3a530-2296-4e4a-83f7-98c4b6f116a6
# ╠═4d1aae3f-7976-41da-aa14-fff4c5574741
# ╟─7bdf5b4e-ea3b-4f79-9e16-ef6ea63d354a
# ╟─d8d1fb74-0018-4685-a283-e768ae877fe4
# ╟─5f73f4e6-4eff-41b9-b68d-3baa5e77e924
