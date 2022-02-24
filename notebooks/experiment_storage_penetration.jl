### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ db59921e-e998-11eb-0307-e396d43191b5
begin
	import Pkg
	Pkg.activate();
	using Random, Distributions
	using Convex, ECOS, Gurobi
	using Plots
	using PlutoUI
	using JLD
	using LinearAlgebra
	using LightGraphs, SimpleWeightedGraphs, GraphPlot
end;

# ╔═╡ 0aac9a3f-a477-4095-9be1-f4babe1e2803
begin
	using Revise
	using CarbonNetworks
end

# ╔═╡ 571d1cff-7311-4db8-8ac3-9e10afefaf18
using LaTeXStrings

# ╔═╡ 9af03cdd-8cb8-4221-b962-05e6ae1634cc
using HDF5

# ╔═╡ c39005df-61e0-4c08-8321-49cc5fe71ef3
md"""
## Description
"""

# ╔═╡ 0510ec8c-2f1b-4704-bb59-cd8a67ef0dc5
md"""
## *TODO*
"""

# ╔═╡ 44275f74-7e7c-48d5-80a0-0f24609ef327
md"""
## Loading
"""

# ╔═╡ a32d6a56-8da8-44b0-b659-21030692630a
begin
	ECOS_OPT = () -> ECOS.Optimizer(verbose=false)
	# GUROBI_ENV = Gurobi.Env()
	# GUROBI_OPT = Convex.MOI.OptimizerWithAttributes(
	# 	() -> Gurobi.Optimizer(GUROBI_ENV), "LogToConsole" => false
	# )
	# ECOS_OPT_2 = Convex.MOI.OptimizerWithAttributes(
	# 	ECOS.Optimizer, "maxit"=> 100, "reltol"=>1e-6, "LogToConsole"=>false
	# 	)
	OPT = ECOS_OPT
	δ = 1e-4
end;

# ╔═╡ 257a6f74-d3c3-42eb-8076-80d26cf164ca
theme(:default, label=nothing, 
		tickfont=(:Times, 8), guidefont=(:Times, 8), legendfont=(:Times, 8), titlefont=(:Times,8), framestyle=:box)

# ╔═╡ 113e61a9-3b21-48d0-9854-a2fcce904e8a
xticks_hr = [0, 6, 12, 18, 24]

# ╔═╡ 9bd515d4-c7aa-4a3d-a4fb-28686290a134
md"""
## Generate data
"""

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

# ╔═╡ 75dfaefd-abec-47e2-acc3-c0ff3a01048e
md"""
### Network
"""

# ╔═╡ 7d738cd3-1edf-4348-aa16-5ade0f87ae6f
#define which nodes have renewable generation
renew_nodes = [12, 14, 15, 18]

# ╔═╡ 63db8d8b-8229-47e0-a309-75813120aaab
placed_gen_nodes = [9, 11, 6, 8, 28, 7] 

# ╔═╡ 35dccacb-8780-4f37-8dc5-29064a84d49c
# allocating low/med capacity to some edges. Most likely, most of those allocations are completely wrong
begin
	add_low = [1, 6, 13, 21, 25, 26, 33, 30, 32, 9, 3]
	add_med = [35, 7, 12, 38, 36, 11]
	level_low = .4
	level_med = .5
end;

# ╔═╡ f999d732-14b3-4ac5-b803-3df7a96ef898
begin
	Random.seed!(2)# Random.seed!(2) - for n_low_capacity = 10 - .4
	net, d_peak, _ = load_synthetic_network("case30.m")
	n = length(d_peak)
	
	if placed_gen_nodes != nothing
		net.B = I(n)[:, placed_gen_nodes]
	end
	_, l_no_renew = size(net.B)

	gen_nodes = findall(vec(sum(net.B, dims = 2)) .> 0)
	
	# Limit line flows and generation capacities
	capacities = ones(size(net.A)[2]) * .6
	capacities[[29, 31]] .= .4 # minimum necessary, edges around gen 22
	capacities[add_low] .= level_low
	capacities[add_med] .= level_med

	net.pmax = net.pmax .* capacities 
	
	net.gmax *= 2.
	
	# Add generators for each renewable source
	net.B = [net.B I(n)[:, renew_nodes]]
	net.gmax = [net.gmax; zeros(length(renew_nodes))]
	net.fq = [net.fq; zeros(length(renew_nodes))]
	net.fl = [net.fl; zeros(length(renew_nodes))]
	
	"Network loaded."
end


# ╔═╡ 0d42b50d-993e-4eea-9025-2b7479bb3b0e
begin
	# color 1: non rewable nodes
	# color 2: renewable nodes
	# color 3 = generator nodes -- TODO
	pal = palette(:tab20)
	graph_colors = [pal[6], pal[8], pal[16]]
	# graph_colors = [:green, :blue, :red]
	node_cols = []
	for k in 1:n
		if k in renew_nodes
			push!(node_cols, graph_colors[1])
		elseif k in gen_nodes
			push!(node_cols, graph_colors[2])
		else
			push!(node_cols, graph_colors[3])
		end
	end
	# node_cols = [k in renew_nodes ? graph_colors[2] : graph_colors[1] for k in 1:n]

end

# ╔═╡ 23690382-3d30-46e3-b26a-a30875be78ec
begin
	#edges with low capacity
	Random.seed!(17)
	
	G = SimpleWeightedGraph(n)
	for j in 1:size(net.A, 2)
		inds = findall(x -> x != 0, net.A[:, j])
		if net.pmax[j] > δ 
			add_edge!(G, inds[1], inds[2], net.pmax[j])
		end
	end
	
	edge_colors = [
		k < .6 ? colorant"orange" : colorant"lightgray" for k in capacities
		]
	
	Gplot = gplot(G, nodelabel=1:n, nodefillc=node_cols, edgestrokec=edge_colors)
	Gplot
end

# ╔═╡ c92e091e-1111-4bda-8d76-543ad1dd8121
#nodes with demand
begin
	Random.seed!(17)
	node_demand_color = [
		d_peak[k] == 0 ? colorant"lightgray" : colorant"orange" for k in 1:n
	]
	Gplot_demand = gplot(G, nodelabel=1:n, nodefillc=node_demand_color)
	Gplot_demand
end

# ╔═╡ 379598be-2c42-4c29-8a24-91b74592da0f
# printing origins and destinations of all edges
begin
	i = 1
	for e in edges(G)
	println(i)
	println(e)
	i+=1
	end
end

# ╔═╡ 39cdac4b-87bb-441c-99e3-402b40bef7d3
_, m, l = get_problem_dims(net);

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
	
	plt2 = plot(renew_data[:, renew_labels.=="solar"], lw=2, palette=:Reds, title="Renewable Data", xlabel="Hour", ylabel="Energy [?]", xticks=xticks_hr)
	plot!(renew_data[:, renew_labels.=="wind"], lw=2, palette=:Greens, title="Renewable Data", xlabel="Hour", ylabel="Energy [?]", xticks=xticks_hr)
	
	plt_time_series = plot(plt1, plt2, layout=(2, 1), size=(600, 300))
	
	savefig(plt_time_series, "../img/storage_penetration_time_series.png")
	plt_time_series
end

# ╔═╡ c82ef027-740a-49b1-93d2-1554c411a896
renewable_penetration = 0.1

# ╔═╡ 9873abd2-fc2f-4a35-8689-51d558debbf0
md"""the number of nodes that have zero demand is $(sum(d_peak.==0))"""

# ╔═╡ 0239e1da-caf5-4593-af1b-5d1e8d2f2b3e
begin
	Random.seed!(1)
	demand_profiles = rand(1:n_demand, n)
	
	d_dyn = [
		1.0 * d_peak .* demand_data[t, demand_profiles]
		for t in 1:T
	]
	
	total_demands = sum(d_dyn)
end;

# ╔═╡ 9e5a1672-d452-42f5-ba5a-a2fa0b1eaada
n_renew = size(renew_data, 2)

# ╔═╡ 522e95d5-15a8-47ec-a79c-c4cc17cf86fd
begin
	Random.seed!(3)
	# here we allocate renewable generation to every node. 
	# we probably want to sparsify it
	# renew_profiles = vcat(rand(1:n_renew, length(renew_nodes)))
	renew_profiles = [1, 2, 3, 4]

	@assert length(renew_profiles) == length(renew_nodes)

	get_renew = (t, renew_profile) -> 
		renew_profile == 0 ? 0.0 : renew_data[t, renew_profile]
	
	renew_gmax = [
		get_renew.(t, renew_profiles) .* total_demands[renew_nodes] * renewable_penetration for t in 1:T
	]
	
	dyn_gmax = [deepcopy(net.gmax) for t in 1:T]
	for t in 1:T
		dyn_gmax[t][end-length(renew_nodes)+1:end] .= max.(δ, renew_gmax[t])
	end
end;

# ╔═╡ bfa4a8f9-cfcd-4e22-b2dc-751226f3a73c
begin
	net_d_dyn = hcat(d_dyn...)' - ((I(n)[:, renew_nodes]) * hcat(renew_gmax...))'
	Random.seed!(7)
	plot(size=(600, 200))
	plot!(net_d_dyn[:, rand(1:n, 10)], lw=2, c=repeat(renew_profiles, outer=T)')
	xlabel!("Time")
	ylabel!("Net demand")
end

# ╔═╡ a8ccbc8e-24e6-4214-a179-4edf3cf26dad
md"""
### Carbon emissions data
"""

# ╔═╡ 496135ec-f720-4d43-8239-d75cc7616f58
md"""Emissions rates:"""

# ╔═╡ 806819d5-7b40-4ca4-aa1e-f1cf0a9a7f3f
begin
# From deChalendar - Tracking emissions codebase
# UNK is 2017 average US power grid intensity according to Schivley 2018
# unit is kg / MWh
	EMISSIONS_FACTORS = Dict(
		"WAT" => 4,
        "NUC" => 16,
        "SUN" => 46,
        "NG" => 469,
        "WND" => 12,
        "COL" => 1000,
        "OIL" => 840,
        "OTH" => 439,
        "UNK" => 439,
        "BIO" => 230,
        "GEO" => 42,
	)
	
	gen_type = ["NG", "NUC", "COL", "OIL", "NG", "COL"]
	renew_tags = []
	rl = renew_labels[renew_profiles]
	for l in rl
		if l == "solar"
			push!(renew_tags, "SUN")
		elseif l == "wind"
			push!(renew_tags, "WND")
		else 
			println("WRONG TAG")
		end
	end

	tags = vcat(gen_type, renew_tags)

	emissions_rates = [EMISSIONS_FACTORS[tag] for tag in tags]
	emissions_rates
end

# ╔═╡ c8f71644-9371-443b-b976-1734cc7ae583
md"""
## What is the impact of charging/discharging efficiency on storage dynamics?
"""

# ╔═╡ 57a37eb0-a547-428c-9f8c-e5f3e30f5260
# Parameters
begin
	c_rate = .25 #charging rate
	s_rel = 0.0 #storage penetration
	mef_times = 1:24
	η_c = 1
	η_d = η_c
	node_storage = 0
end;

# ╔═╡ c5f31bf6-5cca-4b0e-86f0-b44321fec874
# Simply solving the problem
begin
	println("Solving first problem")
	mefs = zeros(
			n, T, length(mef_times)
		)
	
	# Construct dynamic network
	C = sum(d_dyn) .* (s_rel + δ)/T .+ δ
	
	if node_storage >= 1 # to update
		nodes_on = zeros(n)
		nodes_on[node_storage] = 1
		C = C.* nodes_on
	end
	P = C*c_rate


	net_dyn = make_dynamic(net, T, P, C, η_c);
|
	# Construct and solve OPF problem
	opf_dyn = DynamicPowerManagementProblem(net_dyn, d_dyn)
	solve!(opf_dyn, OPT, verbose=false)

	# if opf_dyn.problem.status != Convex.MOI.OPTIMAL
		@show opf_dyn.problem.status
	# end

	# Compute MEFs
	mefs_ = compute_mefs(opf_dyn, net_dyn, d_dyn, emissions_rates)

	for ind_t in 1:T
		mefs[:, :, ind_t] .= mefs_[ind_t];
	end
	status_ = opf_dyn.problem.status
	status_
end

# ╔═╡ 15245a9a-8a42-4ec1-ba2c-9bf95d6c85fa
status_

# ╔═╡ af6d4ef0-f59f-42be-af36-7cf447478e4c
begin
	s_vals = zeros(n, T+1)
	g_vals = zeros(l, T)
	p_vals = zeros(m, T)
	d_vals = zeros(n, T)
	for t in 1:T
    	s_vals[:, t+1] = evaluate(opf_dyn.s[t])./C
		g_vals[:, t] = evaluate(opf_dyn.g[t])./dyn_gmax[t]
		p_vals[:, t] = evaluate(opf_dyn.p[t])./net.pmax
		d_vals[:, t] = d_dyn[t]
	end
end

# ╔═╡ 06498d79-f20a-4691-8c72-8f1f962e6a6f
#node to be displayed here below
node_single = 21

# ╔═╡ 49aea0fc-1e0b-4e9c-8d7b-5d3727e95e1f
begin
	subplts = []
		
	
	
	plot_p = plot(abs.(p_vals)', xlabel="t", ylabel="p(t)/pmax", ylim=(-0.1, 1.1));
	push!(subplts, plot_p)
	
	t_axis = [t for t in 0:T]
	s_subplt = plot(t_axis, s_vals', ylim=(-.1, 1.1))
	title!("Relative SOC \n η = $η_c")
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
	
	
	crt_mefs = mefs[node_single, :, :]
	lim = max(abs(maximum(crt_mefs)), abs(minimum(crt_mefs)));
	clims = (-lim, lim)
	subplt = heatmap(crt_mefs, 
		c=:balance, colorbar=true,
		xlabel="Consumption Time",
		title="$(100*s_rel)% storage, η=$η_c, node=$node_single", 
		xticks=xticks_hr, yticks=xticks_hr, 
		clim=clims
	)
		
	plot!(ylabel="Emissions Time")
	# plot!(colorbar=true)
	push!(subplts, subplt)
	
	
	lay = @layout [Plots.grid(1, 3); Plots.grid(1, 2)]
	
	plots = plot(subplts..., 
		layout=lay
		# size=(650, 200), 
		# bottom_margin=8Plots.pt
	)

	
	plots
end

# ╔═╡ 871ca595-9d9e-4f06-9119-c2b07bfdeb04
begin
	tt = 17
	mef_bar_plt = bar(
		[sum(mefs[k, tt, :]) for k in 1:n]
			)
	xlabel!("Node number")
	ylabel!("MEF [kg/MWh]")
	mef_bar_plt
end

# ╔═╡ 4d056537-3ec5-47d1-9f51-a603b8ae91ff
# data to save
MEFS_to_save = [sum(mefs[k, tt, :]) for k in 1:n]

# ╔═╡ 0d99fc04-0353-4170-a23e-f21460ceaf7e
interesting_nodes = [10, 17, 19, 21, 22, 23, 24]

# ╔═╡ 38e73213-d399-43e5-80e8-851b6cf3299d
begin
	ps = evaluate(opf_dyn.p[tt])./net.pmax
	idx = Int.(round.(abs.(100*ps))) .+1
	cmap = colormap("Blues", 101) 
	edge_flow_colors = [cmap[id] for id in idx]
end;

# ╔═╡ f15393ae-5c04-473e-ae28-ead5554896a3
begin
	Random.seed!(17)
	g_edge_flow_plot = gplot(
		G, nodelabel=1:n, nodefillc=node_cols, edgestrokec=edge_flow_colors
	)
	g_edge_flow_plot
end

# ╔═╡ 83e265f4-0a62-4c53-9203-6a7768c6220e
g_edge_flow_plot

# ╔═╡ 1bd70b8a-90ff-45f3-8ef4-dee0340a02f0
md"""
## What is the impact of congestion on MEFs
"""

# ╔═╡ a1517984-dfdb-4be7-94cd-16fbf2dbc453
begin

	mef_congestion = Dict()
	congestion = Dict()
	factors = [1 1.5 2 2.5 3 5 10]
	for f in factors
		println("Solving problem with line capacities at $f x")
		net_crt = deepcopy(net_dyn)
		net_crt.pmax *= f
		# Construct and solve OPF problem
		opf_crt = DynamicPowerManagementProblem(net_crt, d_dyn)
		solve!(opf_crt, OPT, verbose=false)

		@show opf_crt.problem.status

		mefs_crt = compute_mefs(opf_crt, net_crt, d_dyn, emissions_rates)
		mef_congestion[f] = mefs_crt
		congestion[f] = [abs.(evaluate(opf_crt.p[t])./net_crt.pmax[t]) for t in 1:T]
	end


end

# ╔═╡ 58368b3b-c042-4f37-81f0-f8dd0388e5aa
# prepare plotting data
begin
	# x = congestion level, over one time step (mean, min, max)
	x_mean, x_min, x_max = [], [], []
	# y = different statistics of mefs (median, mean, min, max)
	y_mean, y_median, y_min, y_max = [], [], [], []
	# num congested lines
	n_cong = []
	
	for f in factors
		for t in 1:T
			mefs_cong = mef_congestion[f][t]
			cong_level = congestion[f][t]
			push!(x_mean, mean(cong_level))
			push!(x_min, minimum(cong_level))
			push!(x_max, maximum(cong_level))
			push!(y_mean, mean(mefs_cong))
			push!(y_median, median(mefs_cong))
			push!(y_min, minimum(mefs_cong))
			push!(y_max, maximum(mefs_cong))
			push!(n_cong, sum(cong_level .>= .999))
		end
	end
end

# ╔═╡ 0ffb55c7-1371-45ac-9823-406eb86bd8e4
begin
	scatter(x_mean, y_mean, label = "mean")
	scatter!(x_mean, y_min, label = "min")
	scatter!(x_mean, y_max, label = "max")
	plot!(ylims = (-1000, 7000))
	plot!(x_mean, maximum(emissions_rates) .* ones(length(x_mean)))
	xlabel!("Average congestion level across the network")
	ylabel!("MEF")
end

# ╔═╡ 794b4672-0d30-4df6-9fa1-2d68d34c2502
begin
	scatter(x_max, y_mean, label = "mean")
	scatter!(x_max, y_min, label = "min")
	scatter!(x_max, y_max, label = "max")
	plot!(x_max, maximum(emissions_rates) .* ones(length(x_mean)))
	plot!(ylims = (-1000, 2000))
	xlabel!("Max congestion level across the network")
	ylabel!("MEF")
end

# ╔═╡ 2bde6957-46e5-4697-9451-40b135f8523c
begin
	# scatter(n_cong, y_mean)
	scatter(n_cong, y_max, label = "max", legend=:topleft)
	scatter!(n_cong, y_min, label = "min")
	xlabel!("Number of congested lines")
	ylabel!("MEFs")
end

# ╔═╡ a4e2b01f-6ea9-4129-954b-3404b0e16e8f
begin
	function get_ribbon_params(y, n_cong)
		ymid = []
		dev_up, dev_down = [], []
		x = sort(unique(n_cong))
		for k in x
			idx = findall(n_cong .== k)
			ymid_ = median(y[idx])
			dev_up_ = maximum(y[idx] .- ymid_)
			dev_down_ = maximum(ymid_ .- y[idx])
			push!(ymid, ymid_)
			push!(dev_up, dev_up_)
			push!(dev_down, dev_down_)
		end

		return ymid, dev_up, dev_down, x
	end
end

# ╔═╡ cf7430e5-c46b-4db7-adf0-8058f9ac8640
#ribbon plot
begin
	ymid1, devup1, devdown1, x = get_ribbon_params(y_max, n_cong)
	ymid2, devup2, devdown2, _ = get_ribbon_params(y_min, n_cong)
	plot(x, ymid1, ribbon = (devdown1, devup1), lw = 3, label = "Max", legend=:topleft)
	plot!(x, ymid2, ribbon = (devdown2, devup2), lw=3, label = "Min")
	xlabel!("Number of congested lines")
	ylabel!("MEFs")
	# plot!(ylims = (-2000, 10000))
end

# ╔═╡ 30ec492c-8d21-43f6-bb09-32810494f21e
md"""
## How does storage penetration affect MEFs?
"""

# ╔═╡ 856a78d9-7b4c-453b-b73b-c81eee014e52
RUN_BIG_CELL2 = true

# ╔═╡ 98a0d7c5-b1a8-4ebe-bb73-7ca88b475592
storage_penetrations = [0.0, 0.05, 0.10]

# ╔═╡ 71cd6842-7df3-4b0a-ae23-c49404ddf523
# η_vals = [0.95, 0.99, 1.] 
η_vals = [1.]

# ╔═╡ bbbb358c-e645-4989-bed3-73d9217f7447
md"""
The below cell is where computation of MEFs for different storage penetrations/different η vals happens
"""

# ╔═╡ 6f08828b-4c4b-4f50-bd40-35805a37aae0
begin
	if RUN_BIG_CELL2
	println("--------------------------------------------------------")
	println("Recomputing results for different storage pens")
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
			println("s = $s_rel, η = $η")
			# Construct dynamic network
			C = sum(d_dyn) * (s_rel + δ) .+ δ
			P = C * c_rate
			net_dyn = make_dynamic(net, T, P, C, dyn_gmax, η, η)

			# Construct and solve OPF problem
			opf_dyn = DynamicPowerManagementProblem(net_dyn, d_dyn)
			solve!(opf_dyn, OPT, verbose=false)
			
			if opf_dyn.problem.status != Convex.MOI.OPTIMAL
				@show opf_dyn.problem.status
			end
				
			# Compute MEFs
			mefs = compute_mefs(opf_dyn, net_dyn, d_dyn, emissions_rates)
			for ind_t in 1:length(mef_times)
				results[ind_η][:, :, ind_t, ind_s] .= mefs[ind_t]
			end

			meta[ind_η][ind_s] = (opf=opf_dyn, net=net_dyn)
		end
	end
	end
		
end

# ╔═╡ 15293269-580f-4251-be36-4be6ba8c5a46
md"""Influence of charging efficiency η on the total emissions of the system"""

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
		plot!(storage_penetrations, total_emissions, lw=2, label="η=$crt_η", ls=:dash, markershape=:circle, ms=4)
	end
	
	
	
	plot!(size=(650, 300), xlabel="storage capacity (% total demand)", legend=:bottomright)
	plot!(ylabel="co2 emissions", xlim=(0, Inf))
	plot!(bottom_margin=5Plots.mm, left_margin=5Plots.mm)
	
	savefig(plt_total_emissions, "../img/storage_penetration_total_emissions.png")
	plt_total_emissions
end

# ╔═╡ 0740dc70-a532-4818-b09d-b3b8d60fa6ba
total_mefs = [sum(results[i], dims=2)[:, 1, :, :] for i in 1:length(η_vals)];

# ╔═╡ f26187fb-d4b2-4f0d-8a80-5d831e0de6c3
md"""
Define the charge/discharge efficiency to use in the future...?
"""

# ╔═╡ 4925c50b-12c0-4217-94de-bdcc72c01ccf
idx_η = 1

# ╔═╡ 75d956fc-bcf6-40fe-acd5-b5eef0fc7902
crt_η_ = η_vals[idx_η];

# ╔═╡ 4d9a4b36-6b3d-4836-8501-7f46cd7ab5cc
md"""
The current value for η is $crt_η_
"""

# ╔═╡ c6ee857d-8019-4c4f-bb07-a370a88ea3cf
md"""
MEFs as a function of time, for different charge/discharge efficiencies
"""

# ╔═╡ 11b97da8-0e0d-48dd-a442-4ffa655bec61
interesting_nodes

# ╔═╡ 6186798f-6711-4222-94bb-f53b2d0fad7d
begin
	subplots_mef_storage = Dict()
	curves_ = [1, 3]
	
	for (ind_plt, i) in enumerate(interesting_nodes)
		subplots_mef_storage[i] = plot(
			mef_times, total_mefs[idx_η][i, :, curves_], 
			lw=4, alpha=0.8,
			xlim=(1, 24),
			ylim=(-1000, 4000), 
			yticks = [0, 2000, 4000], 
			xticks=xticks_hr,  
			xlabel=L"t_c", 
			ylabel="MEF", 
			legend=:topright
		)
		# plot!(legend=nothing)
		# plot!(
			# 
			# labels=storage_penetrations[curves]', 
		# )
		# title!("Node = $i")
		
		# ind_plt in [1, 4] && plot!(ylabel="Δco2 (lbs) / Δmwh")
		# ind_plt in 4:6 && plot!(xlabel="hour")
		# ind_plt in [3] && plot!(legend=:topright)
		# xlabel!(L"t_c")
		# ylabel!("MEF")
		# plot!(legend=:topright)
		
		# subplots_mef_storage[i] = plt
	end

	# plot(subplots_mef_storage..., layout=(2, Int(length(interesting_nodes)/2)), leftmargin=4Plots.mm, size=(800, 400))
end

# ╔═╡ d27ef0d8-70b2-4897-9000-8fa70b1862fc
# begin
	
# 	highlighted_node = 21
# 	nd = findall(interesting_nodes .== highlighted_node)[1]
# 	plt_dynamic_mef = plot(subplots_mef_storage[nd], xticks=xticks_hr)
# 	plot!(size=(600, 200), legend=:outertopright)
# 	plot!(title="node $(interesting_nodes[nd]), η=$crt_η_", bottommargin=3Plots.mm)
# 	plot!(ylabel="Δco2 (kg) / Δmwh", xlabel="hour")
	
# 	savefig(plt_dynamic_mef, 
# 		"../img/storage_penetration_dynamic_mef_$(highlighted_node).png")
# 	plt_dynamic_mef
# end

# ╔═╡ 420919bc-f217-4357-bdbb-83e25e83ba56
subplots_mef_storage[21]

# ╔═╡ 5f77f4a9-ff5a-4515-9016-bf36571225c7
subplots_mef_storage[23]

# ╔═╡ 1da34733-fee3-42e1-b5e0-cac3f5f196c9
nodes_heatmaps = [21, 23]

# ╔═╡ 57643580-0a1d-4ad5-ba21-533dcbd73c2f
# data to save mef_storage
begin
data_total_mef = ([k for k in mef_times], [total_mefs[idx_η][n, :, :] for n in nodes_heatmaps])
end

# ╔═╡ f7e0d09c-40bf-4936-987a-a3bcadae5487
begin
	plt_emissions_heatmap = Dict()
	heatmap_subplts = Dict()

	clims_hm = Dict()
	clims_hm[21] = (-500, 500)
	clims_hm[23] = (-3000, 3000)
	for node in nodes_heatmaps
	# node = node_matrix# we focus on a single node
	
	plots = []
	lims = []
	heatmap_subplts[node] = Dict()

	# for idx_η in 1:length(η_vals)
	# 	crt_η_ = η_vals[idx_η]
	# 	for s_idx in 1:length(storage_penetrations)
	# 		crt_results = results[idx_η][node, :, :, s_idx]'
	# 		lim_crt = max(abs(minimum(crt_results)), abs(maximum(crt_results)));
	# 		push!(lims, lim_crt)
	# 	end
	# end

	# clims = (-maximum(lims), maximum(lims))
	# clims = (-800, 800)
		
	for idx_η in 1:length(η_vals)
		crt_η_ = η_vals[idx_η]
		heatmap_subplts[node][idx_η] = Dict()
		for s_idx in 1:length(storage_penetrations)
			crt_results = results[idx_η][node, :, :, s_idx]'
			
			# lim = max(abs(minimum(crt_results)), abs(maximum(crt_results)));
			# clims = (-lim, lim)
			
			subplt = heatmap(crt_results,
				c=:balance, #https://docs.juliaplots.org/latest/generated/colorschemes/
				clim=clims_hm[node], 
				colorbar=false,
				xlabel=L"t_c",
				ylabel=L"t_e",
				title="$(100*storage_penetrations[s_idx])% storage, η=$crt_η_",
				xticks=xticks_hr, yticks=xticks=xticks_hr
			)
			
			# s_idx == 1 && plot!(ylabel=L"t_e")
			s_idx == 3 && plot!(
				colorbar=true, 
				colorbar_tickfontsize=8,
				colorbar_ticks = [-500, 0, 500]
				# colorbar_title = "MEF
			)
			
			heatmap_subplts[node][idx_η][s_idx] = subplt
		end
	end
	
	# @layout l = [grid(3,3)]
	layout_hm = @layout [a{0.3w} b{0.3w} c]
	vec_hm = 
	plt_emissions_heatmap[node] = plot(
		[heatmap_subplts[node][1][si] for si in 1:length(storage_penetrations)]..., 
		# layout=Plots.grid(
			# length(η_vals), length(storage_penetrations), widths=[.29, 0.29, 0.42]
			# ), 
		layout = layout_hm,
		size=(650, 450/3*length(η_vals)), 
		bottom_margin=8Plots.pt
	)
	end
	# savefig(plt_emissions_heatmap, 
		# "../img/storage_penetration_emissions_heatmap.png")
	# plt_emissions_heatmap
		
end

# ╔═╡ 5a99452e-c842-4f6a-ab71-36df7ccefaaf
plt_emissions_heatmap[21]

# ╔═╡ 415741a9-bc70-4d2c-9210-3fce9dd1331b
plt_emissions_heatmap[23]

# ╔═╡ 52c889e4-753c-447c-a9e1-862750b3643f
nn = round.(results[idx_η][10, :, :, 1]'[16:20, 16:20])

# ╔═╡ aff80d55-df50-4d4b-aba4-e62f3c7ec10e
diag(crt_mefs)

# ╔═╡ 62f66995-bd02-4b6f-8eb8-6aeae5436713
md"""
*Question* : why don't I get *exactly* the same values from both matrices. Some diagonal elements are substantially different -- I have to make sure I did not leave any mistake pending
"""

# ╔═╡ 60c89e52-bcb4-41ed-9490-95c7ad7c2288
#data to save
hm_to_save = ([results[1][21, :, :, k]' for k in 1:length(storage_penetrations)], [results[1][23, :, :, k]' for k in 1:length(storage_penetrations)])

# ╔═╡ 59f3559b-aabe-42d7-9975-5fcc0b3de978
md"""
## Plotting the (total?) mefs per node

!!! make sure that the order of cons_time and t_display is appropriate here
"""

# ╔═╡ edabacdd-8d25-4d64-9d4a-ecf1263ac02e
md"""
## Sensitivity analysis
"""

# ╔═╡ 3c5edbc5-8fc7-4d09-98a9-85f2efb699a8
node_sens = 3 # 21

# ╔═╡ 67ad2564-fb20-4a71-a084-0145e8ed24bc
cons_time = 17;

# ╔═╡ c7deae02-3dad-4335-9449-a7e8f8bd5b4f
cons_time

# ╔═╡ acddad02-84ee-480f-a65f-716a4c34710c
begin
	MEF_bar = bar(
		[sum(results[idx_η][k, cons_time, :, 1]) for k in 1:n]
		)
	xlabel!("Node")
	# ylabel!("MEF [kgCO2/MWh]")
	ylabel!("MEF")
end

# ╔═╡ bd116217-0e1c-45a0-9239-e239dc2d639b
s_idx_ = 1 #index of storage penetration

# ╔═╡ a1e23c58-6d7b-4a69-8e33-411a7c051d37
diag(results[idx_η][node_single, :, :, s_idx_])

# ╔═╡ e5806501-044e-4667-a9b2-5d3417a7a49d
storage_penetrations[s_idx_]

# ╔═╡ 5365b74f-595f-4ade-a7af-e8dba53b84f7
md"""
Reference (as in computed) values
"""

# ╔═╡ a9b770e0-b137-40f7-b59a-35ad355b98bd
ref_mefs = results[idx_η][node_sens, :, :, s_idx_]';

# ╔═╡ 956e963a-97af-495e-9475-181322ac2e0c
ref_mefs[:, cons_time]

# ╔═╡ 4aed3df5-441b-445b-9277-a38690eb8603
begin
npoints = 10
ε = 1e-2
end;

# ╔═╡ c9b41436-e0a0-4e57-908f-b45e42122e63
md"""
The cell below perturbs demand at a given time and then we will plot the different variales as a function of the demand, trying to understand the emergence of those patterns
"""

# ╔═╡ 110f3329-c847-47f1-8427-ee959adc8745
RUN_CELL_SENSITIVITY = true

# ╔═╡ 91f7d63c-9e30-4fd4-ab39-9fbf58d101dc
begin
	if RUN_CELL_SENSITIVITY
	println("---------------------------")
	println("Running sensitivity analysis")

		
	# size of the matrices are
	# 2npoints+1: number of different values of demand for which we solve the problem
	# n: number of nodes in the graph
	# l: number of generators (= length(emissions_rates))
	# T: the time horizon

	# println("initial value of the demand:")
	# println(d_dyn[cons_time][node])
	ref_val = deepcopy(d_dyn[cons_time][node_sens])
	if ref_val > 0 
			perturb_vals = [ref_val * (1+i*ε) for i in -npoints:npoints]
			x_axis_vals = [1+i*ε for i in -npoints:npoints]
			idx_ref = npoints+1
			idx_105 = findall(x_axis_vals.==1.05)[1]
			idx_95 = findall(x_axis_vals.==.95)[1]
			idx_110 = findall(x_axis_vals .== 1.1)[1]
			idx_90 = findall(x_axis_vals .== .9)[1]
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
	mefs_sensitivity = zeros(L, T) #first index is perturbation, second is emissions time
		
	net_crt = meta[idx_η][s_idx_].net
		
	for k in 1:L
		d_crt = deepcopy(d_dyn)
		d_crt[cons_time][node_sens] = perturb_vals[k]
		opf_ = DynamicPowerManagementProblem(net_crt, d_crt)
		solve!(opf_, OPT, verbose=false)
		if opf_.problem.status != Convex.MOI.OPTIMAL
			@show opf_.problem.status
		end

		mefs_ = compute_mefs(opf_, net_crt, d_crt, emissions_rates)

		for t in 1:T
			s_sensitivity[k, :, t] = evaluate(opf_.s[t])
			g_sensitivity[k, :, t] = evaluate(opf_.g[t])
			# emissions sensitivity at 100% of the demand
			E_sensitivity[k, :, t] = evaluate(opf_.g[t]).*emissions_rates
			mefs_sensitivity[k, :] = mefs_[cons_time][node_sens, :]
		end
		# println(d_dyn[cons_time][node])
		# println(d_crt[cons_time][node])
		# println(ref_val)
	end
	end
end

# ╔═╡ 77943ac8-36fe-4a13-a36d-db957780d869
begin #E_ref is the total emissions at a given time
	E_ref = zeros(T)
	
	for t in 1:T
		E_ref[t] = evaluate(meta[idx_η][s_idx_].opf.g[t])' * emissions_rates
	end

end

# ╔═╡ 4fd2833c-6c23-4009-8734-980d3dd08c91
md"""
What is the value of emissions when there is no perturbation? 
"""

# ╔═╡ e0f5c93c-e1dd-4a9e-baf1-cbb8daf540dc
md""" *these values should be equal?* """

# ╔═╡ 6fcd6e19-58c3-462d-964f-8cd3127b47a4
sum(E_sensitivity[npoints+1, :, :], dims=1)

# ╔═╡ 2973af52-0bd0-4ba8-855d-297427627e22
E_ref[:]

# ╔═╡ b85b85d0-e1bc-4fc9-81cf-3792b55e3684
t_display=17

# ╔═╡ b674af27-307b-4dbb-8a75-a54bde1f123d
t_display

# ╔═╡ 506e9360-2c25-4ea7-830b-68b4a6bf9026
md"""
Emissions time: $t_display|
"""

# ╔═╡ 30511293-8ba5-486e-956b-e9f2a1ed0505
begin
	γ = 1e-4
	Δ = .05
	ylims = (1-Δ, 1+Δ)
	plt_s = plot(
		x_axis_vals, 
		[s_sensitivity[:, k, t_display]/(s_sensitivity[idx_ref, k, t_display]+γ) for k in 1:n], ylim=ylims
	)
	title!("Storage at time $t_display")
	xlabel!("Change in demand at node $node_sens at time $cons_time")
	ylabel!("Change in storage at all nodes at time $t_display")
	
	plt_E = plot(
		x_axis_vals, 
		[E_sensitivity[:, k, t_display]./(E_sensitivity[idx_ref, k, t_display]+γ) for k in 1:length(emissions_rates)], ylim=ylims
		)
	title!("Emissions at time $t_display")
	xlabel!("Change in demand at node $node_sens at time $cons_time")
	ylabel!("Change in emissions at all generators at time $t_display")
	
	plt_g = plot(
		x_axis_vals, 
		[g_sensitivity[:, k, t_display]./(g_sensitivity[idx_ref, k, t_display]+γ) for k in 1:length(emissions_rates)], ylim=ylims
		)
	title!("Generators at time $t_display")
	xlabel!("Change in demand at node $node_sens at time $cons_time")
	ylabel!("Change in generation at all generators at time $t_display")

	norm_E = sum(E_sensitivity[:, :, t_display], dims=2)./sum(E_sensitivity[idx_ref, :, t_display])
	plt_E_tot = plot(
		x_axis_vals, norm_E
		, ylim=(0.98, 1.08)
		)
	# xlabel!("Change in demand at node $node_sens at time $cons_time")
	# ylabel!("Change in total emissions")
	xlabel!(L"\Delta d/d")
	ylabel!(L"\Delta E/E")
	
	#adding the theoretical curve for the sensitivity
	E_th = (
		sum(E_sensitivity[idx_ref, :, t_display]) .+ (perturb_vals.-ref_val) .* mefs_sensitivity[idx_ref, t_display]
		)./sum(E_sensitivity[idx_ref, :, t_display])

	E_th_105 = (
		sum(E_sensitivity[idx_105, :, t_display]) .+ (perturb_vals.-ref_val * x_axis_vals[idx_105]) .* mefs_sensitivity[idx_105, t_display]
		)./sum(E_sensitivity[idx_ref, :, t_display])

E_th_95 = (
		sum(E_sensitivity[idx_95, :, t_display]) .+ (perturb_vals.-ref_val * x_axis_vals[idx_95]) .* mefs_sensitivity[idx_95, t_display]
		)./sum(E_sensitivity[idx_ref, :, t_display])
	E_th_110 = (
		sum(E_sensitivity[idx_110, :, t_display]) .+ (perturb_vals.-ref_val * x_axis_vals[idx_110]) .* mefs_sensitivity[idx_110, t_display]
		)./sum(E_sensitivity[idx_ref, :, t_display])

	E_th_90 = (
		sum(E_sensitivity[idx_90, :, t_display]) .+ (perturb_vals.-ref_val * x_axis_vals[idx_90]) .* mefs_sensitivity[idx_90, t_display]
		)./sum(E_sensitivity[idx_ref, :, t_display])
	
	plot!(x_axis_vals, E_th, ls=:dash, c=:orange)
	plot!(x_axis_vals, E_th_105, ls=:dash, c=:green)
	plot!(x_axis_vals, E_th_95, ls=:dash, c=:firebrick)
	plot!(x_axis_vals, E_th_110, ls=:dash, c=:violetred)

	scatter!([x_axis_vals[idx_ref]], [norm_E[idx_ref]], c=:orange)
	scatter!([x_axis_vals[idx_105]], [norm_E[idx_105]], c=:green)
	scatter!([x_axis_vals[idx_95]], [norm_E[idx_95]], c=:firebrick)
	scatter!([x_axis_vals[idx_110]], [norm_E[idx_110]], c=:violetred)
	title!("Total emissions at time $t_display")
	
	@show ref_mefs[t_display, cons_time]
	@show t_display
	@show cons_time
	@show ref_val
	
	plot([plt_s, plt_E, plt_g, plt_E_tot]..., size = (650, 650), lw = 3)
	
end

# ╔═╡ 49dc5403-f19b-458d-b9d5-f2baf2e68d17
# data to save
begin
x_sens = x_axis_vals
y_sens_exp = vec(sum(E_sensitivity[:, :, t_display], dims=2)./sum(E_sensitivity[idx_ref, :, t_display]))
y_sens_100 = E_th
y_sens_105 = E_th_105
y_sens_110 = E_th_110
y_sens_90 = E_th_90
y_sens_95 = E_th_95

data_sensitivity = (
	x_sens, y_sens_exp, y_sens_90, y_sens_95, y_sens_100, y_sens_105, y_sens_110
)
end;

# ╔═╡ a1d7b77f-14e4-4a8a-806f-ebf70d4f1e3c
# change in generators, renewable vs non renewable
begin
rel_g_not_renew = sum(g_sensitivity[:, 1:l_no_renew, t_display], dims=2)./sum(g_sensitivity[idx_ref, 1:l_no_renew, t_display].*(1+γ));
rel_g_renew = sum(g_sensitivity[:, l_no_renew+1:end, t_display], dims=2)./sum(g_sensitivity[idx_ref, l_no_renew+1:end, t_display].*(1+γ));
	
plt_curt = plot(
	x_axis_vals, 
	rel_g_not_renew, 
	label="Non renewable"
		)
	plot!(
		x_axis_vals, 
		rel_g_renew,
		label="Renewable"
	)
	title!("Generators at time $t_display")
	xlabel!("Change in demand at node $node_sens at time $cons_time")
	ylabel!("Change in generation at all generators at time $t_display")

end

# ╔═╡ c291accc-8774-4319-a7b5-a5129e699ec0
begin
plot_gt = plot()
for gt in unique(tags)
	g_idx = findall(tags.==gt)
	rel_g_crt = sum(g_sensitivity[:, g_idx, t_display], dims=2)./sum(g_sensitivity[idx_ref, g_idx, t_display].*(1+γ));
	plot!(
		x_axis_vals, 
		rel_g_crt,
		label=gt
	)
end


	title!("Generators at time $t_display")
	xlabel!("Change in demand at node $node_sens at time $cons_time")
	ylabel!("Change in generation at all generators at time $t_display")
	plot_gt
end

# ╔═╡ e9e1f2b7-bbbc-4f7c-9997-1b3ee1796c14
#plot of total emissions in theory (predicted by total mefs) vs practice
begin
	# plot of total Emissions over T
	plt_E_tot_T = plot(
		x_axis_vals, vec(sum(E_sensitivity[:, :, :],
		dims=(2,3)))./vec(sum(E_sensitivity[idx_ref, :, :], dims = (1, 2))), 
		ylim = (.99, 1.01)
		)
	xlabel!(L"\Delta d/d")
	ylabel!(L"\Delta E/E")
		
		#adding the theoretical curve for the sensitivity
		E_th_T = (
			sum(E_sensitivity[idx_ref, :, :], dims = (1, 2)) .+
			(perturb_vals.-ref_val) .* ref_mefs[t_display, cons_time]
		)./sum(E_sensitivity[idx_ref, :, :], dims = (1, 2))
		plot!(x_axis_vals, vec(E_th_T), ls=:dash)
	
end

# ╔═╡ d8d1fb74-0018-4685-a283-e768ae877fe4
md"""
## Complete figure
"""

# ╔═╡ e14a1d29-477d-4ed5-908f-f436f00b7fa2
begin

lw = 2
fs = 8
l_top = @layout [
		a{.6w} b{.4w}
]
plt_top = plot(
		MEF_bar, plt_E_tot, #plt_E_tot_T,
		layout = l_top, size=(800, 200), lw=2, 
		legend=:outertopright, title = ["($i)" for j in 1:1, i in 1:10], titleloc = :right
	)

l_ = @layout[
	c{.21w} d{.21w} e{.21w} f
]

plts_21 = [
	subplots_mef_storage[21], heatmap_subplts[21][1][1], heatmap_subplts[21][1][2], heatmap_subplts[21][1][3]
]
plt_middle = plot(
		plts_21...,
		layout = l_, size=(800, 150), lw=lw, fs = fs, 
		legend=:outertopright, title = ["($i)" for j in 1:1, i in 3:10], titleloc = :right
	)

l_ = @layout[
	g{.21w} h{.21w} i{.21w} j
]
plts_23 = [
	subplots_mef_storage[23], heatmap_subplts[23][1][1], heatmap_subplts[23][1][2], heatmap_subplts[23][1][3]
]
plt_bottom = plot(
		plts_23..., 
		layout = l_, size=(800, 150), lw=lw, fs=fs,
		legend=:outertopright, title = ["($i)" for j in 1:1, i in 7:10], titleloc = :right
	)

l_fig = @layout[
	g{.5h}
	h{.25h}
	i
]
Fig = plot(
	[plt_top, plt_middle, plt_bottom]..., layout = l_fig, size=(650, 400),
	lw=lw, fs=fs,
)
	
#save
savefig(Fig, 
	"../img/Fig_storage.pdf")
Fig
end

# ╔═╡ 062ffc7d-86da-48b0-bb63-8aa16c4bb5b7


# ╔═╡ ec65009f-cda6-4874-be4a-2326c1c46300
# subplots_mef_storage[21], plt_emissions_heatmap[21],
# subplots_mef_storage[23], plt_emissions_heatmap[23],

# ╔═╡ 28ad54e9-2cee-4f99-93e9-40f23471ed94
md"""
## HDF5
"""

# ╔═╡ 2407370d-6cdb-4c36-a8d1-fa8c1e515c61
let
fname = "/Users/lucasfuentes/sensitivity/results/fig2_data_sensitivity_node_$(node_sens)"
fid = h5open(fname, "w")
fid["x"] = data_sensitivity[1]
fid["y_exp"] = data_sensitivity[2]
fid["y_90"] = data_sensitivity[3]
fid["y_95"] = data_sensitivity[4]
fid["y_100"] = data_sensitivity[5]
fid["y_105"] = data_sensitivity[6]
fid["y_110"] = data_sensitivity[7]
close(fid)
end

# ╔═╡ bec5bfdc-5d4a-41c9-87ea-4429cc12d4ae
begin
fname = "/Users/lucasfuentes/sensitivity/results/fig2_data_1"
fid = h5open(fname, "w")
	
fid["MEFs"] = MEFS_to_save
fid["MEF_vs_t_x"] = data_total_mef[1]
fid["MEF_vs_t_y1"] = data_total_mef[2][1]
fid["MEF_vs_t_y2"] = data_total_mef[2][2]
fid["storage_pens"] = storage_penetrations

for k in 1:length(storage_penetrations)
	for ni in nodes_heatmaps
		fid["hm_$(ni)_$k"] = Array(results[idx_η][ni, :, :, k]')
	end
end
	
close(fid)
end

# ╔═╡ a2a1a718-e0e5-4731-ba9b-07f18beee578
heatmap(results[idx_η][21, :, :, 1], c=:balance, clim=(-700, 700))

# ╔═╡ e94c4b92-ceec-412f-adfb-6e9a7344ca39
md"""
## Future work

How do MEFs behave and converge to a given value?
"""

# ╔═╡ c76f2ebe-ce41-47fd-b31a-2851aca53567
let
	plot()
	for ct in 1:T
		plot!(
			storage_penetrations, 
			[sum(results[idx_η][node, ct, :, j]) for j in 1:length(storage_penetrations)], 
			lw=3, ls=:dash, markershape=:circle
			)
	end
	plot!()
	xlabel!("storage penetration")
	ylabel!("Total MEF")
end

# ╔═╡ 25063860-6109-46e7-9dd5-a7fc0c12159e
let
	plot()
	for t in 1:T
		plot!(storage_penetrations, [results[idx_η][node, cons_time, t, j] for j in 1:length(storage_penetrations)], lw=3, ls=:dash, markershape=:circle)
	end
	plot!(size=(500, 500))
	xlabel!("storage penetration")
	ylabel!("MEF")
	title!("MEF at node $node_matrix and consumption time $cons_time for different emsissions times")
end

# ╔═╡ 98876bb5-1370-453f-80fa-0557da20d9ff
md"""# Old"""

# ╔═╡ 16094734-d7a1-4bb2-990b-8de1c367b134
md"""## How does charging efficiency affect mefs"""

# ╔═╡ 4a8c7752-6de8-4ea7-aafe-702b17507185
storage_pen = s_rel

# ╔═╡ 19b6abf5-6951-4492-8f17-f76df29f9289
RUN_BIG_CELL1 = false

# ╔═╡ 47e2e177-3701-471f-ae3c-38276ec69370
begin
	if RUN_BIG_CELL1
	println("----------------------------------------------------------")
	println("Recomputing results for different η")
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

		if opf_dyn.problem.status != Convex.MOI.OPTIMAL
			@show opf_dyn.problem.status
		end

		# Compute MEFs
		mefs = compute_mefs(opf_dyn, net_dyn, d_dyn, emissions_rates)
		println("...MEFs computed")
		for ind_t in 1:length(mef_times)
			results_η[:, :, ind_t, ind_η] .= mefs[ind_t]
		end

		meta_η[ind_η] = (opf=opf_dyn, net=net_dyn)
	end
	end
		

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

# ╔═╡ Cell order:
# ╟─c39005df-61e0-4c08-8321-49cc5fe71ef3
# ╠═0510ec8c-2f1b-4704-bb59-cd8a67ef0dc5
# ╟─44275f74-7e7c-48d5-80a0-0f24609ef327
# ╠═db59921e-e998-11eb-0307-e396d43191b5
# ╠═571d1cff-7311-4db8-8ac3-9e10afefaf18
# ╠═9af03cdd-8cb8-4221-b962-05e6ae1634cc
# ╠═0aac9a3f-a477-4095-9be1-f4babe1e2803
# ╠═a32d6a56-8da8-44b0-b659-21030692630a
# ╠═257a6f74-d3c3-42eb-8076-80d26cf164ca
# ╠═113e61a9-3b21-48d0-9854-a2fcce904e8a
# ╟─9bd515d4-c7aa-4a3d-a4fb-28686290a134
# ╟─1bd72281-4a7f-44f4-974d-632e9d0aaf28
# ╠═0c786da1-7f44-40af-b6d6-e0d6db2242b2
# ╠═5b80ca83-0719-437f-9e51-38f2bed02fb4
# ╠═cfcba5ad-e516-4223-860e-b1f18a6449ba
# ╟─75dfaefd-abec-47e2-acc3-c0ff3a01048e
# ╟─7d738cd3-1edf-4348-aa16-5ade0f87ae6f
# ╟─63db8d8b-8229-47e0-a309-75813120aaab
# ╠═35dccacb-8780-4f37-8dc5-29064a84d49c
# ╠═f999d732-14b3-4ac5-b803-3df7a96ef898
# ╠═15245a9a-8a42-4ec1-ba2c-9bf95d6c85fa
# ╠═0d42b50d-993e-4eea-9025-2b7479bb3b0e
# ╠═83e265f4-0a62-4c53-9203-6a7768c6220e
# ╠═23690382-3d30-46e3-b26a-a30875be78ec
# ╠═c92e091e-1111-4bda-8d76-543ad1dd8121
# ╟─379598be-2c42-4c29-8a24-91b74592da0f
# ╠═39cdac4b-87bb-441c-99e3-402b40bef7d3
# ╟─34d4bd62-6be2-4089-8caa-1a8715bee433
# ╠═b7476391-30b9-4817-babf-7c9078531ee7
# ╠═c82ef027-740a-49b1-93d2-1554c411a896
# ╟─9873abd2-fc2f-4a35-8689-51d558debbf0
# ╠═0239e1da-caf5-4593-af1b-5d1e8d2f2b3e
# ╠═9e5a1672-d452-42f5-ba5a-a2fa0b1eaada
# ╠═522e95d5-15a8-47ec-a79c-c4cc17cf86fd
# ╠═bfa4a8f9-cfcd-4e22-b2dc-751226f3a73c
# ╟─a8ccbc8e-24e6-4214-a179-4edf3cf26dad
# ╟─496135ec-f720-4d43-8239-d75cc7616f58
# ╠═806819d5-7b40-4ca4-aa1e-f1cf0a9a7f3f
# ╟─c8f71644-9371-443b-b976-1734cc7ae583
# ╠═57a37eb0-a547-428c-9f8c-e5f3e30f5260
# ╠═c5f31bf6-5cca-4b0e-86f0-b44321fec874
# ╟─af6d4ef0-f59f-42be-af36-7cf447478e4c
# ╠═06498d79-f20a-4691-8c72-8f1f962e6a6f
# ╟─49aea0fc-1e0b-4e9c-8d7b-5d3727e95e1f
# ╠═4d056537-3ec5-47d1-9f51-a603b8ae91ff
# ╠═871ca595-9d9e-4f06-9119-c2b07bfdeb04
# ╠═0d99fc04-0353-4170-a23e-f21460ceaf7e
# ╟─38e73213-d399-43e5-80e8-851b6cf3299d
# ╟─f15393ae-5c04-473e-ae28-ead5554896a3
# ╠═1bd70b8a-90ff-45f3-8ef4-dee0340a02f0
# ╠═a1517984-dfdb-4be7-94cd-16fbf2dbc453
# ╠═58368b3b-c042-4f37-81f0-f8dd0388e5aa
# ╠═0ffb55c7-1371-45ac-9823-406eb86bd8e4
# ╠═794b4672-0d30-4df6-9fa1-2d68d34c2502
# ╠═2bde6957-46e5-4697-9451-40b135f8523c
# ╠═a4e2b01f-6ea9-4129-954b-3404b0e16e8f
# ╠═cf7430e5-c46b-4db7-adf0-8058f9ac8640
# ╟─30ec492c-8d21-43f6-bb09-32810494f21e
# ╠═856a78d9-7b4c-453b-b73b-c81eee014e52
# ╠═98a0d7c5-b1a8-4ebe-bb73-7ca88b475592
# ╠═71cd6842-7df3-4b0a-ae23-c49404ddf523
# ╟─bbbb358c-e645-4989-bed3-73d9217f7447
# ╠═6f08828b-4c4b-4f50-bd40-35805a37aae0
# ╟─15293269-580f-4251-be36-4be6ba8c5a46
# ╟─cd5fe410-6264-4381-b19f-31d050bc3930
# ╠═0740dc70-a532-4818-b09d-b3b8d60fa6ba
# ╠═f26187fb-d4b2-4f0d-8a80-5d831e0de6c3
# ╠═4925c50b-12c0-4217-94de-bdcc72c01ccf
# ╟─4d9a4b36-6b3d-4836-8501-7f46cd7ab5cc
# ╠═75d956fc-bcf6-40fe-acd5-b5eef0fc7902
# ╟─c6ee857d-8019-4c4f-bb07-a370a88ea3cf
# ╠═11b97da8-0e0d-48dd-a442-4ffa655bec61
# ╠═6186798f-6711-4222-94bb-f53b2d0fad7d
# ╠═d27ef0d8-70b2-4897-9000-8fa70b1862fc
# ╠═420919bc-f217-4357-bdbb-83e25e83ba56
# ╠═5f77f4a9-ff5a-4515-9016-bf36571225c7
# ╠═1da34733-fee3-42e1-b5e0-cac3f5f196c9
# ╠═57643580-0a1d-4ad5-ba21-533dcbd73c2f
# ╠═c7deae02-3dad-4335-9449-a7e8f8bd5b4f
# ╠═b674af27-307b-4dbb-8a75-a54bde1f123d
# ╠═f7e0d09c-40bf-4936-987a-a3bcadae5487
# ╠═5a99452e-c842-4f6a-ab71-36df7ccefaaf
# ╠═415741a9-bc70-4d2c-9210-3fce9dd1331b
# ╠═52c889e4-753c-447c-a9e1-862750b3643f
# ╠═aff80d55-df50-4d4b-aba4-e62f3c7ec10e
# ╠═a1e23c58-6d7b-4a69-8e33-411a7c051d37
# ╠═62f66995-bd02-4b6f-8eb8-6aeae5436713
# ╠═60c89e52-bcb4-41ed-9490-95c7ad7c2288
# ╟─59f3559b-aabe-42d7-9975-5fcc0b3de978
# ╠═acddad02-84ee-480f-a65f-716a4c34710c
# ╟─edabacdd-8d25-4d64-9d4a-ecf1263ac02e
# ╠═3c5edbc5-8fc7-4d09-98a9-85f2efb699a8
# ╠═67ad2564-fb20-4a71-a084-0145e8ed24bc
# ╠═bd116217-0e1c-45a0-9239-e239dc2d639b
# ╠═e5806501-044e-4667-a9b2-5d3417a7a49d
# ╟─5365b74f-595f-4ade-a7af-e8dba53b84f7
# ╠═a9b770e0-b137-40f7-b59a-35ad355b98bd
# ╠═956e963a-97af-495e-9475-181322ac2e0c
# ╠═4aed3df5-441b-445b-9277-a38690eb8603
# ╟─c9b41436-e0a0-4e57-908f-b45e42122e63
# ╠═110f3329-c847-47f1-8427-ee959adc8745
# ╠═91f7d63c-9e30-4fd4-ab39-9fbf58d101dc
# ╠═77943ac8-36fe-4a13-a36d-db957780d869
# ╟─4fd2833c-6c23-4009-8734-980d3dd08c91
# ╟─e0f5c93c-e1dd-4a9e-baf1-cbb8daf540dc
# ╠═6fcd6e19-58c3-462d-964f-8cd3127b47a4
# ╠═2973af52-0bd0-4ba8-855d-297427627e22
# ╟─506e9360-2c25-4ea7-830b-68b4a6bf9026
# ╠═b85b85d0-e1bc-4fc9-81cf-3792b55e3684
# ╠═30511293-8ba5-486e-956b-e9f2a1ed0505
# ╠═49dc5403-f19b-458d-b9d5-f2baf2e68d17
# ╠═a1d7b77f-14e4-4a8a-806f-ebf70d4f1e3c
# ╠═c291accc-8774-4319-a7b5-a5129e699ec0
# ╠═e9e1f2b7-bbbc-4f7c-9997-1b3ee1796c14
# ╟─d8d1fb74-0018-4685-a283-e768ae877fe4
# ╠═e14a1d29-477d-4ed5-908f-f436f00b7fa2
# ╠═062ffc7d-86da-48b0-bb63-8aa16c4bb5b7
# ╠═ec65009f-cda6-4874-be4a-2326c1c46300
# ╟─28ad54e9-2cee-4f99-93e9-40f23471ed94
# ╠═2407370d-6cdb-4c36-a8d1-fa8c1e515c61
# ╠═bec5bfdc-5d4a-41c9-87ea-4429cc12d4ae
# ╠═a2a1a718-e0e5-4731-ba9b-07f18beee578
# ╟─e94c4b92-ceec-412f-adfb-6e9a7344ca39
# ╟─c76f2ebe-ce41-47fd-b31a-2851aca53567
# ╟─25063860-6109-46e7-9dd5-a7fc0c12159e
# ╟─98876bb5-1370-453f-80fa-0557da20d9ff
# ╠═16094734-d7a1-4bb2-990b-8de1c367b134
# ╠═4a8c7752-6de8-4ea7-aafe-702b17507185
# ╠═19b6abf5-6951-4492-8f17-f76df29f9289
# ╠═47e2e177-3701-471f-ae3c-38276ec69370
# ╠═457c1959-94fa-4267-8645-3ed1409cd0a0
# ╟─a7e75e49-5031-4cc4-b96e-6227277ec3ba
# ╟─d9617524-76c3-447d-9b94-0a690f83a7b9
# ╟─fb43471f-6aed-4fd6-a9e0-b165f6c77003
# ╟─9aded04b-e55f-4ebd-97c4-90c3adf62547
