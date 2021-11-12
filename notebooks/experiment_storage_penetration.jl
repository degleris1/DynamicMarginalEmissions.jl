### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

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

# ╔═╡ c39005df-61e0-4c08-8321-49cc5fe71ef3
md"""
## Description
"""

# ╔═╡ 5867a4eb-470a-4a8a-84e1-6f150de1dcde
md"""
## ToDo
- Make sure there is renewable curtailment
"""

# ╔═╡ 751a51bb-97c6-4608-8195-7ed465b9eb7c
md"""
## Things to figure out
- Storage patterns
- Why is a given timestep more significant than other timesteps
- Be careful about the value of the time, when storage is and is not...
- why do I need to increase `p` so significantly to get something that is feasible? 
"""

# ╔═╡ 0510ec8c-2f1b-4704-bb59-cd8a67ef0dc5
md"""
## *TODO*
- Update pmax: I think this is the big one: how do you allocate the pmax so that you have interesting phenomena
- update renewable penetration
- PICK REALISTIC VALUES FOR THE SYSTEM and scale the values properly. This is what you need!
- handle renewable location: done
- handle storage location
- make sure we know exactly how to flip the mef matrix, and we avoid confusion
- add colors for the flows, to see which one is congested vs non congeste
- play with the location of renewables to get something more interesting
- make sure renewable generation is properly allocated
- another way to generate interesting results is to increase demand instead of changing line capacities
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

# ╔═╡ 75dfaefd-abec-47e2-acc3-c0ff3a01048e
md"""
### Network
"""

# ╔═╡ f70aa9bd-4956-4e5a-8b3d-a1ad5b0b54d0
md"""
*Change the pmax and how it is handled*

We probably want a few lines that are low capacity and many that are high capacity. 
"""

# ╔═╡ 7d738cd3-1edf-4348-aa16-5ade0f87ae6f
#define which nodes have renewable generation

# renew_nodes = [1, 8, 15, 26]
# renew_nodes = [4, 12, 16, 28]
# renew_nodes = [15, 18, 19, 20]
renew_nodes = [12, 14, 15, 18]

# ╔═╡ 63db8d8b-8229-47e0-a309-75813120aaab
placed_gen_nodes = [9, 11, 6, 8, 28, 7] # 2 negative mefs
# placed_gen_nodes = [1, 2, 3, 4, 5, 6]
# placed_gen_nodes = [19, 20, 23, 6, 9, 10]

# ╔═╡ 1ff66201-1d03-4c2e-95e5-d3387b0ddd57
md"""
isolated edges - they probably need high capacities
- (25, 26): 34
- (27, 29): 39
- (27, 30): 40
- (9, 11): 14
- (12, 13): 16
"""

# ╔═╡ 35dccacb-8780-4f37-8dc5-29064a84d49c
begin
	# add_low = [1, 39, 33, 6, 9, 36]
	add_low = [1, 6, 13, 21, 25, 26, 33, 30, 32, 9, 3]
	add_med = [35, 7, 12, 38, 36, 11]
	add_zero = []

	
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
	# n_low_capacity = 10
	capacities = ones(size(net.A)[2]) * .6
	# capacities[1:n_low_capacity] .= .4 # .4
	# capacities = Random.shuffle(capacities)
	capacities[[29, 31]] .= .4 # minimum necessary, edges around gen 22
	
	capacities[add_low] .= level_low

	capacities[add_med] .= level_med
	capacities[add_zero] .= δ


	
	# idx_low = [5, 19, 23, 26, 29, 31, 32]
	# capacities[idx_low] .= .4
	# add_low = [1, 2, 3, 4, 6, 7, 8, 11, 12, 15, 9, 20, 21, 30, 33, 37, 28]
	# capacities[add_low] .= .4
	
# 	idx_high = [14, 16, 34, 39, 40] # important edges
# 	capacities[idx_high] .= 1.
# 	add_high = [1, 4, 20, 25, 30, 32, 6, 10, 12, 18, 15, 8, 9, 11, 23, 24, 33, 35, 28, 36, 31]
# 	capacities[add_high] .= 1.
	
	# n_high_capacity = 15
	# capacities = ones(size(net.A)[2])*.4
	# capacities[1:n_high_capacity] .= .9
	# capacities = Random.shuffle(capacities)
	
	# add_low_capacity = 4
	# new_low = rand(1:length(capacities), add_low_capacity)
	# capacities[new_low] .= .4
	net.pmax = net.pmax .* capacities 
	
	net.gmax *= 2.
	
	# Remove lines - Lucas: not sure whether we need to remove the lines - TBD
	# net.pmax[[2, 3, 5, 11, 12, 20, 22, 26, 27, 30|, 39]] .= δ #used to be *= δ
	
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

# ╔═╡ 402a32b8-6b86-4701-a84b-eeb224e64a81
md""" you could also walk mthrough the nodes and make sure that not all edges getting to a node are low capacity?"""

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
	
	plt2 = plot(renew_data[:, renew_labels.=="solar"], lw=2, palette=:Reds, title="Renewable Data", xlabel="Hour", ylabel="Energy [?]", xticks=xticks_hr)
	plot!(renew_data[:, renew_labels.=="wind"], lw=2, palette=:Greens, title="Renewable Data", xlabel="Hour", ylabel="Energy [?]", xticks=xticks_hr)
	
	plt_time_series = plot(plt1, plt2, layout=(2, 1), size=(600, 300))
	
	savefig(plt_time_series, "../img/storage_penetration_time_series.png")
	plt_time_series
end

# ╔═╡ c82ef027-740a-49b1-93d2-1554c411a896
renewable_penetration = 0.1 # set properly!

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

# ╔═╡ 73875562-ad97-4880-8694-0de3d088372a
renewable_penetration

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

# ╔═╡ 2b634b6d-dc54-4ccb-8b41-7a881a018184
renew_gmax

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

# ╔═╡ e8ee5cbb-4afc-4737-b006-90071f6138cd
begin	
	# [:NG, :NG_NEW, :NUC, :PEAK, :COL, :COL], roughly
	emissions_rates = [0.8e3, 0.7e3, 0.9e3, 1.4e3, 1.8e3, 1.8e3]
	emissions_rates = [emissions_rates; zeros(length(renew_nodes))]
end;

# ╔═╡ c8f71644-9371-443b-b976-1734cc7ae583
md"""
## What is the impact of charging/discharging efficiency on storage dynamics?
"""

# ╔═╡ 57a37eb0-a547-428c-9f8c-e5f3e30f5260
c_rate = .25 #charging rate

# ╔═╡ 491da4e6-03e6-49d3-907d-43ddcffdfb42
s_rel = .0 #storage penetration

# ╔═╡ 355ebed2-42e6-41bb-b1db-a72d1aaae56f
mef_times = 1:24

# ╔═╡ a9d16b8a-a488-4a9e-bb6f-4ee05873e86b
begin
	η_c = 1
	η_d = η_c
end;

# ╔═╡ 6cc82262-cc73-4483-b97e-664d5093d69c
@show η_c, η_d

# ╔═╡ 1098efd7-8477-4749-8ed5-b62044bebdc7
node_storage = 0

# ╔═╡ c5f31bf6-5cca-4b0e-86f0-b44321fec874
begin
	println("Solving first problem")
	mefs = zeros(
			n, T, length(mef_times)
		)
	
	# Construct dynamic network
	C = sum(d_dyn) .* (s_rel + δ)/T
	
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
node_single = 30

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
mef_bar_plt
end

# ╔═╡ 0ace8911-c28f-47cd-8117-581b0cc5e1d7
mef_bar_plt

# ╔═╡ 8eb901fe-aa93-4d72-bc9d-194e199db013
diag(crt_mefs)

# ╔═╡ 5ae84f68-a0e6-47bb-9206-cd1384200581
md""" adding colors to the edges based on their flow
I can create a `cgrad` or a `palette` from julia colorschemes
"""

# ╔═╡ e33f9e0f-7c4c-465e-9793-ed790d3256b2
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

# ╔═╡ 3f9eb091-059c-44a5-9b50-ae3cabe24060
md"""
## How does charging efficiency affect MEFs
"""

# ╔═╡ 4a8c7752-6de8-4ea7-aafe-702b17507185
storage_pen = s_rel

# ╔═╡ 6e6b15b1-7685-4a20-9d94-dd703caa2fe9
η_vals = [0.95, 0.99, 1.] 

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

# ╔═╡ 856a78d9-7b4c-453b-b73b-c81eee014e52
RUN_BIG_CELL2 = false

# ╔═╡ 98a0d7c5-b1a8-4ebe-bb73-7ca88b475592
storage_penetrations = [0.0, 0.05, 0.10]

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
			C = sum(d_dyn) * (s_rel + δ)
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
The cursor below allows choosing for a specific value of charge/discharge efficiency
"""

# ╔═╡ 19f4e0cc-0c93-42dc-8ee4-17f52d4e5e90

# @bind idx_η Slider(1:1:length(η_vals))


# ╔═╡ 4925c50b-12c0-4217-94de-bdcc72c01ccf
idx_η = 3

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

# ╔═╡ 1da34733-fee3-42e1-b5e0-cac3f5f196c9
@bind node_matrix Slider(1:1:n)

# ╔═╡ 6fc320b1-b60d-4f49-89ab-bf029ead6b55
md"""
MEF as a function of consumption and emission times for node $node_matrix
"""

# ╔═╡ 32cd894b-ee5d-44b7-8983-82a4b72524a8
begin
plot([d_dyn[i][node_matrix] for i in 1:T], size=(300, 300))
xlabel!("Demand at node $node_matrix")
ylabel!("Time")
end

# ╔═╡ f7e0d09c-40bf-4936-987a-a3bcadae5487
let
	#node = interesting_nodes[highlighted_node]
	node = node_matrix# we focus on a single node
	
	plots = []
	heatmap_subplts = []
	for idx_η in 1:length(η_vals)
		crt_η_ = η_vals[idx_η]
	for s_idx in 1:length(storage_penetrations)
		
		crt_results = results[idx_η][node, :, :, s_idx]'
		# @show maximum(crt_results)
		lim = max(abs(minimum(crt_results)), abs(maximum(crt_results)));
		clims = (-lim, lim)
		subplt = heatmap(crt_results,
			c=:balance, #https://docs.juliaplots.org/latest/generated/colorschemes/
			clim=clims, 
			colorbar=false,
			xlabel="Consumption Time",
			title="$(100*storage_penetrations[s_idx])% storage, η=$crt_η_", 					xticks=xticks_hr, yticks=xticks=xticks_hr
		)
		
		s_idx == 1 && plot!(ylabel="Emissions Time")
		s_idx == 3 && plot!(colorbar=true)
		
		push!(heatmap_subplts, subplt)
	end
	

	end
	# @layout l = [grid(3,3)]
	plt_emissions_heatmap = plot(
		heatmap_subplts..., 
		layout=Plots.grid(
			length(η_vals), length(storage_penetrations), widths=[.29, 0.29, 0.42]
			), 
		size=(650, 650/3*length(η_vals)), 
		bottom_margin=8Plots.pt
	)
	# savefig(plt_emissions_heatmap, 
		# "../img/storage_penetration_emissions_heatmap.png")
	# plt_emissions_heatmap
	
	
	
	
end

# ╔═╡ 52c889e4-753c-447c-a9e1-862750b3643f
nn = round.(results[idx_η][10, :, :, 1]'[16:20, 16:20])

# ╔═╡ aff80d55-df50-4d4b-aba4-e62f3c7ec10e
crt_mefs

# ╔═╡ 62f66995-bd02-4b6f-8eb8-6aeae5436713
md"""
*Question* : why don't I get *exactly* the same values from both matrices. Some diagonal elements are substantially different -- I have to make sure I did not leave any mistake pending
"""

# ╔═╡ 498c3c62-797c-4ea9-b8cb-61f2a91b269e
node_matrix


# ╔═╡ 57ea231c-86e7-44d3-9093-836c9820f88c
node_single

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
node = node_matrix

# ╔═╡ 67ad2564-fb20-4a71-a084-0145e8ed24bc
@bind cons_time Slider(1:1:T)

# ╔═╡ c7deae02-3dad-4335-9449-a7e8f8bd5b4f
cons_time

# ╔═╡ acddad02-84ee-480f-a65f-716a4c34710c
begin
	bar(
		[sum(results[idx_η][k, cons_time, :, 1]) for k in 1:n]
		)
	xlabel!("Node")
	ylabel!("Total MEF")
end

# ╔═╡ b1f2aba6-5b6b-443c-84ab-21c4d2017a07
md"""Cons time = $cons_time"""

# ╔═╡ bd116217-0e1c-45a0-9239-e239dc2d639b
s_idx_ = 1

# ╔═╡ a1e23c58-6d7b-4a69-8e33-411a7c051d37
results[idx_η][node_matrix, :, :, s_idx_]

# ╔═╡ e5806501-044e-4667-a9b2-5d3417a7a49d
storage_penetrations[s_idx_]

# ╔═╡ 5365b74f-595f-4ade-a7af-e8dba53b84f7
md"""
Reference (as in computed) values
"""

# ╔═╡ a9b770e0-b137-40f7-b59a-35ad355b98bd
ref_mefs = results[idx_η][node, :, :, s_idx_]';

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
RUN_CELL_SENSITIVITY = false

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
	ref_val = deepcopy(d_dyn[cons_time][node])
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
		
	net_crt = meta[idx_η][s_idx_].net
		
	for k in 1:L
		d_crt = deepcopy(d_dyn)
		d_crt[cons_time][node] = perturb_vals[k]
		opf_ = DynamicPowerManagementProblem(net_crt, d_crt)
		solve!(opf_, OPT, verbose=false)
		if opf_.problem.status != Convex.MOI.OPTIMAL
			@show opf_.problem.status
		end

		for t in 1:T
			s_sensitivity[k, :, t] = evaluate(opf_.s[t])
			g_sensitivity[k, :, t] = evaluate(opf_.g[t])
			E_sensitivity[k, :, t] = evaluate(opf_.g[t]).*emissions_rates
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
@bind t_display Slider(1:1:T)

# ╔═╡ b674af27-307b-4dbb-8a75-a54bde1f123d
t_display

# ╔═╡ 506e9360-2c25-4ea7-830b-68b4a6bf9026
md"""
Emissions time: $t_display|
"""

# ╔═╡ 30511293-8ba5-486e-956b-e9f2a1ed0505
let
	γ = 1e-4
	Δ = .02
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
		sum(E_sensitivity[idx_ref, :, t_display]) .+ (perturb_vals.-ref_val) .* ref_mefs[t_display, cons_time]
		)./sum(E_sensitivity[idx_ref, :, t_display])
	plot!(x_axis_vals, E_th, ls=:dash)
	title!("Total emissions at time $t_display")
	
	@show ref_mefs[t_display, cons_time]
	@show t_display
	@show cons_time
	@show ref_val
	
	plot([plt_s, plt_E, plt_g, plt_E_tot]..., size = (650, 650), lw = 3)
	
end

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

# ╔═╡ 16094734-d7a1-4bb2-990b-8de1c367b134


# ╔═╡ Cell order:
# ╟─c39005df-61e0-4c08-8321-49cc5fe71ef3
# ╠═5867a4eb-470a-4a8a-84e1-6f150de1dcde
# ╠═751a51bb-97c6-4608-8195-7ed465b9eb7c
# ╠═0510ec8c-2f1b-4704-bb59-cd8a67ef0dc5
# ╟─44275f74-7e7c-48d5-80a0-0f24609ef327
# ╠═db59921e-e998-11eb-0307-e396d43191b5
# ╠═0aac9a3f-a477-4095-9be1-f4babe1e2803
# ╠═a32d6a56-8da8-44b0-b659-21030692630a
# ╠═257a6f74-d3c3-42eb-8076-80d26cf164ca
# ╠═113e61a9-3b21-48d0-9854-a2fcce904e8a
# ╟─9bd515d4-c7aa-4a3d-a4fb-28686290a134
# ╟─75dfaefd-abec-47e2-acc3-c0ff3a01048e
# ╠═f70aa9bd-4956-4e5a-8b3d-a1ad5b0b54d0
# ╠═7d738cd3-1edf-4348-aa16-5ade0f87ae6f
# ╠═63db8d8b-8229-47e0-a309-75813120aaab
# ╠═1ff66201-1d03-4c2e-95e5-d3387b0ddd57
# ╠═f999d732-14b3-4ac5-b803-3df7a96ef898
# ╠═35dccacb-8780-4f37-8dc5-29064a84d49c
# ╟─15245a9a-8a42-4ec1-ba2c-9bf95d6c85fa
# ╟─0d42b50d-993e-4eea-9025-2b7479bb3b0e
# ╠═83e265f4-0a62-4c53-9203-6a7768c6220e
# ╟─0ace8911-c28f-47cd-8117-581b0cc5e1d7
# ╟─23690382-3d30-46e3-b26a-a30875be78ec
# ╟─c92e091e-1111-4bda-8d76-543ad1dd8121
# ╟─379598be-2c42-4c29-8a24-91b74592da0f
# ╠═39cdac4b-87bb-441c-99e3-402b40bef7d3
# ╠═402a32b8-6b86-4701-a84b-eeb224e64a81
# ╟─1bd72281-4a7f-44f4-974d-632e9d0aaf28
# ╠═0c786da1-7f44-40af-b6d6-e0d6db2242b2
# ╠═5b80ca83-0719-437f-9e51-38f2bed02fb4
# ╠═cfcba5ad-e516-4223-860e-b1f18a6449ba
# ╟─34d4bd62-6be2-4089-8caa-1a8715bee433
# ╠═b7476391-30b9-4817-babf-7c9078531ee7
# ╠═c82ef027-740a-49b1-93d2-1554c411a896
# ╟─9873abd2-fc2f-4a35-8689-51d558debbf0
# ╠═0239e1da-caf5-4593-af1b-5d1e8d2f2b3e
# ╠═9e5a1672-d452-42f5-ba5a-a2fa0b1eaada
# ╠═2b634b6d-dc54-4ccb-8b41-7a881a018184
# ╠═73875562-ad97-4880-8694-0de3d088372a
# ╠═522e95d5-15a8-47ec-a79c-c4cc17cf86fd
# ╠═bfa4a8f9-cfcd-4e22-b2dc-751226f3a73c
# ╟─a8ccbc8e-24e6-4214-a179-4edf3cf26dad
# ╠═e8ee5cbb-4afc-4737-b006-90071f6138cd
# ╟─c8f71644-9371-443b-b976-1734cc7ae583
# ╠═57a37eb0-a547-428c-9f8c-e5f3e30f5260
# ╠═491da4e6-03e6-49d3-907d-43ddcffdfb42
# ╠═355ebed2-42e6-41bb-b1db-a72d1aaae56f
# ╠═a9d16b8a-a488-4a9e-bb6f-4ee05873e86b
# ╠═6cc82262-cc73-4483-b97e-664d5093d69c
# ╠═1098efd7-8477-4749-8ed5-b62044bebdc7
# ╠═c5f31bf6-5cca-4b0e-86f0-b44321fec874
# ╟─af6d4ef0-f59f-42be-af36-7cf447478e4c
# ╠═06498d79-f20a-4691-8c72-8f1f962e6a6f
# ╟─49aea0fc-1e0b-4e9c-8d7b-5d3727e95e1f
# ╠═871ca595-9d9e-4f06-9119-c2b07bfdeb04
# ╠═8eb901fe-aa93-4d72-bc9d-194e199db013
# ╠═5ae84f68-a0e6-47bb-9206-cd1384200581
# ╠═f15393ae-5c04-473e-ae28-ead5554896a3
# ╠═e33f9e0f-7c4c-465e-9793-ed790d3256b2
# ╟─3f9eb091-059c-44a5-9b50-ae3cabe24060
# ╠═4a8c7752-6de8-4ea7-aafe-702b17507185
# ╠═6e6b15b1-7685-4a20-9d94-dd703caa2fe9
# ╠═19b6abf5-6951-4492-8f17-f76df29f9289
# ╠═47e2e177-3701-471f-ae3c-38276ec69370
# ╠═457c1959-94fa-4267-8645-3ed1409cd0a0
# ╟─dfa2a03b-6925-4be0-aeac-076c4cf25969
# ╟─a7e75e49-5031-4cc4-b96e-6227277ec3ba
# ╟─d9617524-76c3-447d-9b94-0a690f83a7b9
# ╟─fb43471f-6aed-4fd6-a9e0-b165f6c77003
# ╟─9aded04b-e55f-4ebd-97c4-90c3adf62547
# ╟─30ec492c-8d21-43f6-bb09-32810494f21e
# ╠═856a78d9-7b4c-453b-b73b-c81eee014e52
# ╠═98a0d7c5-b1a8-4ebe-bb73-7ca88b475592
# ╟─bbbb358c-e645-4989-bed3-73d9217f7447
# ╠═6f08828b-4c4b-4f50-bd40-35805a37aae0
# ╟─cd5fe410-6264-4381-b19f-31d050bc3930
# ╠═0740dc70-a532-4818-b09d-b3b8d60fa6ba
# ╟─f26187fb-d4b2-4f0d-8a80-5d831e0de6c3
# ╠═19f4e0cc-0c93-42dc-8ee4-17f52d4e5e90
# ╠═4925c50b-12c0-4217-94de-bdcc72c01ccf
# ╟─4d9a4b36-6b3d-4836-8501-7f46cd7ab5cc
# ╟─75d956fc-bcf6-40fe-acd5-b5eef0fc7902
# ╟─c6ee857d-8019-4c4f-bb07-a370a88ea3cf
# ╠═6186798f-6711-4222-94bb-f53b2d0fad7d
# ╟─d27ef0d8-70b2-4897-9000-8fa70b1862fc
# ╟─6fc320b1-b60d-4f49-89ab-bf029ead6b55
# ╟─1da34733-fee3-42e1-b5e0-cac3f5f196c9
# ╟─32cd894b-ee5d-44b7-8983-82a4b72524a8
# ╠═c7deae02-3dad-4335-9449-a7e8f8bd5b4f
# ╠═b674af27-307b-4dbb-8a75-a54bde1f123d
# ╠═f7e0d09c-40bf-4936-987a-a3bcadae5487
# ╠═52c889e4-753c-447c-a9e1-862750b3643f
# ╠═aff80d55-df50-4d4b-aba4-e62f3c7ec10e
# ╠═a1e23c58-6d7b-4a69-8e33-411a7c051d37
# ╠═62f66995-bd02-4b6f-8eb8-6aeae5436713
# ╠═498c3c62-797c-4ea9-b8cb-61f2a91b269e
# ╠═57ea231c-86e7-44d3-9093-836c9820f88c
# ╟─59f3559b-aabe-42d7-9975-5fcc0b3de978
# ╠═acddad02-84ee-480f-a65f-716a4c34710c
# ╟─edabacdd-8d25-4d64-9d4a-ecf1263ac02e
# ╠═3c5edbc5-8fc7-4d09-98a9-85f2efb699a8
# ╟─b1f2aba6-5b6b-443c-84ab-21c4d2017a07
# ╟─67ad2564-fb20-4a71-a084-0145e8ed24bc
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
# ╟─d8d1fb74-0018-4685-a283-e768ae877fe4
# ╟─5f73f4e6-4eff-41b9-b68d-3baa5e77e924
# ╟─e94c4b92-ceec-412f-adfb-6e9a7344ca39
# ╟─c76f2ebe-ce41-47fd-b31a-2851aca53567
# ╟─25063860-6109-46e7-9dd5-a7fc0c12159e
# ╠═16094734-d7a1-4bb2-990b-8de1c367b134
