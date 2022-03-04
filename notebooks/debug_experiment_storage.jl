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

# ╔═╡ fc7f1535-18ae-4f76-ae99-06839728360f
using LightGraphs, SimpleWeightedGraphs, GraphPlot

# ╔═╡ efc01f26-6b99-4a70-b34b-b022fb0d6f5b
md"""
## TODO
- Conduct the battery placement analysis? i.e. depending on where you place the battery you have the emergence of squares depending on where it is? (i.e. close nodes should be more sensitive to the presence of the battery?)

- *make sure* that the orders in which you query the `mef` matrix is correct (cons time vs t_display etc.)

- Introduce larger discrepancies in the line capacities, right now it's a uniform scaling

- Introduce renewable generation and locational storage and see how the system evolves with that...

- look at curtailment of renewables

- look at congestion

- revisit how the total demand is sized

- they say the problem is solved to optimality, but we see non empty batteries at the last timestep...?

- why does storage location update every time? 

- if you add some transmission cost, how do the mefs scale with distance? 

- change the way in which line capacity limits are set, because it seems to be too strong a parameter for now... and for small networks you can probably change them independently

- add transmission costs?
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

# ╔═╡ 538bf950-6a1f-431a-b4ea-aefcf833bded
DIFFERENT_LINEAR_COSTS = false

# ╔═╡ db4f89b9-0ffc-41d9-a46e-1eb0b6fb6f83
# Define the properties of the network
begin
	n_ = 3
	l_ = 3
end;

# ╔═╡ 365ed8da-0f55-4d98-a570-6bc68f328cc2
begin
node_storage = 2 # node that has storage. if 0, all nodes do.
node_renewable = 1 # node that has renewable generation. if 0 all nodes do
node_demand = 0 # node that has the demand. if 0 all nodes do
end;

# ╔═╡ c058136b-e536-41be-8c7a-28a8c51f3b22
md"""
Put the storage in a single location: $node_storage. 

If 0, then all the nodes have storage. 
"""

# ╔═╡ 05405ddb-7046-400b-9ae6-4899f37de6a3
@bind percent_renewable Slider(0:0.1:1)

# ╔═╡ a2932123-66a5-4d29-a794-3026d8984aff
md"""Percent of total generation capacity provided by renewable: $(100*percent_renewable) %"""

# ╔═╡ fb3e2e73-8b0b-43ba-a102-39ad9599941f
begin
	if LOAD_RANDOM_GRAPH
		A, B, cq, cl, _, gmax, pmax, _, _ = generate_random_data(n_, l_, 1)
		m_, _ = size(B)
		β = rand(Uniform(0, 1), m_)
		F = make_pfdf_matrix(A, β)
		
		if ~ QUAD_COSTS
			cq = [zeros(l_)]
		end
		
		if node_renewable > 0
			cq[1][node_renewable] = 0
			cl[1][node_renewable] = 0
			gmax[1][node_renewable]  = percent_renewable/(1-percent_renewable) * (
				sum([gmax[1][k] for k in 1:l_]) - gmax[1][node_renewable]
			)
		end
		
		net = PowerNetwork(cq[1], cl[1], pmax[1], gmax[1], A, B, F)
	else
		net, _, _ = load_synthetic_network("case14.m");
	end
	
	emissions_rates = rand(Exponential(1), l_)
	
	if node_renewable > 0
		emissions_rates[node_renewable] = 0
	end
	
end;

# ╔═╡ 75be5683-85a8-45c3-925c-7a53dbf9ae7d
net.gmax

# ╔═╡ ca2263db-d63a-4f26-9643-0aefa77b4d7c
emissions_rates

# ╔═╡ 0b806e4a-73a9-48f1-993c-875d0c9a01c3
pmax_ref = net.pmax;

# ╔═╡ f8072762-643c-485d-86d7-caf5a54c7d1e
n, m, l = get_problem_dims(net);

# ╔═╡ 058f7c12-d237-4130-aeb2-99279953d3f8
begin
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

# ╔═╡ 4a0c8a47-4aea-41ef-8a45-e2ea5789721f
spens = 0:0.1:3;

# ╔═╡ b77196c7-987d-4175-9508-88c11cedbc3c
@bind spen Slider(spens)

# ╔═╡ 2085073b-4770-4815-9e38-c36b850ab8d4
md"""
spen = $spen
"""

# ╔═╡ 01ec04b7-e317-43b9-86b8-8a7e9ecbcdea
η = 1

# ╔═╡ c9a5eeec-389b-4dbf-9c57-3ddd7ebff264
begin #alphas are currently not used
α_max = 30
step = .1
αs = 1:step:α_max;
end;

# ╔═╡ c35434eb-c850-4650-870e-b4b719ff1a9b
@bind αi Slider(1:1:length(αs))

# ╔═╡ 76e07e5e-8022-40e6-84aa-84bb050c7545
begin
	p1 = 500
	p2 = 500
	p3 = 1
end;

# ╔═╡ 6198b321-6573-4100-8166-6048c5fa2980
α = αs[αi];

# ╔═╡ 5c07c744-639b-4f56-af9a-a4b2366864bb
md"""α = $α"""

# ╔═╡ 377a4e9e-8d8c-4a76-9965-8df0a21dcf9c
# TODO : revisit how the demand is sized
# TODO: revisit how the power line constraints are set, too

begin
	pmax_mat = zeros(m, length(αs));
	for i in 1:length(αs)
		pmax_mat[:, i] = αs[i]*ones(m);
	end
	
	d_dyn = [rand(Bernoulli(0.8), n) .* rand(Exponential(2), n) for _ in 1:T];
	if node_demand > 0
		d_dyn = [zeros(n) for _ in 1:T];
		for t in 1:T
			d_dyn[t][node_demand] = rand(Exponential(2));
		end
	end
	for i in 1:T
		if node_renewable > 0
			d_dyn[i] = d_dyn[i]/sum(d_dyn[i]) * (
					sum(net.gmax) - net.gmax[node_renewable]) * rand(Uniform(.5, 1.)
					);
		else
			d_dyn[i] = d_dyn[i]/sum(d_dyn[i]) * (
					sum(net.gmax)) * rand(Uniform(.5, 1.)
					);
		end
	end

	
end;

# ╔═╡ cadc304e-24d3-4aa7-8858-24255effaa13
d_dyn

# ╔═╡ 6726672d-ab82-4364-bcea-bec33034bdac
begin
	
	mefs = zeros(
			n, T, length(mef_times)
		)
	println("Solving problem with parameters:")
	println("s = $spen, η = $η, α=$α, renewable_pc = $(100*percent_renewable)%")
	
	# Construct dynamic network
	C = sum(d_dyn) .* (spen + δ)/T
	if node_storage >= 1
		nodes_on = zeros(n)
		nodes_on[node_storage] = 1
		C = C.* nodes_on
	end
	P = 100*C 
	# net.pmax = pmax_mat[:, αi];
	net.pmax = [p1, p2, p3];

	net_dyn = make_dynamic(net, T, P, C, η);
	if DIFFERENT_LINEAR_COSTS
		net_dyn.fl = [rand(Exponential(2), l) .* net_dyn.fl[1] for _ in 1:T]
	end

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

# ╔═╡ 79eb4f29-1457-41b4-96c3-9bfbe7a54208
begin
	
	subplts_nodes = []
	lim = max(abs(minimum(mefs)), abs(maximum(mefs)));
	clims = (-lim, lim)
	for node_ in 1:n
		crt_map = mefs[node_, :, :]
		subplt = heatmap(crt_map, 
		c=:balance, colorbar=false,
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

# ╔═╡ c1f12b9d-60fe-4c10-b74b-a63dfd965577
md"""Cons time = $cons_time"""

# ╔═╡ c62e4fe2-2971-4479-af09-11b5467b2fee
@bind t_display Slider(1:1:T)

# ╔═╡ 5bd2052c-8f4b-4969-8864-8906ece79df5
md"""
Emissions time: $t_display
"""

# ╔═╡ 971d508d-c862-44ac-a0a4-a2fd9aa9156b
@bind node Slider(1:1:n)

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
	
	plots = plot(subplts..., 
		layout=lay
		# size=(650, 200), 
		# bottom_margin=8Plots.pt
	)

	
	plots
	
	
	
	
end

# ╔═╡ 2b8ccbdc-bbbe-4a0a-a71b-52fad873bc70
md""" node = $node"""

# ╔═╡ c6d3c265-d3e9-4f4e-b5e7-23a3762beb49
mefs[node, :, cons_time]

# ╔═╡ c9e229f9-8ce6-4cd4-a0d7-4f53f2ced8ec
begin
	println("Running sensitivity analysis")
	# size of the matrices are
	# 2npoints+1: number of different values of demand for which we solve the problem
	# n: number of nodes in the graph
	# l: number of generators (= length(emissions_rates))
	# T: the time horizon

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
	
	for k in 1:L
		d_crt = deepcopy(d_dyn)
		d_crt[cons_time][node] = perturb_vals[k]

		opf_ = DynamicPowerManagementProblem(net_dyn, d_crt)
		solve!(opf_, OPT, verbose=false)
		if opf_.problem.status != Convex.MOI.OPTIMAL
			@show opf_.problem.status
		end

		for t in 1:T
			s_sensitivity[k, :, t] = evaluate(opf_.s[t])
			g_sensitivity[k, :, t] = evaluate(opf_.g[t])
			E_sensitivity[k, :, t] = evaluate(opf_.g[t]).*emissions_rates
		end

	end
end

# ╔═╡ c418c302-2bd5-40c5-aca6-e0302372903a
md"""node = $node"""

# ╔═╡ fa2d971a-38d6-498e-b1fd-71e8ef420e7c
net.pmax

# ╔═╡ 418c2024-64c8-4dd3-ac1f-0d71ff99edd8
evaluate(opf_dyn.p[1])

# ╔═╡ 1182f86d-0d85-4174-8cfe-144bf2e38a0f
net.gmax

# ╔═╡ 8c57c424-d4bc-4ff4-95a2-b301561d8aeb
mefs

# ╔═╡ 2aa296d4-7f47-4585-8a01-c70e9dddd2ac
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

# ╔═╡ 4a43aedb-b17e-459b-aadb-86528d1bb7c4
md"""Visualizing total mefs"""

# ╔═╡ 96a33baa-9cb2-446c-a850-076ccf503735
let
bar_plt = bar([mefs[k, t_display, cons_time] for k in 1:n])
xlabel!("node")
ylabel!("nodal mefs at time $t_display")
	
bar_plt2 = bar([sum(mefs[k, :, cons_time]) for k in 1:n])
xlabel!("node")
ylabel!("total nodal mefs")
	
plot(bar_plt, bar_plt2)

end

# ╔═╡ 2061e505-62b3-4e50-b91d-2852fcfc24ad
md"""
Is there renewable curtailment? As in -- can they not use all the renewable generation because of congestion? 
"""

# ╔═╡ 438f18ea-5a4c-4916-91c2-9be5f2e8fcc1
let
	if node_renewable > 0
		plt1 = plot(g_vals[node_renewable, :], ylim=(0, 1.1))
		xlabel!("Time")
		ylabel!("Generation of renewable")
		title!("Relative")

		plt2 = plot(g_vals[node_renewable, :] * net.gmax[node_renewable])
		xlabel!("Time")
		ylabel!("Generation of renewable")
		title!("Absolute")
		
		plot(plt1, plt2)
	end
end

# ╔═╡ 9d98607d-d6a8-48ba-996b-aab22be9b3a9
net.gmax

# ╔═╡ d2e380ba-6b1a-4602-ab87-b4926eaa8125
md"""
## Evolution of mefs as a function of storage penetration
"""

# ╔═╡ ca0e5b0a-f6c4-492d-9f44-e903e256e4d3
md"""Cons time = $cons_time"""

# ╔═╡ 20f0a35f-3aec-4b72-8ee2-eaca1e98dac0
md"""Em time = $t_display"""

# ╔═╡ 75197a33-b583-4a5e-bc2f-20e091ad0200
md"""Node = $node"""

# ╔═╡ 76a36e77-bbcd-42ad-a3f7-5a1b45d1a767
RUN_CELL_STORAGE = true

# ╔═╡ 4587a3ba-f478-4ae9-ba2a-15ea8755fa75
mefs_spens = zeros(n, length(mef_times), T, length(spens));

# ╔═╡ e034f565-9926-491b-b7b0-a09bf873f0b2
let
	if RUN_CELL_STORAGE
		
		println("Running cell Storage")

		for k in 1:length(spens)
			C = sum(d_dyn) .* (spens[k] + δ)/T
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
				mefs_spens[:, :, ind_t, k] .= mefs_[ind_t];
			end
		end
	end
end

# ╔═╡ 6bfe841c-2f05-4e30-9edb-2181b64c92c1
begin
	
	plot()
	for t in mef_times
		plot!(spens, [mefs_spens[node, cons_time, t, k] for k in 1:length(spens)], lw=3, ls=:dash, markershape=:circle)
	end
	plot!()
	xlabel!("Storage pen")
	ylabel!("MEFs")
	title!("MEFs as a function of storage pen for constime = $cons_time")
end

# ╔═╡ a873fab5-d047-4670-8136-b7dbe9363a19
md"""Above plot seems to show that at a given consumption time, the mefs converge"""

# ╔═╡ f889110b-add4-4a71-b951-6cfd8879c030
md"""
How about total mefs? 
"""

# ╔═╡ 22330790-1627-4555-b4cd-0f554c000385
let
	plot()
	for ct in 1:T # looping over consumption times
		plot!(
			spens, 
			[sum(mefs_spens[node, ct, :, k]) for k in 1:length(spens)], 
			lw=3, ls=:dash, markershape=:circle
			)
	end
	plot!()
	xlabel!("Storage Pen")
	ylabel!("Total mef")
	title!("Total MEF as a function of storage pen for node $node")
end

# ╔═╡ 3114e5b2-c969-4bb6-99e8-432690a7152b
mefs_spens[node, :, cons_time, end]


# ╔═╡ c5a9c260-f51c-439e-8d40-806444c6e243
mefs[node, :, cons_time]

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
# ╠═efc01f26-6b99-4a70-b34b-b022fb0d6f5b
# ╠═fc7f1535-18ae-4f76-ae99-06839728360f
# ╠═2498bfac-3108-11ec-2b8b-7fb26f96afbb
# ╠═6a260a4f-9f96-464b-b101-1127e6ec48fe
# ╠═f94d2b5b-779a-4de0-9753-c077bc925fa1
# ╠═1be301f4-31fe-44e9-895a-49bb1eec512f
# ╠═6b706efa-3343-4b9a-bd2f-1bf263707836
# ╠═7d4a1f79-2a2f-4179-bdb6-ea394b7ca5fb
# ╟─5aef540d-719e-4212-9e2d-40c75cb685a7
# ╠═06070be9-fae4-41b1-bdd3-066f7e785439
# ╠═538bf950-6a1f-431a-b4ea-aefcf833bded
# ╠═db4f89b9-0ffc-41d9-a46e-1eb0b6fb6f83
# ╟─c058136b-e536-41be-8c7a-28a8c51f3b22
# ╠═365ed8da-0f55-4d98-a570-6bc68f328cc2
# ╟─a2932123-66a5-4d29-a794-3026d8984aff
# ╟─05405ddb-7046-400b-9ae6-4899f37de6a3
# ╠═75be5683-85a8-45c3-925c-7a53dbf9ae7d
# ╠═ca2263db-d63a-4f26-9643-0aefa77b4d7c
# ╠═cadc304e-24d3-4aa7-8858-24255effaa13
# ╠═fb3e2e73-8b0b-43ba-a102-39ad9599941f
# ╠═058f7c12-d237-4130-aeb2-99279953d3f8
# ╠═0b806e4a-73a9-48f1-993c-875d0c9a01c3
# ╠═f8072762-643c-485d-86d7-caf5a54c7d1e
# ╟─65b961b0-cc87-450e-bedb-43c5f38aafc5
# ╟─73dac19a-3663-4f25-ac96-24e3a9997d3e
# ╠═4a0c8a47-4aea-41ef-8a45-e2ea5789721f
# ╟─2085073b-4770-4815-9e38-c36b850ab8d4
# ╠═b77196c7-987d-4175-9508-88c11cedbc3c
# ╠═01ec04b7-e317-43b9-86b8-8a7e9ecbcdea
# ╠═c9a5eeec-389b-4dbf-9c57-3ddd7ebff264
# ╠═5c07c744-639b-4f56-af9a-a4b2366864bb
# ╠═c35434eb-c850-4650-870e-b4b719ff1a9b
# ╠═76e07e5e-8022-40e6-84aa-84bb050c7545
# ╠═6198b321-6573-4100-8166-6048c5fa2980
# ╠═377a4e9e-8d8c-4a76-9965-8df0a21dcf9c
# ╠═6726672d-ab82-4364-bcea-bec33034bdac
# ╠═488db7c4-a048-47bd-b2bf-430d7f5664ae
# ╟─f59c420b-a1fe-4847-8115-a33166be54aa
# ╟─2b8ccbdc-bbbe-4a0a-a71b-52fad873bc70
# ╠═c6d3c265-d3e9-4f4e-b5e7-23a3762beb49
# ╠═79eb4f29-1457-41b4-96c3-9bfbe7a54208
# ╟─bdd10c3a-5b1f-4f63-860c-268b9c4f035d
# ╟─c1f12b9d-60fe-4c10-b74b-a63dfd965577
# ╟─124bc7c3-bb3b-4b36-bb7e-4f3f52b7de58
# ╠═c9e229f9-8ce6-4cd4-a0d7-4f53f2ced8ec
# ╟─5bd2052c-8f4b-4969-8864-8906ece79df5
# ╟─c62e4fe2-2971-4479-af09-11b5467b2fee
# ╟─c418c302-2bd5-40c5-aca6-e0302372903a
# ╟─971d508d-c862-44ac-a0a4-a2fd9aa9156b
# ╠═fa2d971a-38d6-498e-b1fd-71e8ef420e7c
# ╠═418c2024-64c8-4dd3-ac1f-0d71ff99edd8
# ╠═1182f86d-0d85-4174-8cfe-144bf2e38a0f
# ╠═8c57c424-d4bc-4ff4-95a2-b301561d8aeb
# ╠═2aa296d4-7f47-4585-8a01-c70e9dddd2ac
# ╟─4a43aedb-b17e-459b-aadb-86528d1bb7c4
# ╟─96a33baa-9cb2-446c-a850-076ccf503735
# ╟─2061e505-62b3-4e50-b91d-2852fcfc24ad
# ╠═438f18ea-5a4c-4916-91c2-9be5f2e8fcc1
# ╠═9d98607d-d6a8-48ba-996b-aab22be9b3a9
# ╟─d2e380ba-6b1a-4602-ab87-b4926eaa8125
# ╟─ca0e5b0a-f6c4-492d-9f44-e903e256e4d3
# ╟─20f0a35f-3aec-4b72-8ee2-eaca1e98dac0
# ╟─75197a33-b583-4a5e-bc2f-20e091ad0200
# ╟─76a36e77-bbcd-42ad-a3f7-5a1b45d1a767
# ╠═4587a3ba-f478-4ae9-ba2a-15ea8755fa75
# ╟─e034f565-9926-491b-b7b0-a09bf873f0b2
# ╟─6bfe841c-2f05-4e30-9edb-2181b64c92c1
# ╟─a873fab5-d047-4670-8136-b7dbe9363a19
# ╟─f889110b-add4-4a71-b951-6cfd8879c030
# ╟─22330790-1627-4555-b4cd-0f554c000385
# ╠═3114e5b2-c969-4bb6-99e8-432690a7152b
# ╠═c5a9c260-f51c-439e-8d40-806444c6e243
# ╟─9f8ef7ba-8e3d-4ac8-bb7f-44be3dee822d
# ╟─2892cc44-a2a6-4e36-b9fc-a2d5e4e564bc
# ╟─0ed86b74-66e6-428e-84cc-5afe96e49060
# ╟─1132387c-80b1-416d-b167-909be1e3e3ab