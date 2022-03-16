### A Pluto.jl notebook ###
# v0.18.1

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
	using Convex, ECOS
	using Plots
	using PlutoUI
	using LinearAlgebra
	using LightGraphs, SimpleWeightedGraphs
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

# ╔═╡ 0f9bfc53-8a1a-4e25-a82e-9bc4dc0a11fc
md"""This notebook aims at illustrating the method developed for computing marginal emissions, as well as the codebase built around it. 

Note this is work in progress. 
"""

# ╔═╡ 44275f74-7e7c-48d5-80a0-0f24609ef327
md"""
## Loading
"""

# ╔═╡ a32d6a56-8da8-44b0-b659-21030692630a
begin
	ECOS_OPT = () -> ECOS.Optimizer()
	OPT = ECOS_OPT
	δ = 1e-4
end;

# ╔═╡ 113e61a9-3b21-48d0-9854-a2fcce904e8a
xticks_hr = [0, 6, 12, 18, 24]

# ╔═╡ 935dabe1-467f-4c36-bdff-4cb6807b672f
theme(:default, label=nothing, 
		tickfont=(:Times, 8), guidefont=(:Times, 8), titlefont=(:Times,8), framestyle=:box)

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

# ╔═╡ cfcba5ad-e516-4223-860e-b1f18a6449ba
begin
	plt1 = plot(
		demand_data[:, rand(1:26, 5)], lw=2, palette=:Blues, title="Demand Data", xlabel="Hour", ylabel="Energy [?]", xticks=xticks_hr)
	
	plt1
end

# ╔═╡ 75dfaefd-abec-47e2-acc3-c0ff3a01048e
md"""
### Network
"""

# ╔═╡ e87dbd09-8696-43bd-84e0-af17517584dd
md"""
Instantiating a random network. 
"""

# ╔═╡ b63b9168-4213-4b95-bc2a-8888c032e9c5
sample(1:300, 4)

# ╔═╡ 6888f84d-daa9-4cfd-afc8-5aac00aeecab
begin
n = 6 # number of nodes
l = 6 # number of generators
T = 5 # number of timesteps
ns = 3 # number of storage nodes

net_dyn, _ = generate_network(n, l, T, ns)
end;

# ╔═╡ 176064fe-023b-49ca-ac06-f1a4e4be046c
begin
	n_demand = 2 # number of nodes with demand
	_, sd = size(demand_data)
	# constructing the demand vector
	idx_d  = rand(1:sd, n_demand)
	idx_ns = rand(1:n, n_demand)
	d_dyn =  [zeros(n) for _ in 1:T]
	for t in 1:T 
		d_dyn[t][idx_ns] = demand_data[t, idx_d]
	end
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
	type_idx = sample(1:length(EMISSIONS_FACTORS), l)
	keys_ = [kk for kk in keys(EMISSIONS_FACTORS)]
	emissions_rates = [EMISSIONS_FACTORS[keys_[ii]] for ii in type_idx]

end;

# ╔═╡ 91e525b2-a276-4641-bad0-7dd6be026790
mef_times = 1:T

# ╔═╡ cad822b2-4a01-4cc2-b557-fc1693d7a06e
s_rel = .1

# ╔═╡ 6f08828b-4c4b-4f50-bd40-35805a37aae0
begin
	results = zeros(n, T, length(mef_times))
	
	# Construct dynamic network
	η = .95

	# Construct and solve OPF problem
	opf_dyn = DynamicPowerManagementProblem(net_dyn, d_dyn)
	solve!(opf_dyn, OPT, verbose=false)

	println("Status of problem 1")
	@show opf_dyn.problem.status
end

# ╔═╡ d1d93fe7-79df-4253-86c2-31277ba792fc
begin
	# Compute MEFs
	mefs = compute_mefs(opf_dyn, net_dyn, d_dyn, emissions_rates)
	for ind_t in 1:length(mef_times)
		results[:, :, ind_t] .= mefs[ind_t]
	end
end

# ╔═╡ 15293269-580f-4251-be36-4be6ba8c5a46
md"""Influence of charging efficiency η on the total emissions of the system"""

# ╔═╡ 0740dc70-a532-4818-b09d-b3b8d60fa6ba
total_mefs = reshape(sum(results, dims=2), (n, length(mef_times))) 

# ╔═╡ 6870576e-4a46-44c4-978e-223fb4be96bc
@bind plot_id Slider(1:n)

# ╔═╡ 6186798f-6711-4222-94bb-f53b2d0fad7d
begin
	subplots_mef_storage = Dict()
	
	for (ind_plt, i) in enumerate(1:n)
		subplots_mef_storage[i] = plot(
			mef_times, total_mefs[i, :], 
			lw=4, alpha=0.8,
			# xlim=(1, 24),
			# ylim=(-1000, 4000), 
			# yticks = [0, 2000, 4000], 
			# xticks=xticks_hr,  
			xlabel=L"t_c", 
			ylabel="MEF", 
			legend=:topright
		)
	end
	plot(subplots_mef_storage[plot_id])

end

# ╔═╡ f7e0d09c-40bf-4936-987a-a3bcadae5487
begin
	plt_emissions_heatmap = Dict()
	heatmap_subplts = Dict()

	for node in 1:n
	# node = node_matrix# we focus on a single node
	
	plots = []
	lims = []

	crt_results = results[node, :, :]

	clim_ = 1.2*maximum(abs.(crt_results))
	subplt = heatmap(crt_results,
		c=:balance, #https://docs.juliaplots.org/latest/generated/colorschemes/
		clim=(-clim_, clim_), 
		colorbar=false,
		xlabel=L"t_c",
		ylabel=L"t_e"
	)
	
	# s_idx == 1 && plot!(ylabel=L"t_e")
	plot!(
		colorbar=true, 
		colorbar_tickfontsize=8,
		colorbar_ticks = [-500, 0, 500]
		# colorbar_title = "MEF
	)
	
	heatmap_subplts[node] = subplt
	end
		
end

# ╔═╡ 58373196-41ee-4d24-8aea-28b632d1c900
heatmap_subplts[plot_id]

# ╔═╡ edabacdd-8d25-4d64-9d4a-ecf1263ac02e
md"""
## Sensitivity analysis
"""

# ╔═╡ 3c5edbc5-8fc7-4d09-98a9-85f2efb699a8
node_sens = 2

# ╔═╡ 67ad2564-fb20-4a71-a084-0145e8ed24bc
cons_time = 4

# ╔═╡ 5365b74f-595f-4ade-a7af-e8dba53b84f7
md"""
Reference (as in computed) values
"""

# ╔═╡ a9b770e0-b137-40f7-b59a-35ad355b98bd
ref_mefs = results[node_sens, :, :]';

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

# ╔═╡ 91f7d63c-9e30-4fd4-ab39-9fbf58d101dc
begin
	println("Running sensitivity analysis")

	
	ref_val = deepcopy(d_dyn[cons_time][node_sens])
	if ref_val > 0 
			perturb_vals = [ref_val * (1+i*ε) for i in -npoints:npoints]
			x_axis_vals = [1+i*ε for i in -npoints:npoints]
			idx_ref = npoints+1
			idx_105 = findall(x_axis_vals.==1.05)[1]
			idx_95 = findall(x_axis_vals.==.95)[1]
			idx_110 = findall(x_axis_vals .== 1.1)[1]
			idx_90 = findall(x_axis_vals .== .9)[1]
			plot_sens_flag = true
	else # Demand is zero
			perturb_vals = [i*ε for i in 0:npoints]
			x_axis_vals = perturb_vals
			idx_ref = 1
			plot_sens_flag = false
	end
	L = length(perturb_vals)
	E_sensitivity = zeros(L, length(emissions_rates), T);
	s_sensitivity = zeros(L, ns, T)
	g_sensitivity = zeros(L, l, T);
	mefs_sensitivity = zeros(L, T)
		
	for k in 1:L
		d_crt = deepcopy(d_dyn)
		d_crt[cons_time][node_sens] = perturb_vals[k]
		opf_ = DynamicPowerManagementProblem(net_dyn, d_crt)
		solve!(opf_, OPT, verbose=false)
		if opf_.problem.status != Convex.MOI.OPTIMAL
			@show opf_.problem.status
		end

		println("Computing MEFs")
		mefs_ = compute_mefs(opf_, net_dyn, d_crt, emissions_rates)
		
		for t in 1:T
			@show evaluate(opf_.s[t])
			@show ns
			@show size(s_sensitivity)
			s_sensitivity[k, :, t] .= evaluate(opf_.s[t])
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

# ╔═╡ 77943ac8-36fe-4a13-a36d-db957780d869
begin #E_ref is the total emissions at a given time
	E_ref = zeros(T)
	
	for t in 1:T
		E_ref[t] = evaluate(opf_dyn.g[t])' * emissions_rates
	end

end

# ╔═╡ 4fd2833c-6c23-4009-8734-980d3dd08c91
md"""
What is the value of emissions when there is no perturbation? 
"""

# ╔═╡ 6fcd6e19-58c3-462d-964f-8cd3127b47a4
sum(E_sensitivity[npoints+1, :, :], dims=1)

# ╔═╡ 2973af52-0bd0-4ba8-855d-297427627e22
E_ref[:]

# ╔═╡ b85b85d0-e1bc-4fc9-81cf-3792b55e3684
e_time = 1

# ╔═╡ 30511293-8ba5-486e-956b-e9f2a1ed0505
begin
	γ = 1e-4
	Δ = .05
	ylims = (1-Δ, 1+Δ)
	plt_s = plot(
		x_axis_vals, 
		[s_sensitivity[:, k, e_time]/(s_sensitivity[idx_ref, k, e_time]+γ) for k in 1:ns], ylim=ylims
	)
	title!("Storage at time $e_time")
	xlabel!("Change in demand at node $node_sens at time $cons_time")
	ylabel!("Change in storage at all nodes at time $e_time")
	
	plt_E = plot(
		x_axis_vals, 
		[E_sensitivity[:, k, e_time]./(E_sensitivity[idx_ref, k, e_time]+γ) for k in 1:length(emissions_rates)], ylim=ylims
		)
	title!("Emissions at time $e_time")
	xlabel!("Change in demand at node $node_sens at time $cons_time")
	ylabel!("Change in emissions at all generators at time $e_time")
	
	plt_g = plot(
		x_axis_vals, 
		[g_sensitivity[:, k, e_time]./(g_sensitivity[idx_ref, k, e_time]+γ) for k in 1:length(emissions_rates)], ylim=ylims
		)
	title!("Generators at time $e_time")
	xlabel!("Change in demand at node $node_sens at time $cons_time")
	ylabel!("Change in generation at all generators at time $e_time")

	norm_E = sum(E_sensitivity[:, :, e_time], dims=2)./sum(E_sensitivity[idx_ref, :, e_time])
	plt_E_tot = plot(
		x_axis_vals, norm_E
		, ylim=(0.90, 1.1)
		)
	# xlabel!("Change in demand at node $node_sens at time $cons_time")
	# ylabel!("Change in total emissions")
	xlabel!(L"\Delta d/d")
	ylabel!(L"\Delta E/E")
	
	#adding the theoretical curve for the sensitivity
	E_th = (
		sum(E_sensitivity[idx_ref, :, e_time]) .+ (perturb_vals.-ref_val) .* mefs_sensitivity[idx_ref, e_time]
	)./sum(E_sensitivity[idx_ref, :, e_time])

	if plot_sens_flag
	E_th_105 = (
		sum(E_sensitivity[idx_105, :, e_time]) .+ (perturb_vals.-ref_val * x_axis_vals[idx_105]) .* mefs_sensitivity[idx_105, e_time]
	)./sum(E_sensitivity[idx_ref, :, e_time])

E_th_95 = (
		sum(E_sensitivity[idx_95, :, e_time]) .+ (perturb_vals.-ref_val * x_axis_vals[idx_95]) .* mefs_sensitivity[idx_95, e_time]
)./sum(E_sensitivity[idx_ref, :, e_time])
	E_th_110 = (
		sum(E_sensitivity[idx_110, :, e_time]) .+ (perturb_vals.-ref_val * x_axis_vals[idx_110]) .* mefs_sensitivity[idx_110, e_time]
	)./sum(E_sensitivity[idx_ref, :, e_time])

	E_th_90 = (
		sum(E_sensitivity[idx_90, :, e_time]) .+ (perturb_vals.-ref_val * x_axis_vals[idx_90]) .* mefs_sensitivity[idx_90, e_time]
	)./sum(E_sensitivity[idx_ref, :, e_time])
	end
	
	plot!(x_axis_vals, E_th, ls=:dash, c=:orange)
	
	if plot_sens_flag
	plot!(x_axis_vals, E_th_105, ls=:dash, c=:green)
	plot!(x_axis_vals, E_th_95, ls=:dash, c=:firebrick)
	plot!(x_axis_vals, E_th_110, ls=:dash, c=:violetred)
	end

	scatter!([x_axis_vals[idx_ref]], [norm_E[idx_ref]], c=:orange)

	if plot_sens_flag
	scatter!([x_axis_vals[idx_105]], [norm_E[idx_105]], c=:green)
	scatter!([x_axis_vals[idx_95]], [norm_E[idx_95]], c=:firebrick)
	scatter!([x_axis_vals[idx_110]], [norm_E[idx_110]], c=:violetred)
	end
	title!("Total emissions at time $e_time")
	
	@show ref_mefs[e_time, cons_time]
	@show e_time
	@show cons_time
	@show ref_val
	
	plot([plt_s, plt_E, plt_g, plt_E_tot]..., size = (650, 650), lw = 3)
	
end

# ╔═╡ Cell order:
# ╟─c39005df-61e0-4c08-8321-49cc5fe71ef3
# ╟─0f9bfc53-8a1a-4e25-a82e-9bc4dc0a11fc
# ╟─44275f74-7e7c-48d5-80a0-0f24609ef327
# ╠═db59921e-e998-11eb-0307-e396d43191b5
# ╠═0aac9a3f-a477-4095-9be1-f4babe1e2803
# ╠═a32d6a56-8da8-44b0-b659-21030692630a
# ╠═113e61a9-3b21-48d0-9854-a2fcce904e8a
# ╠═935dabe1-467f-4c36-bdff-4cb6807b672f
# ╟─9bd515d4-c7aa-4a3d-a4fb-28686290a134
# ╟─1bd72281-4a7f-44f4-974d-632e9d0aaf28
# ╠═0c786da1-7f44-40af-b6d6-e0d6db2242b2
# ╠═cfcba5ad-e516-4223-860e-b1f18a6449ba
# ╟─75dfaefd-abec-47e2-acc3-c0ff3a01048e
# ╟─e87dbd09-8696-43bd-84e0-af17517584dd
# ╠═b63b9168-4213-4b95-bc2a-8888c032e9c5
# ╠═6888f84d-daa9-4cfd-afc8-5aac00aeecab
# ╠═176064fe-023b-49ca-ac06-f1a4e4be046c
# ╟─a8ccbc8e-24e6-4214-a179-4edf3cf26dad
# ╟─496135ec-f720-4d43-8239-d75cc7616f58
# ╠═806819d5-7b40-4ca4-aa1e-f1cf0a9a7f3f
# ╠═91e525b2-a276-4641-bad0-7dd6be026790
# ╠═cad822b2-4a01-4cc2-b557-fc1693d7a06e
# ╠═6f08828b-4c4b-4f50-bd40-35805a37aae0
# ╠═d1d93fe7-79df-4253-86c2-31277ba792fc
# ╟─15293269-580f-4251-be36-4be6ba8c5a46
# ╠═0740dc70-a532-4818-b09d-b3b8d60fa6ba
# ╠═6870576e-4a46-44c4-978e-223fb4be96bc
# ╟─6186798f-6711-4222-94bb-f53b2d0fad7d
# ╟─f7e0d09c-40bf-4936-987a-a3bcadae5487
# ╠═58373196-41ee-4d24-8aea-28b632d1c900
# ╟─edabacdd-8d25-4d64-9d4a-ecf1263ac02e
# ╠═3c5edbc5-8fc7-4d09-98a9-85f2efb699a8
# ╠═67ad2564-fb20-4a71-a084-0145e8ed24bc
# ╟─5365b74f-595f-4ade-a7af-e8dba53b84f7
# ╠═a9b770e0-b137-40f7-b59a-35ad355b98bd
# ╠═956e963a-97af-495e-9475-181322ac2e0c
# ╠═4aed3df5-441b-445b-9277-a38690eb8603
# ╟─c9b41436-e0a0-4e57-908f-b45e42122e63
# ╟─91f7d63c-9e30-4fd4-ab39-9fbf58d101dc
# ╠═77943ac8-36fe-4a13-a36d-db957780d869
# ╟─4fd2833c-6c23-4009-8734-980d3dd08c91
# ╠═6fcd6e19-58c3-462d-964f-8cd3127b47a4
# ╠═2973af52-0bd0-4ba8-855d-297427627e22
# ╠═b85b85d0-e1bc-4fc9-81cf-3792b55e3684
# ╠═30511293-8ba5-486e-956b-e9f2a1ed0505
