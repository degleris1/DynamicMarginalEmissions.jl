### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ db59921e-e998-11eb-0307-e396d43191b5
begin
	import Pkg
	Pkg.activate(joinpath(@__DIR__, "../scratch"));
	using Random
	using Convex, ECOS
	using Plots
	using PlutoUI
	using LinearAlgebra
	using LightGraphs, SimpleWeightedGraphs
end;

# ╔═╡ 0aac9a3f-a477-4095-9be1-f4babe1e2803
begin
	using Revise
	using DynamicMarginalEmissions
end

# ╔═╡ 0107b2ad-5e74-49ce-9231-26caed98c5c3
using LaTeXStrings

# ╔═╡ 0f7ed93f-72bf-41de-ba6d-f46926d9fc46
using Distributions

# ╔═╡ c39005df-61e0-4c08-8321-49cc5fe71ef3
md"""
## Description
"""

# ╔═╡ 0f9bfc53-8a1a-4e25-a82e-9bc4dc0a11fc
md"""This notebook aims at illustrating the method developed for computing marginal emissions, as well as the codebase built around it. 
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

# ╔═╡ 935dabe1-467f-4c36-bdff-4cb6807b672f
theme(:default, label=nothing, 
		tickfont=(:Times, 8), guidefont=(:Times, 8), titlefont=(:Times,8), framestyle=:box)

# ╔═╡ 9bd515d4-c7aa-4a3d-a4fb-28686290a134
md"""
## Load data
"""

# ╔═╡ 1bd72281-4a7f-44f4-974d-632e9d0aaf28
md"""
### Demand time series
"""

# ╔═╡ 0c786da1-7f44-40af-b6d6-e0d6db2242b2
demand_data = .5 * load_demand_data("2021_07_01", normalize_rows=true);

# ╔═╡ cfcba5ad-e516-4223-860e-b1f18a6449ba
begin
	plt1 = plot(
		demand_data[:, rand(1:26, 5)], lw=2, palette=:Blues, title="Demand Data", xlabel="Hour", ylabel="Demand")
	
	plt1
end

# ╔═╡ 75dfaefd-abec-47e2-acc3-c0ff3a01048e
md"""
### Network data
"""

# ╔═╡ c31093bf-4bd6-4642-81cf-377895799de7
begin
	θ, d, net, β = load_synthetic_network("case14.m")
	θ.pmax = θ.pmax/4
end;

# ╔═╡ d225147d-658b-4f3e-87ca-be7c33aefe33
#dimensions of the problem
begin
n, m = size(θ.A)
k = size(θ.B, 2)
end;

# ╔═╡ e2bb5fc1-32e1-4dd5-95ec-f39c4ee7fa14
# user inputs
begin
	T = 23 # time horizon
	ns = 2 # number of batteries in the network
end;

# ╔═╡ d6c7b48a-bb28-413d-9ebc-f3675351c302
begin
	dd = demand_data[1:T, rand(1:size(demand_data,2), n)];
	demand_schedule = [dd[k, :] for k in 1:size(dd, 1)] * .8;
end;

# ╔═╡ 40095334-a083-4abc-945d-d1efa1d285e8
begin # randomly assign ns storage nodes to the network
	S = zeros(n, ns);
	idx_s = rand(1:n, ns)
	for j in 1:ns
		S[idx_s[j], j] = 1.
	end
end

# ╔═╡ fcd7896a-cc91-4b97-9210-f87362965cbd
# battery parameters
begin
	# penetration as % of total generation capacity
	battery_penetration = .1
	C = battery_penetration * sum(θ.gmax)/ns * ones(ns) # Battery capacities
	P = 0.25 * C # Battery ramping rate
end;

# ╔═╡ f37f2ec5-b3a4-412c-9c10-da2cde64ec38
# make network a dynamic network
θ_dyn = make_dynamic(θ, T, S, P, C);

# ╔═╡ 806819d5-7b40-4ca4-aa1e-f1cf0a9a7f3f
begin
	gen_types = ["COL", "SUN", "NG", "WND", "NG"]
	emissions_rates = [EMISSIONS_FACTORS[gen_type] for gen_type in gen_types];
end;

# ╔═╡ 4fbf7a22-040f-4558-8a41-b581d2b98e94
begin

params = Dict(
	:time_horizon => T,
	:demand_schedule => demand_schedule,
	
	:node_line_incidence_matrix => θ_dyn.A, 
	:node_generator_incidence_matrix => θ_dyn.B,
	:node_battery_incidence_matrix => S,

	:generator_costs_quadratic => θ_dyn.fq,
	:generator_costs_linear => θ_dyn.fl,
	:generator_capacities => θ_dyn.gmax,
	:generator_emissions_rates => emissions_rates,
	:generator_ramping_rates => nothing,

	:line_susceptances => β,
	:line_capacities => θ_dyn.pmax,

	:battery_capacities => C,
	:battery_max_powers => P,
	)
end;

# ╔═╡ 7d4b9df7-2181-4af7-a947-9a81297fe907
lmes, tot_lmes,  g, meta = dispatch_and_compute_emissions_rates(; params...);

# ╔═╡ a8d82fd9-b7bd-4f8a-9763-7ed684766e58
md"""# Some Plots"""

# ╔═╡ 25659734-b76d-4498-826a-18b5d751b37c
begin
	plot(transpose(hcat(demand_schedule...)))
	xlabel!("Time step")
	ylabel!("Nodal demand")
end

# ╔═╡ fb910c5e-5b81-4103-9829-efab4f920f54
begin
	plot(transpose(tot_lmes))
	xlabel!("Time step")
	ylabel!("Total nodal LME")
end

# ╔═╡ edabacdd-8d25-4d64-9d4a-ecf1263ac02e
md"""
## Sensitivity analysis
"""

# ╔═╡ 380103e3-d893-4fc1-ae38-110a21cd863d
# user input
begin
	node_sens = 4 # node on which we run the sensitivity analysis
	cons_time = 20 # consumption time
	e_time = 20 # emissions time
end;

# ╔═╡ 4aed3df5-441b-445b-9277-a38690eb8603
begin
npoints = 10
ε = 1e-2
end;

# ╔═╡ 91f7d63c-9e30-4fd4-ab39-9fbf58d101dc
begin # running sensitivity analysis
	
	d_dyn = deepcopy(demand_schedule)
	# local copy of the demand schedule
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

	#Initiate the local variables
	L = length(perturb_vals)
	E_sensitivity = zeros(L, k, T);
	g_sensitivity = zeros(L, k, T)
		
	for j in 1:L
		d_crt = deepcopy(d_dyn)
		d_crt[cons_time][node_sens] = perturb_vals[j]
		opf_ = DynamicPowerManagementProblem(θ_dyn, d_crt)
		solve!(opf_, OPT, verbose=false)
		if ~(opf_.problem.status in [Convex.MOI.OPTIMAL, Convex.MOI.ALMOST_OPTIMAL])
			@show opf_.problem.status
		end
		
		for t in 1:T
			g_sensitivity[j, :, t] = evaluate(opf_.g[t])
			E_sensitivity[j, :, t] = evaluate(opf_.g[t]).*emissions_rates
		end
	end
end

# ╔═╡ 83c5e43f-585d-4588-941a-7d336d8ddd13
begin
	plot(lmes[node_sens, :, cons_time])
	xlabel!("Emission Time")
	ylabel!("LME")
	title!("LME of node $(node_sens) at time $(cons_time)")
end

# ╔═╡ e2f33432-cccd-4d00-b108-fe8ca7fd8d64
begin
	norm_E = sum(E_sensitivity[:, :, e_time], dims=2)./sum(E_sensitivity[idx_ref, :, e_time])
	plt_E_tot = plot(
		x_axis_vals, norm_E, 
		ylim=(0.95, 1.05), 
		linewidth=3
	)
	xlabel!(L"\Delta d/d")
	ylabel!(L"\Delta E/E")
	
	#adding the theoretical curve for the sensitivity
	E_th = (
		sum(E_sensitivity[idx_ref, :, e_time]) .+ (perturb_vals.-ref_val) .* lmes[node_sens, e_time, cons_time]
	)./sum(E_sensitivity[idx_ref, :, e_time])

	
	plot!(x_axis_vals, E_th, ls=:dash, c=:orange, linewidth=3)

	scatter!([x_axis_vals[idx_ref]], [norm_E[idx_ref]], c=:orange)

	title!("Total emissions at time $e_time")
end

# ╔═╡ Cell order:
# ╟─c39005df-61e0-4c08-8321-49cc5fe71ef3
# ╟─0f9bfc53-8a1a-4e25-a82e-9bc4dc0a11fc
# ╟─44275f74-7e7c-48d5-80a0-0f24609ef327
# ╠═db59921e-e998-11eb-0307-e396d43191b5
# ╠═0aac9a3f-a477-4095-9be1-f4babe1e2803
# ╠═0107b2ad-5e74-49ce-9231-26caed98c5c3
# ╠═0f7ed93f-72bf-41de-ba6d-f46926d9fc46
# ╠═a32d6a56-8da8-44b0-b659-21030692630a
# ╠═935dabe1-467f-4c36-bdff-4cb6807b672f
# ╟─9bd515d4-c7aa-4a3d-a4fb-28686290a134
# ╟─1bd72281-4a7f-44f4-974d-632e9d0aaf28
# ╠═0c786da1-7f44-40af-b6d6-e0d6db2242b2
# ╟─cfcba5ad-e516-4223-860e-b1f18a6449ba
# ╟─75dfaefd-abec-47e2-acc3-c0ff3a01048e
# ╠═c31093bf-4bd6-4642-81cf-377895799de7
# ╠═d225147d-658b-4f3e-87ca-be7c33aefe33
# ╠═d6c7b48a-bb28-413d-9ebc-f3675351c302
# ╠═e2bb5fc1-32e1-4dd5-95ec-f39c4ee7fa14
# ╠═40095334-a083-4abc-945d-d1efa1d285e8
# ╠═fcd7896a-cc91-4b97-9210-f87362965cbd
# ╠═f37f2ec5-b3a4-412c-9c10-da2cde64ec38
# ╠═806819d5-7b40-4ca4-aa1e-f1cf0a9a7f3f
# ╠═4fbf7a22-040f-4558-8a41-b581d2b98e94
# ╠═7d4b9df7-2181-4af7-a947-9a81297fe907
# ╠═a8d82fd9-b7bd-4f8a-9763-7ed684766e58
# ╟─25659734-b76d-4498-826a-18b5d751b37c
# ╠═fb910c5e-5b81-4103-9829-efab4f920f54
# ╟─edabacdd-8d25-4d64-9d4a-ecf1263ac02e
# ╠═380103e3-d893-4fc1-ae38-110a21cd863d
# ╠═4aed3df5-441b-445b-9277-a38690eb8603
# ╠═91f7d63c-9e30-4fd4-ab39-9fbf58d101dc
# ╟─83c5e43f-585d-4588-941a-7d336d8ddd13
# ╟─e2f33432-cccd-4d00-b108-fe8ca7fd8d64
