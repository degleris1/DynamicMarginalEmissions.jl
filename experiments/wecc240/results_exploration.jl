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

# ╔═╡ 0ae79723-a5cf-4508-b41d-9622948185a9
using Pkg; Pkg.activate("")

# ╔═╡ 5303b439-2bbb-4a04-b17e-7df6f2983493
using CSV, DataFrames, TOML

# ╔═╡ 668445dc-2437-421f-9251-b4044e5849f6
# We need all of these (including the weird ones like PooledArrays) 
# to make the BSON file load
using Convex, InlineStrings, PooledArrays, MathOptInterface, BSON

# ╔═╡ 64f8e88a-dfbf-4d25-b40e-af688e9e9f00
using SparseArrays

# ╔═╡ 32e5f26a-9b2f-4fc0-a0cd-1a5f101f0db9
using StatsBase: mean

# ╔═╡ e19f3dbe-b54a-45c3-b496-cf762f821ed5
using Statistics

# ╔═╡ 31aedaa9-2d5c-4ddf-acfa-8b836a252f70
using CarbonNetworks

# ╔═╡ 4d52befa-c24c-4ad0-b975-376b8c8af3d2
using Dates

# ╔═╡ 92a7dddf-a46e-4d42-9ca2-4d8b2e891b50
using Plots

# ╔═╡ a4577627-0a2b-403d-8558-4ccbf0622749
using LinearAlgebra

# ╔═╡ 620e86a8-60ab-4aa3-b313-29ab83ef5f4e
using PlutoUI

# ╔═╡ d1eef849-92ea-49ed-b05e-c4e055a85c2c
using LinearRegression

# ╔═╡ cab761be-ee1b-4002-a187-df72c29d9771
# equivalent to include that will replace it

function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end;

# ╔═╡ 113b99d8-0708-4da8-a8d6-7c60734e4a31
util = ingredients("util.jl")

# ╔═╡ f354a0e5-5e4f-40ab-807d-f9f19e74d378
analysis = ingredients("analysis_utils.jl")

# ╔═╡ f13e5dea-33b9-45c0-876e-42654fe6a8c7
util.FUEL_EMISSIONS

# ╔═╡ 6db70f24-e8ba-461e-8d86-00e9a37b44d3
md"""
## Load data
"""

# ╔═╡ f9fab4fe-baec-4bfd-9d84-ef9caac85f5f
config = TOML.parsefile(joinpath(@__DIR__, "../../config.toml"))

# ╔═╡ d7598abb-2be7-4e3b-af9e-14827ef5a3b0
DATA_DIR = config["data"]["GOOGLE_DRIVE"]

# ╔═╡ 45c73bb3-eecf-4b92-8e91-6a4c83addfdc
RESULTS_DIR = config["data"]["SAVE_DIR"]

# ╔═╡ efc7f458-1bbd-45e9-b05a-a076eee7a90c
run_names = ["july18_static", "july18_dynamic"]

# ╔═╡ 00c39925-e3a4-4d5e-9a00-e950c711648e
paths = [joinpath(DATA_DIR, "results240", name) for name in run_names]

# ╔═╡ e055e21b-166b-41f6-83d7-1ebf5e2626da
cases = [BSON.load(joinpath(p, "case.bson"), @__MODULE__) for p in paths];

# ╔═╡ ede8ffc2-2a2a-4a09-b702-610270f24491
results = map(analysis.load_results, paths);

# ╔═╡ f2197b5f-562a-435a-88a4-fc26415e77a9
n_nodes, _ = size(cases[1][:A]);

# ╔═╡ 24a57eb1-8761-4e03-9f36-3f5823ca703d
md"""
## Todo
- net demand (i.e. substract the renewable generation at each node from demand)
- hour of day
"""

# ╔═╡ 2e5e38fc-2a85-4cdb-9f30-b5100503d10c
md"""
## Total MEF for the dynamic model
"""

# ╔═╡ 083f1edc-06cf-413e-b187-df6d23b51850
alldates = d -> true

# ╔═╡ f0999280-9383-41ba-803b-840dcca60de5
dates = sort(collect(keys(results[2])))

# ╔═╡ 81898fc0-be77-4941-848e-581fb8cfb1ed
# collect total mefs
mefs = analysis.get_nodal_mefs(results[2], alldates; hybrid_mode=false);

# ╔═╡ 38c99b5d-a25a-4c0e-8908-961247d7e722
mefs_obs = analysis.get_nodal_mefs(results[2], alldates; hybrid_mode=false, observed=true);

# ╔═╡ dfbf46a2-8a02-4c77-b757-6a7546a36f1f
hours = [mod(k, 24) for k in 1:size(mefs)[2]];

# ╔═╡ c8c8d968-7e5b-4edf-8fec-2a7625b119a3
# representing total mefs
begin
	plot(transpose(mefs), color=:gray, legend=false)
	ylims!(-100, 1e4)
	xlabel!("Time [hr]")
	ylabel!("MEF []")
	title!("Total mef over time for each node")
end

# ╔═╡ f5e407e6-fa2c-4825-a7fe-6690bf6ad06c
# representing observed mef
begin
	plot(transpose(mefs_obs), color=:gray, legend=false)
	ylims!(-100, 1e4)
	xlabel!("Time [hr]")
	ylabel!("MEF []")
	title!("Observed mef over time for each node")
end

# ╔═╡ b7a82318-2999-42c0-975f-079b28104684
md"""
### For a given node: 
"""

# ╔═╡ 0fdee674-965a-4f9b-9f52-aff2efe935f8
node = 65

# ╔═╡ 8d206091-89b8-49c3-af53-95a262298426
begin
	p_ = plot()
	plot!(mefs[node, :], label="true")
	plot!(mefs_obs[node, :], label="observed")
	xlabel!("Time")
	ylabel!("MEF")
	ylims!(-100, 2000)
	xlims!(200, 400)
	p_
end

# ╔═╡ 7f1ebbfc-4091-4efc-bc56-151e91991bc2
size(mefs[node, :])

# ╔═╡ 7bbd7076-5a16-4914-8118-0d9238f75569
begin
	mefs_node_days = reshape(mefs[node, :], (24, Int(length(mefs[node, :])/24)))
	mean_ = mean(mefs_node_days, dims=2)
	std_ = std(mefs_node_days, dims=2);
	plot(mean_ - std_, fillrange=mean_ .+ std_, fillalpha = .3, label="std range")
	plot!(mean_, label="mean", lw=2)
	xlabel!("Hour of day")
	ylabel!("MEF")
end

# ╔═╡ e8fee23c-dc1a-430a-b896-31f09c4c503e
md"""
## Compute regression-based and observed MEFs
"""

# ╔═╡ 1bed3bbb-7411-441c-8028-f4dc91b0f1aa
md"""
confusion: 
- observed MEF can be defined as the diagonal of the mef matrix
- then can also be defined as the approximation by linear regression at that given hour

Can you extract and analyze both of them? 
"""

# ╔═╡ a1da5ecd-219d-4e4e-b39d-1addddc89f13
demand = analysis.get_demand(results[2]);

# ╔═╡ 38b8e954-6e93-4f2d-9a87-2ebbab8a4071
Δd = hcat([demand[:, i] - demand[:, i-1] for i in 2:size(demand)[2]]...);

# ╔═╡ 576e6a18-eb15-4578-9bda-569c1397a8d3
co2_rates = cases[2][:co2_rates];

# ╔═╡ 112eb399-f273-42a0-a6b6-52ce543bc696
E = analysis.get_total_emissions(results[2], transpose(co2_rates));

# ╔═╡ a003c654-7072-4e3f-9aa8-154c0efa565c
begin
	plot(transpose(E), legend=false)
	xlabel!("Time")
	ylabel!("Total emissions across the network")
end

# ╔═╡ fb6a6791-c704-43f1-9784-1cd17dfa7858
ΔE = [E[i]-E[i-1] for i in 2:length(E)];

# ╔═╡ fd198417-f955-443f-b7bd-40b5106356eb
hour_filt1 = 10 .<= hours[2:end] .<= 20;

# ╔═╡ eece6a95-ee8a-466f-886b-ac3a42d3f91a
.!hour_filt1

# ╔═╡ 09b71d83-53ad-4b7b-ba70-b51ea90aa37e
begin
	scatter(Δd[node, hour_filt1], ΔE[hour_filt1], legend=false)
	scatter!(Δd[node, .!hour_filt1], ΔE[.!hour_filt1], legend=false)
	# mef_node = mefs_reg[node]
	# xx = LinRange(minimum(Δd[node, :]), maximum(Δd[node, :]), 10)
	# plot!(xx, mef_single_node*xx, lw=3)
	
	xlabel!("Δd at node $(node)")
	ylabel!("ΔE")
	title!("ΔE total wrt Δd at node $(node)")
end

# ╔═╡ ec0278e1-0300-4842-80f0-857b414867a9
μ = 1e-4

# ╔═╡ 66ecb9fd-b976-432a-9df3-7a6a8618bc2e
# linear regression on all nodes at once
begin
	X = transpose(Δd); # matrix has many nodes with zero demand; therefore columns have to be filtered so as to use only the nodes with actual mef
	demand_nodes = findall(vec(sum(abs.(X), dims=1)).!==0.0)
	colin_node = 64 # node that is a lin comb of others
	demand_nodes = demand_nodes[demand_nodes .!= colin_node]
	X_demand = X[:, demand_nodes]
	rank_ = rank(X_demand)
	@assert rank_ == size(X_demand)[2]
	mef_opt = inv(X_demand'*X_demand + μ*I(rank_))*X_demand'*ΔE
	mefs_reg = zeros(n_nodes)
	mefs_reg[demand_nodes] .= mef_opt
end;

# ╔═╡ 364f3ff0-75b7-4ba2-a82a-08e92171daa2
# linear regression on a single node
begin
	XX = transpose(Δd)[:, node]; # matrix has many nodes with zero demand; therefore columns have to be filtered so as to use only the nodes with actual mef
	mef_single_node = inv(transpose(XX)*XX)*transpose(XX)*ΔE
end

# ╔═╡ 7db52f69-03f7-4c9d-b746-4312b9c038b2
md"""
The MEF computed by linear regression at node $(node) is $(mefs_reg[node]). 

The MEF if computed on a single node at the same node is $(mef_single_node)
"""

# ╔═╡ 3dc3b990-1e8b-43b7-856a-17667c6be634
# histogram
begin
	bins = -1000:100:3000
	histogram(mefs[node, :], bins=bins, alpha=0.4, label="True")
	histogram!(mefs_obs[node, :], bins=bins, alpha=0.4, label="Observed")
	plot!([mefs_reg[node]], seriestype="vline", label="Regression", lw=3)
	plot!([mean(mefs[node, :])], seriestype="vline", label="Mean True", lw=3)
	plot!([mean(mefs_obs[node, :])], seriestype="vline", label="Mean Observed", lw=3)
	xlabel!("MEF")
	ylabel!("Counts")
end

# ╔═╡ ad960911-9354-4a31-bcbc-f1f0a95482ae
md"""
### Analysis of regression-based mefs estimation
"""

# ╔═╡ 7833642f-5b87-4d1d-b311-72901851a58f
let
	u = median(mefs, dims=2)[demand_nodes]
	p = sortperm(u)
	scatter(u[p], label="median mef true")
	scatter!(median(mefs_obs, dims=2)[demand_nodes][p], label="median mef obs")
	scatter!(mefs_reg[demand_nodes][p], label="mef reg")
	ylims!(0, 2000)
	xlabel!("Node number")
	ylabel!("MEF")
	title!("MEFs of nodes, sorted according to true")
end

# ╔═╡ e42f39c5-9f4e-4e1c-bfdc-f7bc9cd29623
demand_nodes[35]

# ╔═╡ bed124c1-b6ea-4eb6-9d18-da41651d932d
md"""
### Estimate the prediction power of this regression. 
"""

# ╔═╡ 0defd4c1-dbb9-4b66-bca0-3c4605c252bb
md"""
### Compute regression-based mefs at a given hour
"""

# ╔═╡ e8c4a66c-94b7-4b0e-9b43-fcaa0f235733
@bind hr Slider(1:24)

# ╔═╡ 6b5f3833-6168-4d45-9283-0fec34404658
hr

# ╔═╡ f8ebd537-bacc-4e1a-af8f-12f1196536ff
ΔE_hr = [ΔE[hr + k*24] for k in 0:length(results[2])-1];

# ╔═╡ 827932a6-68b6-4b37-8f02-fdfa819fee2c
Δd_hr = hcat([Δd[:, hr + k*24] for k in 0:length(results[2])-1]...);

# ╔═╡ 1790b9e8-dc3c-4e47-8d8e-0a0cdea118fd
begin
	scatter(Δd_hr[node, :], ΔE_hr, legend=false)
	ylabel!("ΔE at given hour")
	xlabel!("Δx at given hour")
end

# ╔═╡ 667d9d4a-88d9-48e2-86bc-63f2a587f1f5
# linear regression on all nodes at once
begin
	X_hr = transpose(Δd_hr); # matrix has many nodes with zero demand; therefore columns have to be filtered so as to use only the nodes with actual mef
	X_hr = X_hr[:, demand_nodes]
	rank_hr = rank(X_hr)
	@assert rank_hr == size(X_hr)[2]
	mef_hr_opt = inv(X_hr'*X_hr + μ*I(rank_hr))*X_hr'*ΔE
	mefs_hr_reg = zeros(n_nodes)
	mefs_hr_reg[demand_nodes] .= mef_hr_opt
end;

# ╔═╡ 49c340c0-64e1-4d82-a24a-f7b83c108c64
rank_hr

# ╔═╡ dc0bd18b-2276-4a05-a7b1-f43acaf21b85
size(X_hr)

# ╔═╡ 56bd49af-7c69-4778-be5b-c47e906b69ed
rank_hr

# ╔═╡ f8f3ae28-e46b-43fb-bc06-a799ea0f6cdf
Δd_hr

# ╔═╡ 305bdb8a-e186-4c5a-b927-b7a406ac260a
md"""
### Compare regression-based MEF vs exact differentiation-based MEFs
"""

# ╔═╡ 701285fc-887d-425c-98a3-1ee09235066f
md"""
Todo: 
- demand map: estimate mef bsed on region-level demand and then map to nodal, compare with results from comptuations
"""

# ╔═╡ 1d5d1c1c-050a-434e-a977-b8d2aea69b26
# get the co2 costs if they were not saved
co2_costs = util.get_costs(data[:case].heat, data[:case].fuel, util.FUEL_EMISSIONS);

# ╔═╡ a34e32d8-3e18-4c22-87e2-032360661498
g_opt = [[results_dyn[d].g[k].value for k in 1:24] for d in 1:365];

# ╔═╡ fa51de9f-6f7a-4173-96d4-18f5f225aa1a
E_tot = [[dot(co2_costs,g_opt[d][k]) for k in 1:24] for d in 1:365];

# ╔═╡ be2f82a2-5c22-41a3-abe4-5a6f9de4e6d7
md"""
We want to compute, at a given hour, the changes in emissions and related them to changes in demand at a given node.
"""

# ╔═╡ 82ad33b1-3719-4aeb-9c82-1f8f1812ed61
node1, node2

# ╔═╡ 7213814b-3dd9-4c8d-9f85-947a33d96e44
# recontsruct the flows over the network
function get_ps(case, results)

	F = CarbonNetworks.make_pfdf_matrix(case.A, case.β)

	ps = []

	for k in 1:length(results)
		p_crt = []
		for j in 1:length(results[k].g)
			g_crt = results[k].g[j]
			d_crt = results[k].d[j]
			p = CarbonNetworks.evaluate(F*(d_crt - case.B*g_crt))
			push!(p_crt, p)
		end
		push!(ps, p_crt)
	end

	return ps
end

# ╔═╡ ff851aef-904f-444a-afb4-f7ef7d7766c4
ps_dyn = get_ps(data_dyn[:case], results_dyn);

# ╔═╡ d9f9dda4-b559-44a4-9677-fbe967322b2b
# get pmax - HARDCODED for now
begin
line_weight = 2.0
line_max = 100.0
Z = 1e3
pmax = line_weight * min.(data_dyn[:case].fmax / Z, line_max);
end;

# ╔═╡ 2ac72e25-e6f8-4468-be71-f64b3797a96c
# gets the number of congested lines at each timestep
congested = [[sum(abs.(ps_dyn[d][h]./pmax) .> .99) for h in 1:24] for d in 1:365];

# ╔═╡ df45b862-8d94-4b71-951c-5d84a7d16af2
# demand at the given hour over the year
d_h = hcat([results_dyn[d].d[hour] for d in 1:365]...);

# ╔═╡ 087a3ba0-14f4-4373-99f1-3dd2ccdb71b9
gmax = hcat([results_dyn[d].gmax[k] for k in 1:24 for d in 1:365]...)

# ╔═╡ 5ba56789-c1bb-4c1b-95cc-4943561434f4
md"""
Plotting all the generator data. 
"""

# ╔═╡ 44674366-5546-44b7-b6cc-b2476e3c8fc6
begin
	f = Figure()
	ax = Axis(f[1,1], xlabel="time", ylabel="gmax/gmax(t=0)")
	ng, _ = size(gmax)
	for gid in 1:ng
		lines!(gmax[gid, :]./gmax[gid, 1], color="gray")
	end
	f
end


# ╔═╡ 7ae73ea6-8878-43cc-9207-fda58fa41d9c
df = util.load_wecc_240_dataset();

# ╔═╡ 812a305a-c9bb-4b57-9621-cadf4b70984c
# what we would expect for the month of jan

let
gmax_ = []
for hour in 1:24
	for day in 1:1
		for month in 1:1
demand_map = util.get_demand_map(hour, day, month, 2004, df.demand)
node_names, node_ids = util.get_node_info(df.branch)

B, gmin, gmax, ramp, heat, fuel = util.get_generator_data(demand_map, node_ids, df.gen, df.heat)

push!(gmax_, gmax)
end
end
end
	
gmax_ = hcat(gmax_...)
	
f = Figure()
ax = Axis(f[1,1], xlabel="time", ylabel="gmax/gmax(t=0)")
ng, _ = size(gmax_)
for gid in 1:ng
lines!(gmax_[gid, :]./gmax_[gid, 1], color="gray")
end
f

end

# ╔═╡ 91b474ba-e955-416c-8bc7-aa2d9b88995f
md"""

## regression-based mef estimation
"""

# ╔═╡ e1b5c93e-9241-442b-a8ce-5c7d91809efc
md"""
### Regressing over each node
"""

# ╔═╡ 847333af-8145-433a-b24c-554af2468da4
md"""
Here we want to estimate the MEF from each node. The problem is that demand is proportionally divided, therefore a regression analysis would probably have mefs lie on an affine set and easily attribute negative MEFs.
"""

# ╔═╡ 13478bcb-c4fc-4532-a1ef-4513b15e3295
# # changes in emissions
# ΔE = [(E_tot[d][hour] - E_tot[d][hour-1]) for d in 1:365];

# ╔═╡ 70723775-1912-48ed-9ae5-d4663d0f81d3
# gather all the demand changes, at every node
# Δd = hcat([
# 	[results_dyn[d].d[hour][n] - results_dyn[d].d[hour-1][n] for d in 1:365] 
# 	for n in 1:length(nodes)
# ]...);

# ╔═╡ 50acbd4b-1a02-4c55-b605-caf07f12bd74
MEFs_reg = LinearRegression.slope(lr_E)

# ╔═╡ 749ef6df-5909-4f40-a5ff-0ed7877de9ab
lines(MEFs_reg)

# ╔═╡ eb4e9485-5365-4c39-b791-3a2137679280
md"""
we see that many nodes have a negative average MEF, which does not make sense from scatter plots
"""

# ╔═╡ 5e77a674-9ece-415e-abd5-2b244e47beb9
md"""
Regressing over demand regions.
"""

# ╔═╡ 870d811d-ece6-4513-8f9a-558585d0d70a
df.demand

# ╔═╡ 18941712-b800-48c3-ad09-b16a0e5a8f9b
demand_map = util.get_demand_map(1, 1, 1, 2004, df.demand)

# ╔═╡ 8d7c5fd9-5325-481a-a572-ca4e8dcc8655
regions = unique(Array(df.participation[:, "Region"]))

# ╔═╡ 2238db2f-3860-4abb-938a-8b771b0ee6d8
# we want to create a vector of demands for all the regions only
begin

	demands_regions = []

	for k in 1:size(df.demand)[1]
		yr = df.demand[k, "Year"]
		month = df.demand[k, "Month"]
		day = df.demand[k, "Day"]
		hr = df.demand[k, "Period"]

		demand_map = util.get_demand_map(hr, day, month, yr, df.demand)

		d_regions = [demand_map[rr] for rr in regions]

		push!(demands_regions, d_regions)
		
	end
	demands_regions = hcat(demands_regions...)
end;

# ╔═╡ 4a1b69b6-1567-4199-9a32-b2019e88366c
begin
	_, nts = size(demands_regions)
	ndays = Int(nts/24)
	idx_hr = [hour+k*24 for k in 0:ndays-1]
end

# ╔═╡ 6f68dfa5-7f43-4ed0-baf0-683b8c14cc5b
begin
	Δd_regions = transpose(hcat([demands_regions[:, idx] - demands_regions[:, idx-1] for idx in idx_hr]...));
	Δd_regions = Δd_regions[1:end-1, :]; # we remove the last day as we assume we did not compute anything for it
end

# ╔═╡ 81a3439e-0a5c-41ab-b4fa-9bda14451282
lr_regions = linregress(Δd_regions, ΔE);

# ╔═╡ cb04852d-6ccf-46b2-94cb-e5faa3f4cb7a
mefs_regions = LinearRegression.slope(lr_regions)

# ╔═╡ c22d2cc0-8ea3-4a06-8283-9af4d13d78eb
lines(mefs_regions)

# ╔═╡ 5cce0ff7-b8ae-4288-84b4-20e2fc50b0b4
region_id = 7

# ╔═╡ 85a13faf-79ce-4e40-a69e-48f65dc2a7e9
begin
	
	Δd_crt = Δd_regions[:, region_id]
	x_crt = LinRange(1, maximum(Δd_crt), 100)
	p=scatter(Δd_crt, ΔE)
	slope_reg = mefs_regions[region_id]
	lines!(x_crt, slope_reg*x_crt, color=:black)
	p
end

# ╔═╡ 34cd82e7-b58f-469e-be91-a1391903d661
md"""
## NOTE: the linear regression really does not work... to inquire. 
"""

# ╔═╡ 90aa489d-3c50-4c4a-ac19-351de901fbe4
md"""
## add subtitle
"""

# ╔═╡ 28591eda-995f-4224-98b8-ca1eab27559e
# only the diagonal terms in the mefs
all_λs_obs_dyn = [[diag(results_dyn[d].λ[n, :, :]) for d in 1:365] for n in 1:length(nodes)];

# ╔═╡ 8e1f0bef-459d-4598-a01a-e59e00f53247
# collect the average observable mef from the dynamic model
begin

	λ_obs_1 = [all_λs_obs_dyn[node1][d][hour] for d in 1:365]
	λ_obs_2 = [all_λs_obs_dyn[node2][d][hour] for d in 1:365]

	λ_total_1 = all_λs_dyn[:, node1]
	λ_total_2 = all_λs_dyn[:, node2]
end;

# ╔═╡ 34253cd5-3049-4651-a5f6-06807a2233bd
let
	f = Figure()
	ax = Axis(f[1,1], xlabel="day", ylabel="number of congested lines")

	lines!([congested[d][hour] for d in 1:365])
	lines!(λ_total_1./1000)

	ax = Axis(f[1,2], xlabel="total demand", ylabel="number congested lines")

	scatter!(vec(sum(d_h, dims=1)), [congested[d][hour] for d in 1:365])
	# ylims!(ax, -2000, 2000)

	ax = Axis(f[2,2], xlabel="# congested lines", ylabel="λ true")

	scatter!([congested[d][hour] for d in 1:365], λ_total_1)
	ylims!(ax, -2000, 2000)


	ax = Axis(f[2,1], xlabel="# congested lines", ylabel="observed mef")

	scatter!([congested[d][hour] for d in 1:365], λ_obs_1)
	ylims!(ax, -2000, 2000)
	

	
	f

end

# ╔═╡ 6f2896ac-ecf6-4256-82e0-0a4070fa14af
begin
	f2 = Figure()
	hist(f2[1, 1], λ_obs_1, bins = -1000:200:1000)
	hist!(f2[1, 1], λ_total_1, bins = -1000:200:1000)
	f2
end

# ╔═╡ 982d36b3-a7b6-4066-9441-9e02596423dd
let
f = Figure()

ax = Axis(f[1, 1], xlabel = "x label", ylabel = "y label",
    title = "Title")
for n = 1:length(nodes)
	lines!(all_λs_dyn[:, n], color=:gray)
end
ylims!(ax, -1000, 4000)
ax.title="All total MEFs over time"
ax.xlabel="Day"
ax.ylabel="MEF"
f
end

# ╔═╡ e834a7d2-5dd9-4995-b3c9-3c771673edaa
let
f = Figure()

ax = Axis(f[1, 1], xlabel = "x label", ylabel = "y label",
    title = "Title")

lines!(λ_obs_1)
lines!(λ_total_1, color=:red)
ylims!(ax, -2000, 2000)
f
end

# ╔═╡ b886963f-e1a9-4782-918c-f5009e57abc5
let
	f = Figure()
	ax = Axis(f[1,1], xlabel="nodal demand", ylabel=" observed mef")

	scatter!(d_h[node1, :], λ_obs_1)
	ylims!(ax, -2000, 2000)

	ax = Axis(f[1,2], xlabel="nodal demand", ylabel="true mef")

	scatter!(d_h[node1, :], λ_total_1)
	ylims!(ax, -2000, 2000)


	ax = Axis(f[2,1], xlabel="total demand", ylabel="observed mef")

	scatter!(vec(sum(d_h, dims=1)), λ_obs_1)
	ylims!(ax, -2000, 2000)
	
	ax = Axis(f[2,2], xlabel="total demand", ylabel="true mef")

	scatter!(vec(sum(d_h, dims=1)), λ_total_1)
	ylims!(ax, -2000, 2000)
	
	f

end

# ╔═╡ bd656131-eb56-467d-8f98-5e8de88266c3
md"""
There clearly are two populations. Therefore you want to: 
- classify the points wrt them being in group 1 or 2
- run different linear regressions based on the group they're in
- compare those different regressions.

the very simple answer for now is to do: 
- smaller than X in group 0
- larger than X is group 1
"""

# ╔═╡ e91052d7-1ca0-4ff5-aae4-1b8ace8f93cf
begin
	cutoff1 = 500
	assign_group(cutoff, x) = x > cutoff ? 1 : 0;

	group1 = assign_group.(cutoff1, λ_obs_1)
end

# ╔═╡ 6377fd8f-7b2b-4720-bf3a-a542075bcedd
lines(group1)

# ╔═╡ 2dd6a4cb-42ab-4c82-be5e-529c763059c2
Bool.(group1)

# ╔═╡ 4d032e76-afe0-4b6b-bed2-bf1df1362dc5
md"""
Rolling window in julia? 
"""

# ╔═╡ 28664572-0078-434f-a4cb-d70eaafe8d4d
begin
lr1_group1 = linregress(xx1[Bool.(group1)], ΔE[Bool.(group1)])
lr1_group0 = linregress(xx1[Bool.(1 .-group1)], ΔE[Bool.(1 .-group1)])
end;

# ╔═╡ 9ca0036b-e68f-40ab-abc7-32fa346c01da
LinearRegression.slope(lr1_group1)

# ╔═╡ 08c75daa-fe09-4839-9390-cd5034eb8125
LinearRegression.slope(lr1_group0)

# ╔═╡ f9153c23-f969-4854-8c92-2437efd5b8ce
begin
	s1 = scatter(xx1[Bool.(group1)], ΔE[Bool.(group1)])
	scatter!(xx1[Bool.(1 .-group1)], ΔE[Bool.(1 .-group1)])
	s1
end

# ╔═╡ 58af890e-c198-4766-a7e3-3024dafe3190
p2=scatter(xx2, ΔE)

# ╔═╡ Cell order:
# ╠═0ae79723-a5cf-4508-b41d-9622948185a9
# ╠═5303b439-2bbb-4a04-b17e-7df6f2983493
# ╠═668445dc-2437-421f-9251-b4044e5849f6
# ╠═64f8e88a-dfbf-4d25-b40e-af688e9e9f00
# ╠═32e5f26a-9b2f-4fc0-a0cd-1a5f101f0db9
# ╠═e19f3dbe-b54a-45c3-b496-cf762f821ed5
# ╠═31aedaa9-2d5c-4ddf-acfa-8b836a252f70
# ╠═4d52befa-c24c-4ad0-b975-376b8c8af3d2
# ╠═92a7dddf-a46e-4d42-9ca2-4d8b2e891b50
# ╠═a4577627-0a2b-403d-8558-4ccbf0622749
# ╟─cab761be-ee1b-4002-a187-df72c29d9771
# ╠═113b99d8-0708-4da8-a8d6-7c60734e4a31
# ╠═f354a0e5-5e4f-40ab-807d-f9f19e74d378
# ╟─f13e5dea-33b9-45c0-876e-42654fe6a8c7
# ╟─6db70f24-e8ba-461e-8d86-00e9a37b44d3
# ╠═f9fab4fe-baec-4bfd-9d84-ef9caac85f5f
# ╠═d7598abb-2be7-4e3b-af9e-14827ef5a3b0
# ╠═45c73bb3-eecf-4b92-8e91-6a4c83addfdc
# ╠═efc7f458-1bbd-45e9-b05a-a076eee7a90c
# ╠═00c39925-e3a4-4d5e-9a00-e950c711648e
# ╠═e055e21b-166b-41f6-83d7-1ebf5e2626da
# ╠═ede8ffc2-2a2a-4a09-b702-610270f24491
# ╠═f2197b5f-562a-435a-88a4-fc26415e77a9
# ╟─24a57eb1-8761-4e03-9f36-3f5823ca703d
# ╟─2e5e38fc-2a85-4cdb-9f30-b5100503d10c
# ╠═083f1edc-06cf-413e-b187-df6d23b51850
# ╠═f0999280-9383-41ba-803b-840dcca60de5
# ╠═81898fc0-be77-4941-848e-581fb8cfb1ed
# ╠═38c99b5d-a25a-4c0e-8908-961247d7e722
# ╠═dfbf46a2-8a02-4c77-b757-6a7546a36f1f
# ╟─c8c8d968-7e5b-4edf-8fec-2a7625b119a3
# ╟─f5e407e6-fa2c-4825-a7fe-6690bf6ad06c
# ╟─b7a82318-2999-42c0-975f-079b28104684
# ╠═0fdee674-965a-4f9b-9f52-aff2efe935f8
# ╠═8d206091-89b8-49c3-af53-95a262298426
# ╠═7f1ebbfc-4091-4efc-bc56-151e91991bc2
# ╟─7bbd7076-5a16-4914-8118-0d9238f75569
# ╟─e8fee23c-dc1a-430a-b896-31f09c4c503e
# ╟─1bed3bbb-7411-441c-8028-f4dc91b0f1aa
# ╠═a1da5ecd-219d-4e4e-b39d-1addddc89f13
# ╠═38b8e954-6e93-4f2d-9a87-2ebbab8a4071
# ╠═576e6a18-eb15-4578-9bda-569c1397a8d3
# ╠═112eb399-f273-42a0-a6b6-52ce543bc696
# ╟─a003c654-7072-4e3f-9aa8-154c0efa565c
# ╠═fb6a6791-c704-43f1-9784-1cd17dfa7858
# ╠═fd198417-f955-443f-b7bd-40b5106356eb
# ╠═eece6a95-ee8a-466f-886b-ac3a42d3f91a
# ╠═09b71d83-53ad-4b7b-ba70-b51ea90aa37e
# ╟─7db52f69-03f7-4c9d-b746-4312b9c038b2
# ╠═ec0278e1-0300-4842-80f0-857b414867a9
# ╠═66ecb9fd-b976-432a-9df3-7a6a8618bc2e
# ╟─364f3ff0-75b7-4ba2-a82a-08e92171daa2
# ╠═3dc3b990-1e8b-43b7-856a-17667c6be634
# ╟─ad960911-9354-4a31-bcbc-f1f0a95482ae
# ╠═7833642f-5b87-4d1d-b311-72901851a58f
# ╠═e42f39c5-9f4e-4e1c-bfdc-f7bc9cd29623
# ╟─bed124c1-b6ea-4eb6-9d18-da41651d932d
# ╟─0defd4c1-dbb9-4b66-bca0-3c4605c252bb
# ╠═620e86a8-60ab-4aa3-b313-29ab83ef5f4e
# ╟─e8c4a66c-94b7-4b0e-9b43-fcaa0f235733
# ╟─6b5f3833-6168-4d45-9283-0fec34404658
# ╠═f8ebd537-bacc-4e1a-af8f-12f1196536ff
# ╠═827932a6-68b6-4b37-8f02-fdfa819fee2c
# ╠═1790b9e8-dc3c-4e47-8d8e-0a0cdea118fd
# ╠═49c340c0-64e1-4d82-a24a-f7b83c108c64
# ╠═dc0bd18b-2276-4a05-a7b1-f43acaf21b85
# ╠═56bd49af-7c69-4778-be5b-c47e906b69ed
# ╠═667d9d4a-88d9-48e2-86bc-63f2a587f1f5
# ╠═f8f3ae28-e46b-43fb-bc06-a799ea0f6cdf
# ╟─305bdb8a-e186-4c5a-b927-b7a406ac260a
# ╟─701285fc-887d-425c-98a3-1ee09235066f
# ╠═1d5d1c1c-050a-434e-a977-b8d2aea69b26
# ╠═a34e32d8-3e18-4c22-87e2-032360661498
# ╠═fa51de9f-6f7a-4173-96d4-18f5f225aa1a
# ╟─be2f82a2-5c22-41a3-abe4-5a6f9de4e6d7
# ╠═82ad33b1-3719-4aeb-9c82-1f8f1812ed61
# ╠═7213814b-3dd9-4c8d-9f85-947a33d96e44
# ╠═ff851aef-904f-444a-afb4-f7ef7d7766c4
# ╠═d9f9dda4-b559-44a4-9677-fbe967322b2b
# ╠═2ac72e25-e6f8-4468-be71-f64b3797a96c
# ╠═34253cd5-3049-4651-a5f6-06807a2233bd
# ╠═df45b862-8d94-4b71-951c-5d84a7d16af2
# ╠═087a3ba0-14f4-4373-99f1-3dd2ccdb71b9
# ╟─5ba56789-c1bb-4c1b-95cc-4943561434f4
# ╠═44674366-5546-44b7-b6cc-b2476e3c8fc6
# ╠═7ae73ea6-8878-43cc-9207-fda58fa41d9c
# ╟─812a305a-c9bb-4b57-9621-cadf4b70984c
# ╟─91b474ba-e955-416c-8bc7-aa2d9b88995f
# ╟─e1b5c93e-9241-442b-a8ce-5c7d91809efc
# ╟─847333af-8145-433a-b24c-554af2468da4
# ╠═13478bcb-c4fc-4532-a1ef-4513b15e3295
# ╠═70723775-1912-48ed-9ae5-d4663d0f81d3
# ╠═50acbd4b-1a02-4c55-b605-caf07f12bd74
# ╠═749ef6df-5909-4f40-a5ff-0ed7877de9ab
# ╟─eb4e9485-5365-4c39-b791-3a2137679280
# ╟─5e77a674-9ece-415e-abd5-2b244e47beb9
# ╠═870d811d-ece6-4513-8f9a-558585d0d70a
# ╠═18941712-b800-48c3-ad09-b16a0e5a8f9b
# ╠═8d7c5fd9-5325-481a-a572-ca4e8dcc8655
# ╠═2238db2f-3860-4abb-938a-8b771b0ee6d8
# ╠═4a1b69b6-1567-4199-9a32-b2019e88366c
# ╠═6f68dfa5-7f43-4ed0-baf0-683b8c14cc5b
# ╠═81a3439e-0a5c-41ab-b4fa-9bda14451282
# ╠═cb04852d-6ccf-46b2-94cb-e5faa3f4cb7a
# ╠═c22d2cc0-8ea3-4a06-8283-9af4d13d78eb
# ╠═5cce0ff7-b8ae-4288-84b4-20e2fc50b0b4
# ╠═85a13faf-79ce-4e40-a69e-48f65dc2a7e9
# ╠═34cd82e7-b58f-469e-be91-a1391903d661
# ╠═90aa489d-3c50-4c4a-ac19-351de901fbe4
# ╠═28591eda-995f-4224-98b8-ca1eab27559e
# ╠═8e1f0bef-459d-4598-a01a-e59e00f53247
# ╠═6f2896ac-ecf6-4256-82e0-0a4070fa14af
# ╠═982d36b3-a7b6-4066-9441-9e02596423dd
# ╠═e834a7d2-5dd9-4995-b3c9-3c771673edaa
# ╠═b886963f-e1a9-4782-918c-f5009e57abc5
# ╟─bd656131-eb56-467d-8f98-5e8de88266c3
# ╠═6377fd8f-7b2b-4720-bf3a-a542075bcedd
# ╠═e91052d7-1ca0-4ff5-aae4-1b8ace8f93cf
# ╠═2dd6a4cb-42ab-4c82-be5e-529c763059c2
# ╠═4d032e76-afe0-4b6b-bed2-bf1df1362dc5
# ╠═d1eef849-92ea-49ed-b05e-c4e055a85c2c
# ╠═28664572-0078-434f-a4cb-d70eaafe8d4d
# ╠═9ca0036b-e68f-40ab-abc7-32fa346c01da
# ╠═08c75daa-fe09-4839-9390-cd5034eb8125
# ╠═f9153c23-f969-4854-8c92-2437efd5b8ce
# ╠═58af890e-c198-4766-a7e3-3024dafe3190
