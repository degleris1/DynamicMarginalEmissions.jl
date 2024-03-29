### A Pluto.jl notebook ###
# v0.19.4

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

# ╔═╡ 8876dc66-7d6d-48ac-94d5-5dd489e81b87
using LaTeXStrings

# ╔═╡ 620e86a8-60ab-4aa3-b313-29ab83ef5f4e
using PlutoUI

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
run_names = ["july18_static", "july18_dynamic", "year18_dynamic"]

# ╔═╡ 00c39925-e3a4-4d5e-9a00-e950c711648e
paths = [joinpath(DATA_DIR, "results240", name) for name in run_names]

# ╔═╡ e055e21b-166b-41f6-83d7-1ebf5e2626da
cases = [BSON.load(joinpath(p, "case.bson"), @__MODULE__) for p in paths];

# ╔═╡ ede8ffc2-2a2a-4a09-b702-610270f24491
results = map(analysis.load_results, paths);

# ╔═╡ 11aef0e7-372d-4b3b-845d-e63bc620b5f9
idx = 3

# ╔═╡ f2197b5f-562a-435a-88a4-fc26415e77a9
n_nodes = size(cases[idx][:A], 1);

# ╔═╡ 6b3a9e01-e14c-488d-bde5-7eae4e8e87e0
n_gens = size(cases[idx][:B], 2)

# ╔═╡ ae72be0b-b54f-47bd-9888-ad297463bb38
md"""
We are analyzing the dataset: $(run_names[idx])
"""

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
dates = sort(collect(keys(results[idx])))

# ╔═╡ 81898fc0-be77-4941-848e-581fb8cfb1ed
# collect total mefs
mefs = analysis.get_nodal_mefs(results[idx], alldates; hybrid_mode=false);

# ╔═╡ 38c99b5d-a25a-4c0e-8908-961247d7e722
mefs_obs = analysis.get_nodal_mefs(results[idx], alldates; hybrid_mode=false, observed=true);

# ╔═╡ dfbf46a2-8a02-4c77-b757-6a7546a36f1f
hours = [h for _ in 1:length(dates) for h in 1:24];

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
node = 200

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

# ╔═╡ 7bbd7076-5a16-4914-8118-0d9238f75569
begin
	mefs_nodes_days = [reshape(mefs[n, :], (24, Int(length(mefs[n, :])/24))) for n in 1:n_nodes]
	mean_ = mean(mefs_nodes_days[node], dims=2)
	std_ = std(mefs_nodes_days[node], dims=2);
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
demand = analysis.get_demand(results[idx]);

# ╔═╡ 66caf393-88ce-4e6c-9a4a-afdfa7b24304
Δd = diff(demand, dims=2);

# ╔═╡ 576e6a18-eb15-4578-9bda-569c1397a8d3
co2_rates = cases[2][:co2_rates];

# ╔═╡ 112eb399-f273-42a0-a6b6-52ce543bc696
E = analysis.get_total_emissions(results[idx], transpose(co2_rates));

# ╔═╡ 9b4e015d-a0b6-4046-b9c7-9f12b93d84fa
length(sort(collect(keys(results[idx]))))

# ╔═╡ 98e1f152-67c7-43d2-b3b5-68bbc322478d
plot([sum(demand[:, i]) for i in 1:size(demand, 2)], label=false, xlabel="Time", ylabel="Demand")

# ╔═╡ a003c654-7072-4e3f-9aa8-154c0efa565c
begin
	plot(E, legend=false)
	xlabel!("Time")
	ylabel!("Total emissions across the network")
end

# ╔═╡ fb6a6791-c704-43f1-9784-1cd17dfa7858
ΔE = diff(E, dims=1);

# ╔═╡ 9323dbf3-6a22-476a-8186-766fe1515c18
function get_R2(preds, y)
	res = preds - y 
	R2 = 1 - res'*res/sum((y .- mean(y)).^2)
	return R2
end;

# ╔═╡ 26d698fe-c2da-4a34-b446-4dd4dfacf861
function plot_preds(preds, y)
	p = scatter(preds, y, legend=false, alpha=.2)
	plot!(
		LinRange(minimum(y), maximum(y), 10), 
		LinRange(minimum(y), maximum(y), 10), 
		ls=:dash, lw = 3
	)
	ylabel!("ΔE real")
	xlabel!("ΔE predicted")
	return p
end;

# ╔═╡ 01805ad9-0d4c-4d8e-8f79-69da710a3d2f
scatter_alpha=.2

# ╔═╡ cab3f6a6-f982-459f-9fa0-8cc7290cd29e
node_filt = [k for k in 1:n_nodes if k!= node]

# ╔═╡ ec0278e1-0300-4842-80f0-857b414867a9
μ = 1

# ╔═╡ ad960911-9354-4a31-bcbc-f1f0a95482ae
md"""
### Analysis of regression-based mefs estimation
"""

# ╔═╡ dede3a80-c551-4673-97d4-66640ac0023b
md"""
## Same analysis, but with two groups
"""

# ╔═╡ fd198417-f955-443f-b7bd-40b5106356eb
hour_filt1 = 10 .<= hours[2:end] .<= 20;

# ╔═╡ 7991857c-eef0-4536-89df-600131e188d6
function get_reg_mefs(X, y, μ)
	return inv(X'*X + μ*I(size(X, 2)))*X'*y
end

# ╔═╡ 66ecb9fd-b976-432a-9df3-7a6a8618bc2e
# linear regression on all nodes at once, over all timesteps lumped together
begin
	X = transpose(Δd); # matrix has many nodes with zero demand; therefore columns have to be filtered so as to use only the nodes with actual mef
	idx_valid = vec(sum(X .=== missing, dims=2).== 0)
	X_demand = X[idx_valid, :]
	demand_nodes = findall(vec(sum(abs.(X_demand), dims=1)).!==0.0)
	X_demand = X_demand[:, demand_nodes]
	rank_ = rank(X_demand)
	mef_opt = get_reg_mefs(X_demand, ΔE[idx_valid], μ)
	mefs_reg = zeros(n_nodes)
	mefs_reg[demand_nodes] .= mef_opt
end;

# ╔═╡ 2e59c950-93e7-48d4-840c-c40f7296614d
begin
	# compute R2 for this regression
	preds = X_demand*mef_opt
	R2 = get_R2(preds, ΔE[idx_valid])
end;

# ╔═╡ 05305b47-b71a-4452-96ca-ef809924d5c2
md"""
The regression has predictive power: R^2 = $(round(R2 * 100)) %
"""

# ╔═╡ 4afd6719-cff9-4891-87a1-82405b07c084
let
	plot_preds(preds, ΔE[idx_valid])
end

# ╔═╡ ded35429-73ce-461c-9144-84741cabb512
res_pred = transpose(Δd)[:, node_filt] * mefs_reg[node_filt];

# ╔═╡ 07c3faea-d89b-4e88-98f6-21423aa8d202
let
	lw = 4
	p = scatter(Δd[node, :], ΔE-res_pred, label=false, alpha=scatter_alpha)
	xx = LinRange(minimum(skipmissing(Δd[node, :])), maximum(skipmissing(Δd[node, :])), 10)
	plot!(xx, mefs_reg[node]*xx, lw=lw, label="Reg MEF=$(round(mefs_reg[node]))")
	plot!(xx, mean(mefs[node, :])*xx, ls=:dash, lw=lw, label="Mean True MEF=$(round( mean(mefs[node, :])))")
	xlabel!("Δd at node $(node)")
	ylabel!(L"$\Delta E -\sum_{j\neq i} \Delta d_j \beta_j$")
	title!("ΔE total wrt Δd at node $(node)")
	p
end

# ╔═╡ 3dc3b990-1e8b-43b7-856a-17667c6be634
# histogram
let
	bins = -1000:100:3000
	histogram(mefs[node, :], bins=bins, alpha=0.4, label="True")
	histogram!(mefs_obs[node, :], bins=bins, alpha=0.4, label="Observed")
	plot!([mefs_reg[node]], seriestype="vline", label="Regression", lw=3)
	plot!([median(mefs[node, :])], seriestype="vline", label="Median True", lw=3)
	plot!([median(mefs_obs[node, :])], seriestype="vline", label="Median Observed", lw=3, color="firebrick")
	xlabel!("MEF")
	ylabel!("Counts")
end

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

# ╔═╡ e3c716f3-0460-4b6b-a49a-f783a23950b8
# linear regression on all nodes at once, over all timesteps lumped together
begin
	X_demand_g1 = X_demand[hour_filt1[idx_valid], :]
	X_demand_g2 = X_demand[.!hour_filt1[idx_valid], :]
	mefs_g1, mefs_g2 = zeros(n_nodes), zeros(n_nodes)
	mefs_g1[demand_nodes] = get_reg_mefs(X_demand_g1, ΔE[idx_valid][hour_filt1[idx_valid]], μ)
	mefs_g2[demand_nodes] = get_reg_mefs(X_demand_g2, ΔE[idx_valid][.!hour_filt1[idx_valid]], μ)
end;

# ╔═╡ 581b23b6-4850-4617-97b9-9775c2cad229
begin
	# get R2 for both subgroups
	preds_g1 = X_demand_g1 * mefs_g1[demand_nodes]
	R2_g1 = get_R2(preds_g1, ΔE[idx_valid][hour_filt1[idx_valid]])

	preds_g2 = X_demand_g2 * mefs_g2[demand_nodes]
	R2_g2 = get_R2(preds_g2, ΔE[idx_valid][.!hour_filt1[idx_valid]])

end;

# ╔═╡ 2e65831e-97e0-4b13-8e00-a8bc7741e026
md"""
The prediction power for both subgroups are: 
- Group 1: R2 = $(round(R2_g1*100)) %
- Group 2: R2 = $(round(R2_g2*100)) %
"""

# ╔═╡ 9d314690-0f1d-46c8-92e8-72de962a580a
# visualizing the prediction power
let
p1 = plot_preds(preds_g1, ΔE[idx_valid][hour_filt1[idx_valid]])
p2 = plot_preds(preds_g2, ΔE[idx_valid][.!hour_filt1[idx_valid]])
plot(p1, p2)
end

# ╔═╡ 7d4a7daf-778e-4040-997a-be1cc95c75c0
# histogram
let
	bins = -1000:100:3000
	histogram(mefs[node, :], bins=bins, alpha=0.4, label="True")
	# histogram!(mefs_obs[node, :], bins=bins, alpha=0.4, label="Observed")
	plot!([mefs_g1[node]], seriestype="vline", label="Regression 1", lw=3)
	plot!([mefs_g2[node]], seriestype="vline", label="Regression 2", lw=3)
	# plot!([median(mefs[node, :])], seriestype="vline", label="Median True", lw=3)
	# plot!([median(mefs_obs[node, :])], seriestype="vline", label="Median Observed", lw=3, color="firebrick")
	xlabel!("MEF")
	ylabel!("Counts")
end

# ╔═╡ 621c8dc3-fd4f-48da-9d49-cdad33941975
md"""
# Extracting regression mef from a subset of node
"""

# ╔═╡ 5b287f18-ab7a-4588-b231-b1b6f1315ddd
th_cor = .99

# ╔═╡ a46a11ed-7433-4b36-baca-11f899755eca
conn_nodes, XX = analysis.get_connected_clusters(X_demand, th_cor);

# ╔═╡ bef12ca3-b2f1-4a41-a343-0ac98b9155e7
let
	heatmap(XX.>=th_cor)
	xlabel!("Node id")
	ylabel!("Node id")
end

# ╔═╡ 45335ae0-b302-4dfb-966f-b6cc0adb5dc7
md"""
we see we have $(length(conn_nodes)) clusters of connected nodes. We are expecting that these nodes have demands that vary nearly equally. 

For instance, let us visualize the demand patterns of a given group
"""

# ╔═╡ 4b34e966-42b8-4d87-9e79-da678ab51b29
plot([length(g) for g in conn_nodes])

# ╔═╡ 7af9c674-ad8b-48e7-b228-42814129ece5
md"""
We see that for all groups, all curves actually collapse on one another. 
"""

# ╔═╡ c0869e4f-5960-40f7-8cf2-043fb7a95a84
md"""
Now, we can summarize each group by a single demand vector, the total demand across this group!
"""

# ╔═╡ 4042ad0a-0fe2-4c50-a5ef-692b24c9eae3
# contains total Δd for each group
X_demand_gp = hcat([sum(X_demand[:, nodes_gp], dims=2) for nodes_gp in conn_nodes]...);

# ╔═╡ 265a5b81-f2a2-43d2-b1d0-007c6c39d104
# linear regression on all nodes at once, over all timesteps lumped together
begin

	mef_gp = get_reg_mefs(X_demand_gp, ΔE[idx_valid], μ)
	# mefs_reg = zeros(n_nodes)
	# mefs_reg[demand_nodes] .= mef_opt
end;

# ╔═╡ dbbc2ca0-9121-40aa-8394-2af3b84bacc0
md"""
The MEF should represent, on average over this cluster, how much emissions change wrt a change in demand. Therefore it should be equal to the average mef in the group? (Maybe weighted average)?
"""

# ╔═╡ db3cf5fa-c088-40ab-bec0-0417ee4c9ba6
@bind gp_id Slider(1:length(conn_nodes))

# ╔═╡ a4179927-2b3a-46bb-9ca1-327dbbc3b44c
idx_valid_demand = vec(sum(demand.===missing, dims=1) .==0);

# ╔═╡ 0a1f8447-6c8b-45ff-875f-74748980d5f4
d_gp = demand[demand_nodes, idx_valid_demand][conn_nodes[gp_id], :];

# ╔═╡ ecec3381-2d36-4e2c-bb07-ee24681eb391
begin
	p = plot()
	for k in 1:size(d_gp, 1)
		plot!(d_gp[k, :]/d_gp[k, 1], alpha=.4, ls=:dash, legend=false)
	end
	xlims!(4000, 4100)
	title!("Group $(gp_id), with $(size(d_gp,1)) nodes.")
	p
end

# ╔═╡ 39c546d6-c267-4b45-86fc-79426ae0d39e
# factor in the weighted average
ω = mean(d_gp./sum(d_gp, dims=1), dims=2)

# ╔═╡ 624f9b72-7792-4b6f-b138-669751a749b3
# weighted average mefs over the region
weighted_mef_gp = ω'*mefs[conn_nodes[gp_id], :]

# ╔═╡ 515c9225-87b4-424e-845c-cd6da662e25a
let
	histogram(vec(weighted_mef_gp), label="Weighted avg MEF")
	plot!([mef_gp[gp_id]], seriestype="vline", label="Regression group", lw=3)
	plot!([median(mefs_reg[demand_nodes][conn_nodes[gp_id]])], seriestype="vline", label="Reg node, median in group", lw=3)
	xlabel!("MEF")
	ylabel!("Counts")
end

# ╔═╡ 8d333cbb-35ba-409e-aa0c-9791ceec722f
let
	histogram(mef_gp, alpha=.4)
	histogram!(mefs_reg[demand_nodes], alpha=.4)
end

# ╔═╡ dc7ebc3f-38bb-4416-a488-7035ad34d1e8
md"""
## Orders of magnitude of estimated MEFS with these methods
"""

# ╔═╡ 6166cad7-1b60-4df2-95d3-1974d297dce0
let
	bins = -1000:100:3000
	alpha=0.4
	histogram(mean(mefs, dims=2), bins=bins, alpha=alpha, label="True")
	histogram!(mefs_reg, bins=bins, alpha=alpha, label="Reg")
	xlabel!("MEF")
	ylabel!("Counts")
	title!("Average values over all timesteps")
end

# ╔═╡ a757ff7e-402d-4ba3-9ecd-f2e1cd539a6e
hours_filt_day = 10 .< (1:24) .< 20

# ╔═╡ 7ee50896-add5-4e32-a72c-5470ce6c7849
true_mefs_g1 = [mean(mean(mefs_[hours_filt_day, :], dims=2)) for mefs_ in mefs_nodes_days]

# ╔═╡ 984d0af5-0afc-4146-b570-0d6a0547f3b8
true_mefs_g2 = [mean(mean(mefs_[.!hours_filt_day, :], dims=2)) for mefs_ in mefs_nodes_days]

# ╔═╡ 09b71d83-53ad-4b7b-ba70-b51ea90aa37e
let
	c1 = :skyblue2 
	c2 = :firebrick
	lw = 4
	
	scatter(Δd[node, hour_filt1], (ΔE-res_pred)[hour_filt1], alpha=scatter_alpha, label="10-20", color=c1)
	scatter!(Δd[node, .!hour_filt1], (ΔE-res_pred)[.!hour_filt1], alpha=scatter_alpha, label="20-10", color=c2)
	xx = LinRange(minimum(skipmissing(Δd[node, :])), maximum(skipmissing(Δd[node, :])), 10)
	plot!(xx, mefs_g1[node]*xx, lw=4, label="MEF=$(round(mefs_g1[node]))", color=c1)
	plot!(xx, mefs_g2[node]*xx, lw=4, label="MEF=$(round(mefs_g2[node]))", color=c2)
	plot!(xx, true_mefs_g1[node]*xx, lw=4, label="MEF=$(round(mefs_g1[node]))", color=c1, ls=:dash)
	plot!(xx, true_mefs_g2[node]*xx, lw=4, ls=:dash, label="MEF=$(round(mefs_g2[node]))", color=c2)
	
	xlabel!("Δd at node $(node)")
	ylabel!("ΔE")
	title!("ΔE total wrt Δd at node $(node)")
end

# ╔═╡ f5f2c2be-b9ff-4b3f-bf3e-65ece46931fa
let
	bins = -1000:100:3000
	alpha=0.4
	histogram(
		true_mefs_g1, 
		bins=bins, alpha=alpha, label="True")
	histogram!(mefs_g1, bins=bins, alpha=alpha, label="Reg")
	xlabel!("MEF")
	ylabel!("Counts")
	title!("Average values over Group 1")
end

# ╔═╡ efaee005-a8f5-48fd-9db6-1169ba2e1aa0
let
	bins = -1000:100:3000
	alpha=0.4
	histogram(
		true_mefs_g2, 
		bins=bins, alpha=alpha, label="True")
	histogram!(mefs_g2, bins=bins, alpha=alpha, label="Reg")
	xlabel!("MEF")
	ylabel!("Counts")
	title!("Average values over Group 2")
end

# ╔═╡ 0defd4c1-dbb9-4b66-bca0-3c4605c252bb
md"""
### Compute regression-based mefs at a given hour
"""

# ╔═╡ e8c4a66c-94b7-4b0e-9b43-fcaa0f235733
@bind hr Slider(1:24)

# ╔═╡ 6b5f3833-6168-4d45-9283-0fec34404658
hr

# ╔═╡ f8ebd537-bacc-4e1a-af8f-12f1196536ff
ΔE_hr = [ΔE[hr + k*24] for k in 0:length(results[idx])-1];

# ╔═╡ 827932a6-68b6-4b37-8f02-fdfa819fee2c
Δd_hr = hcat([Δd[:, hr + k*24] for k in 0:length(results[idx])-1]...);

# ╔═╡ 1790b9e8-dc3c-4e47-8d8e-0a0cdea118fd
begin
	scatter(Δd_hr[node, :], ΔE_hr, legend=false)
	ylabel!("ΔE at given hour")
	xlabel!("Δx at given hour")
end

# ╔═╡ 665867f1-3adc-4de8-8fd1-bfafa3a56ae9
[median(mef_reg[demand_nodes][gp_id])]

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
# ╠═8876dc66-7d6d-48ac-94d5-5dd489e81b87
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
# ╠═11aef0e7-372d-4b3b-845d-e63bc620b5f9
# ╠═f2197b5f-562a-435a-88a4-fc26415e77a9
# ╠═6b3a9e01-e14c-488d-bde5-7eae4e8e87e0
# ╟─ae72be0b-b54f-47bd-9888-ad297463bb38
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
# ╟─8d206091-89b8-49c3-af53-95a262298426
# ╟─7bbd7076-5a16-4914-8118-0d9238f75569
# ╟─e8fee23c-dc1a-430a-b896-31f09c4c503e
# ╟─1bed3bbb-7411-441c-8028-f4dc91b0f1aa
# ╠═a1da5ecd-219d-4e4e-b39d-1addddc89f13
# ╠═66caf393-88ce-4e6c-9a4a-afdfa7b24304
# ╠═576e6a18-eb15-4578-9bda-569c1397a8d3
# ╠═112eb399-f273-42a0-a6b6-52ce543bc696
# ╠═9b4e015d-a0b6-4046-b9c7-9f12b93d84fa
# ╠═98e1f152-67c7-43d2-b3b5-68bbc322478d
# ╟─a003c654-7072-4e3f-9aa8-154c0efa565c
# ╠═fb6a6791-c704-43f1-9784-1cd17dfa7858
# ╠═66ecb9fd-b976-432a-9df3-7a6a8618bc2e
# ╠═9323dbf3-6a22-476a-8186-766fe1515c18
# ╠═2e59c950-93e7-48d4-840c-c40f7296614d
# ╟─05305b47-b71a-4452-96ca-ef809924d5c2
# ╠═26d698fe-c2da-4a34-b446-4dd4dfacf861
# ╟─4afd6719-cff9-4891-87a1-82405b07c084
# ╠═01805ad9-0d4c-4d8e-8f79-69da710a3d2f
# ╠═cab3f6a6-f982-459f-9fa0-8cc7290cd29e
# ╠═ded35429-73ce-461c-9144-84741cabb512
# ╟─07c3faea-d89b-4e88-98f6-21423aa8d202
# ╠═ec0278e1-0300-4842-80f0-857b414867a9
# ╠═3dc3b990-1e8b-43b7-856a-17667c6be634
# ╟─ad960911-9354-4a31-bcbc-f1f0a95482ae
# ╠═7833642f-5b87-4d1d-b311-72901851a58f
# ╟─dede3a80-c551-4673-97d4-66640ac0023b
# ╠═fd198417-f955-443f-b7bd-40b5106356eb
# ╠═e3c716f3-0460-4b6b-a49a-f783a23950b8
# ╠═581b23b6-4850-4617-97b9-9775c2cad229
# ╟─2e65831e-97e0-4b13-8e00-a8bc7741e026
# ╟─9d314690-0f1d-46c8-92e8-72de962a580a
# ╟─09b71d83-53ad-4b7b-ba70-b51ea90aa37e
# ╠═7991857c-eef0-4536-89df-600131e188d6
# ╠═7d4a7daf-778e-4040-997a-be1cc95c75c0
# ╟─621c8dc3-fd4f-48da-9d49-cdad33941975
# ╠═5b287f18-ab7a-4588-b231-b1b6f1315ddd
# ╠═a46a11ed-7433-4b36-baca-11f899755eca
# ╠═bef12ca3-b2f1-4a41-a343-0ac98b9155e7
# ╟─45335ae0-b302-4dfb-966f-b6cc0adb5dc7
# ╠═4b34e966-42b8-4d87-9e79-da678ab51b29
# ╠═0a1f8447-6c8b-45ff-875f-74748980d5f4
# ╟─ecec3381-2d36-4e2c-bb07-ee24681eb391
# ╟─7af9c674-ad8b-48e7-b228-42814129ece5
# ╟─c0869e4f-5960-40f7-8cf2-043fb7a95a84
# ╠═4042ad0a-0fe2-4c50-a5ef-692b24c9eae3
# ╠═265a5b81-f2a2-43d2-b1d0-007c6c39d104
# ╟─dbbc2ca0-9121-40aa-8394-2af3b84bacc0
# ╠═db3cf5fa-c088-40ab-bec0-0417ee4c9ba6
# ╠═a4179927-2b3a-46bb-9ca1-327dbbc3b44c
# ╠═39c546d6-c267-4b45-86fc-79426ae0d39e
# ╠═624f9b72-7792-4b6f-b138-669751a749b3
# ╠═515c9225-87b4-424e-845c-cd6da662e25a
# ╠═8d333cbb-35ba-409e-aa0c-9791ceec722f
# ╟─dc7ebc3f-38bb-4416-a488-7035ad34d1e8
# ╟─6166cad7-1b60-4df2-95d3-1974d297dce0
# ╠═a757ff7e-402d-4ba3-9ecd-f2e1cd539a6e
# ╟─7ee50896-add5-4e32-a72c-5470ce6c7849
# ╟─984d0af5-0afc-4146-b570-0d6a0547f3b8
# ╟─f5f2c2be-b9ff-4b3f-bf3e-65ece46931fa
# ╟─efaee005-a8f5-48fd-9db6-1169ba2e1aa0
# ╟─0defd4c1-dbb9-4b66-bca0-3c4605c252bb
# ╠═620e86a8-60ab-4aa3-b313-29ab83ef5f4e
# ╟─e8c4a66c-94b7-4b0e-9b43-fcaa0f235733
# ╟─6b5f3833-6168-4d45-9283-0fec34404658
# ╠═f8ebd537-bacc-4e1a-af8f-12f1196536ff
# ╠═827932a6-68b6-4b37-8f02-fdfa819fee2c
# ╠═1790b9e8-dc3c-4e47-8d8e-0a0cdea118fd
# ╠═665867f1-3adc-4de8-8fd1-bfafa3a56ae9
