### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

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

# ╔═╡ 0f7e3ce7-9cf2-46ea-926b-43b7601246f7
using StringDistances

# ╔═╡ 7a42f00e-193c-45ea-951f-dcd4e1c1975f
using CairoMakie

# ╔═╡ 5cb1709a-eda0-41b3-8bff-f58c19608be5
using PlutoUI

# ╔═╡ 2d3cf797-4cc2-4aad-bc3e-94f5474e99f9
begin
	using GeoMakie
	using GeoMakie.GeoInterface
	using GeoMakie.GeoJSON
	using Downloads
end

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

# ╔═╡ ae02b617-f2d0-4fa6-86f9-3a6e4088a803
begin
	fnm_static = "wecc240_static_results_initialCode.bson"
	fnm_dynamic = "wecc240_dynamic_results_COND50_initialCode.bson"
end

# ╔═╡ 6de86962-a420-4885-ae7a-18748549c4c2
path = joinpath(DATA_DIR, fnm_static)

# ╔═╡ 2757231c-ef30-417a-87dd-7d155049ba47
data = BSON.load(path, @__MODULE__);

# ╔═╡ 37b3f4ba-9fb0-4285-aa33-f9905414c764
results = data[:results];

# ╔═╡ 8cab03dd-f034-443f-9a60-32aa87d1fde5
path_dyn = joinpath(DATA_DIR, fnm_dynamic)

# ╔═╡ 83e2c12a-fe36-4123-baf3-0e8c1bdecead
data_dyn = BSON.load(path_dyn, @__MODULE__);

# ╔═╡ 704175e9-921e-42ea-877f-35ce07610b8a
results_dyn = data_dyn[:results];

# ╔═╡ e3288d1f-4b66-49d5-9270-57827151a361
nodes = data[:meta].node_names

# ╔═╡ b4f91614-ada2-4961-8913-96855f7ca81b
md"""
## Load substation data
"""

# ╔═╡ 0f0b8da3-30b6-4166-aab8-322ee320e971
"""
Assign a geographical coordinate to nodes in the network. 
"""

# ╔═╡ bfb33d19-ceca-4a3a-87d0-d75960cd4544
begin
wecc_states = uppercase.(["ca", "or", "wa", "nv", "mt", "id", "wy", "ut", "co", "az", "nm"])
substations = let
df = DataFrame(CSV.File(joinpath(DATA_DIR, "substations.csv")))
filter!(r -> !ismissing(r.STATE) && r.STATE in wecc_states, df)
end;
end;

# ╔═╡ a0591816-88a4-4144-9d43-09f1205614af
function get_coords(node)
	 _, i = findnearest(
		 rstrip(node[1:end-10]), 
		 coalesce.(substations.NAME, ""), 
		 DamerauLevenshtein()
	 )
	 return substations.LATITUDE[i], substations.LONGITUDE[i]
end;

# ╔═╡ 5392f525-ecb3-47c7-a32f-73a6b02967df
begin
x = [c[2] for c in get_coords.(nodes)]
y = [c[1] for c in get_coords.(nodes)]
end;

# ╔═╡ d2bacf4a-af37-4ff9-bebb-3dc3d06edd8a
md"""
## Average MEFs
"""

# ╔═╡ 3c9f279c-ca03-4a8d-b9bc-3ad3668b76f7
is_valid = [d.status for d in results] .== MathOptInterface.OPTIMAL;

# ╔═╡ 29d3d70d-03d2-4883-9efd-ff382afef4af
mefs = [v ? d.λ : missing for (d, v) in zip(results, is_valid)];

# ╔═╡ cfcc5416-2038-48c4-a2b0-bcd92b574441
demands = [v ? d.d : missing for (d, v) in zip(results, is_valid)];

# ╔═╡ ad687e9a-9d7b-4990-9772-d9cfafe26421
begin
	node1 = 195
	node2 = 50
end

# ╔═╡ cbc71e2e-0bd1-441c-bf17-c60053a60795
md"""
## Plot!
"""

# ╔═╡ 59316c15-a94c-4c56-a30a-0e6c23629de7
hour = 6

# ╔═╡ 161fcf67-b65b-4661-bc9e-ff714268b444
λ = mean(skipmissing(mefs[hour, :]))[:, 1]

# ╔═╡ 2c0c2056-01a0-48eb-852d-92a172703975
all_λs = reduce(vcat, skipmissing(mefs[hour, :]))[:, 1]

# ╔═╡ a9362a09-277b-4e97-9e66-d6df99f18a70
all_mefs_a = [m[node1] for m in skipmissing(mefs[hour, :])]

# ╔═╡ fe91b3ba-3159-48b0-a3d5-7af7bfe6fc34
all_mefs_b = [m[node2] for m in skipmissing(mefs[hour, :])]

# ╔═╡ 7ffbe1bc-8cc6-4033-a70b-880209164199
function fig_map(fig = Figure(resolution=(450, 300), fontsize=10))
	# Everthing in === is from https://lazarusa.github.io/BeautifulMakie/GeoPlots/geoCoastlinesStatesUS/
	# ===========
    ax = Axis(fig[1,1])
	
	states_url = "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/us-states.json"
    states = Downloads.download(states_url)
    states_geo = GeoJSON.read(read(states, String))
    n = length(GeoInterface.features(states_geo))

    trans = Proj4.Transformation("+proj=longlat +datum=WGS84", "+proj=latlong", 
        always_xy=true) 
	
    # see https://proj.org/operations/projections/index.html for more options 
    # all input data coordinates are projected using this function
    ax.scene.transformation.transform_func[] = Makie.PointTrans{2}(trans)
	
    xlims!(ax, -125, -100)
	ylims!(ax, 32, 51)
    
	# now the plot 
    lines!(ax, GeoMakie.coastlines(), color = :black)
    poly!(ax, states_geo, color=(:lightgray, 0.2), 
		colormap = :plasma, strokecolor = :black, 
        strokewidth = 1, overdraw = true)
	# ============

	
	# and my part
	
	# Edges
	A = data[:case].A
	for j in 1:size(A, 2)
		fr = findfirst(==(-1), A[:, j])
		to = findfirst(==(1), A[:, j])

		lines!(ax, [x[fr]; x[to]], [y[fr]; y[to]], 
			color=(:black, 0.08))
	end

	# Nodes
	sct = scatter!(ax, x, y, markersize=8, marker=:hexagon, color=λ/1e3, 
		colormap=:jet1, colorrange=(-0.25, 1.5))

	# Colorbar
	cb = Colorbar(fig[1, 2], sct, label="Marginal Emissions Rate [ton CO2 / MWh]")

	ax.xticks = [-150]
	ax.yticks = [20]
	ax.title = "WECC 2004: Average Nodal MEFs at Hour $hour"

	cb.tellheight = true
	
    return fig, ax
end

# ╔═╡ 07268e37-5b62-4ab3-8d0d-5bab2286cdbe
fig_map()[1];

# ╔═╡ e5e10f07-1001-4438-b32d-c1f25cce04b1
md"""
## Analyze dynamic data
"""

# ╔═╡ b90eb7df-a78c-4bc5-ae3b-41f62e38da54
total_mef(λ) = sum(λ, dims=1)[1, :][hour]

# ╔═╡ d4d509bd-8f96-4da3-917f-a65acb569953
nodal_mef_dyn(n) = reduce(vcat, [total_mef(results_dyn[d].λ[n, :, :]) for d in 1:365])

# ╔═╡ d1f26911-bd79-4ce6-b0d8-218f8a772840
all_λs_dyn = reduce(hcat, [nodal_mef_dyn(n) for n in 1:length(nodes)]);

# ╔═╡ b53cc8dd-c36e-4cf8-9f1d-473a0e985234
function fig_distr(fig = Figure(resolution=(300, 300)))
	ax = Axis(fig[1, 1])
	hidedecorations!(ax, ticks=false, ticklabels=false, label=false)
	
	kwargs = (strokewidth=1, strokecolor=:black, direction=:y)
	density!(ax, all_mefs_a/1e3; offset=4.0, kwargs...)
	density!(ax, all_mefs_b/1e3; offset=2.0, kwargs...)
	density!(ax, all_λs/1e3; color=(:slategray, 0.7), kwargs...)

	kwargs = (direction=:y,)
	density!(ax, nodal_mef_dyn(node1)/1e3; color=(:red, 0.2), offset=4, kwargs...)
	density!(ax, nodal_mef_dyn(node2)/1e3; color=(:red, 0.2), offset=2, kwargs...)
	density!(ax, reshape(all_λs_dyn, :)/1e3; color=(:red, 0.2), kwargs...)

	kwargs = (linewidth=4, color=:black)
	hlines!(ax, [mean(all_mefs_a)/1e3]; xmin=4/6, kwargs...)
	hlines!(ax, [mean(all_mefs_b)/1e3]; xmin=2/6, xmax=4/6, kwargs...)
	hlines!(ax, [mean(all_λs)/1e3]; xmax=2/6, kwargs...)
	

	ylims!(ax, -0.25, 1.5)
	ax.ylabel = "MEF [ton CO2 / MWh]"
	
	
	xlims!(ax, 0, 6)
	ax.xlabel = "Frequency"
	ax.xticks = [0, 2, 4]
	ax.xtickformat = xs -> ["All", "LA", "SF"]

	return fig, ax
end

# ╔═╡ e1a1acda-1d52-45bd-8257-8b7249318c9b
fig_distr()[1];

# ╔═╡ c6f2eb39-a0e6-44bf-8649-f25ef72961a4
full_figure = let
	fig = Figure(resolution=(650, 300), fontsize=10)

	f1, ax1 = fig_distr(fig[2, 1])
	f2, ax2 = fig_map(fig[2, 2])


	colsize!(fig.layout, 1, Auto(0.5))
	
	ax2.title = ""
	ax1.ylabel = "Marginal Emissions Rate [ton CO2 / MWh]"

	for (label, layout) in zip(["A", "B"], [fig[2, 1], fig[2, 2]])
    	Label(layout[1, 1, TopLeft()], label,
	        textsize = 18,
			font="Noto Sans Bold",
	        padding = (0, 0, 5, 0),
	        halign = :right
		)
	end

	Label(fig[1, 1:2], "WECC 2004: Nodal MEFs at Hour $hour", 
		textsize=12,
		padding=(0, 0, 0, 0),
		valign=:bottom
	)
	rowgap!(fig.layout, 1, 6)



	fig
end

# ╔═╡ 5154fdd8-a58d-4faa-aced-7212ed0dc705
save(joinpath(RESULTS_DIR, "wecc240_full_figure.pdf"), full_figure)

# ╔═╡ dfc765e0-39d3-4ae4-93f0-4f0406f9f358
λ_dyn = mean(all_λs_dyn, dims=1)[1, :]

# ╔═╡ 305bdb8a-e186-4c5a-b927-b7a406ac260a
md"""
### Compare regression-based MEF vs exact differentiation-based MEFs
"""

# ╔═╡ 701285fc-887d-425c-98a3-1ee09235066f
md"""
Todo: 
- Compute total emissions
- Enable regression-based MEF estimation
- Compare with diagonal terms in the matrix
- 
- Compare with analytic/exact MEF computation
"""

# ╔═╡ 1d5d1c1c-050a-434e-a977-b8d2aea69b26
co2_costs = util.get_costs(data[:case].heat, data[:case].fuel, util.FUEL_EMISSIONS)

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

# ╔═╡ 13478bcb-c4fc-4532-a1ef-4513b15e3295
# changes in emissions

ΔE = [(E_tot[d][hour] - E_tot[d][hour-1]) for d in 1:365];

# ╔═╡ 70723775-1912-48ed-9ae5-d4663d0f81d3
# gather all the demand changes, at ever node

Δd = hcat([
	[results_dyn[d].d[hour][n] - results_dyn[d].d[hour-1][n] for d in 1:365] 
	for n in 1:length(nodes)
]...);

# ╔═╡ 1fba9a4a-8090-4fc5-a373-670ed04dfb4e
# get the values for the scatter plot at the two nodes of interest

begin
xx1 = Δd[:, node1]
xx2 = Δd[:, node2]
end;

# ╔═╡ 0df59933-2c53-46ae-88fe-a7fbd1b4b339
lr_E = linregress(Δd, ΔE)

# ╔═╡ 50acbd4b-1a02-4c55-b605-caf07f12bd74
MEFs_reg = LinearRegression.slope(lr_E)

# ╔═╡ 749ef6df-5909-4f40-a5ff-0ed7877de9ab
lines(MEFs_reg)

# ╔═╡ db737f8d-ed35-43c8-84fb-c7329763b3c1
MEFs_reg[node1]

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
# ╟─cab761be-ee1b-4002-a187-df72c29d9771
# ╟─113b99d8-0708-4da8-a8d6-7c60734e4a31
# ╟─f13e5dea-33b9-45c0-876e-42654fe6a8c7
# ╟─6db70f24-e8ba-461e-8d86-00e9a37b44d3
# ╠═f9fab4fe-baec-4bfd-9d84-ef9caac85f5f
# ╠═d7598abb-2be7-4e3b-af9e-14827ef5a3b0
# ╠═45c73bb3-eecf-4b92-8e91-6a4c83addfdc
# ╠═ae02b617-f2d0-4fa6-86f9-3a6e4088a803
# ╠═6de86962-a420-4885-ae7a-18748549c4c2
# ╠═2757231c-ef30-417a-87dd-7d155049ba47
# ╠═37b3f4ba-9fb0-4285-aa33-f9905414c764
# ╠═8cab03dd-f034-443f-9a60-32aa87d1fde5
# ╠═83e2c12a-fe36-4123-baf3-0e8c1bdecead
# ╠═704175e9-921e-42ea-877f-35ce07610b8a
# ╠═e3288d1f-4b66-49d5-9270-57827151a361
# ╟─b4f91614-ada2-4961-8913-96855f7ca81b
# ╟─0f0b8da3-30b6-4166-aab8-322ee320e971
# ╠═0f7e3ce7-9cf2-46ea-926b-43b7601246f7
# ╠═a0591816-88a4-4144-9d43-09f1205614af
# ╠═bfb33d19-ceca-4a3a-87d0-d75960cd4544
# ╠═5392f525-ecb3-47c7-a32f-73a6b02967df
# ╟─d2bacf4a-af37-4ff9-bebb-3dc3d06edd8a
# ╠═3c9f279c-ca03-4a8d-b9bc-3ad3668b76f7
# ╠═29d3d70d-03d2-4883-9efd-ff382afef4af
# ╠═cfcc5416-2038-48c4-a2b0-bcd92b574441
# ╠═161fcf67-b65b-4661-bc9e-ff714268b444
# ╠═2c0c2056-01a0-48eb-852d-92a172703975
# ╠═ad687e9a-9d7b-4990-9772-d9cfafe26421
# ╠═a9362a09-277b-4e97-9e66-d6df99f18a70
# ╠═fe91b3ba-3159-48b0-a3d5-7af7bfe6fc34
# ╟─cbc71e2e-0bd1-441c-bf17-c60053a60795
# ╠═7a42f00e-193c-45ea-951f-dcd4e1c1975f
# ╠═5cb1709a-eda0-41b3-8bff-f58c19608be5
# ╠═2d3cf797-4cc2-4aad-bc3e-94f5474e99f9
# ╠═59316c15-a94c-4c56-a30a-0e6c23629de7
# ╠═07268e37-5b62-4ab3-8d0d-5bab2286cdbe
# ╟─7ffbe1bc-8cc6-4033-a70b-880209164199
# ╠═e1a1acda-1d52-45bd-8257-8b7249318c9b
# ╟─b53cc8dd-c36e-4cf8-9f1d-473a0e985234
# ╟─c6f2eb39-a0e6-44bf-8649-f25ef72961a4
# ╠═5154fdd8-a58d-4faa-aced-7212ed0dc705
# ╟─e5e10f07-1001-4438-b32d-c1f25cce04b1
# ╠═b90eb7df-a78c-4bc5-ae3b-41f62e38da54
# ╠═d4d509bd-8f96-4da3-917f-a65acb569953
# ╠═d1f26911-bd79-4ce6-b0d8-218f8a772840
# ╠═dfc765e0-39d3-4ae4-93f0-4f0406f9f358
# ╟─305bdb8a-e186-4c5a-b927-b7a406ac260a
# ╟─701285fc-887d-425c-98a3-1ee09235066f
# ╠═1d5d1c1c-050a-434e-a977-b8d2aea69b26
# ╠═a34e32d8-3e18-4c22-87e2-032360661498
# ╠═fa51de9f-6f7a-4173-96d4-18f5f225aa1a
# ╟─be2f82a2-5c22-41a3-abe4-5a6f9de4e6d7
# ╠═82ad33b1-3719-4aeb-9c82-1f8f1812ed61
# ╠═1fba9a4a-8090-4fc5-a373-670ed04dfb4e
# ╠═13478bcb-c4fc-4532-a1ef-4513b15e3295
# ╠═70723775-1912-48ed-9ae5-d4663d0f81d3
# ╠═0df59933-2c53-46ae-88fe-a7fbd1b4b339
# ╠═50acbd4b-1a02-4c55-b605-caf07f12bd74
# ╠═749ef6df-5909-4f40-a5ff-0ed7877de9ab
# ╠═db737f8d-ed35-43c8-84fb-c7329763b3c1
# ╠═28591eda-995f-4224-98b8-ca1eab27559e
# ╠═8e1f0bef-459d-4598-a01a-e59e00f53247
# ╠═6f2896ac-ecf6-4256-82e0-0a4070fa14af
# ╠═982d36b3-a7b6-4066-9441-9e02596423dd
# ╟─e834a7d2-5dd9-4995-b3c9-3c771673edaa
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
