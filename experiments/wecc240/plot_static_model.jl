### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 0ae79723-a5cf-4508-b41d-9622948185a9
using Pkg; Pkg.activate("")

# ╔═╡ 79ba1c40-99ab-11ec-161e-836401b0041f
using BSON

# ╔═╡ 5303b439-2bbb-4a04-b17e-7df6f2983493
using CSV, DataFrames

# ╔═╡ 668445dc-2437-421f-9251-b4044e5849f6
# We need all of these (including the weird ones like PooledArrays) 
# to make the BSON file load
using SparseArrays, InlineStrings, PooledArrays, MathOptInterface

# ╔═╡ 32e5f26a-9b2f-4fc0-a0cd-1a5f101f0db9
using StatsBase: mean

# ╔═╡ 7a42f00e-193c-45ea-951f-dcd4e1c1975f
using CairoMakie

# ╔═╡ 0f7e3ce7-9cf2-46ea-926b-43b7601246f7
using StringDistances

# ╔═╡ 2d3cf797-4cc2-4aad-bc3e-94f5474e99f9
begin
	using GeoMakie, Proj4
	using GeoMakie.GeoInterface
	using GeoMakie.GeoJSON
	using Downloads
end

# ╔═╡ 5cb1709a-eda0-41b3-8bff-f58c19608be5
using PlutoUI

# ╔═╡ 6db70f24-e8ba-461e-8d86-00e9a37b44d3
md"""
## Load data
"""

# ╔═╡ 6de86962-a420-4885-ae7a-18748549c4c2
path = "/Users/degleris/Data/carbon_networks/wecc240_static_results.bson"

# ╔═╡ 2757231c-ef30-417a-87dd-7d155049ba47
data = BSON.load(path, @__MODULE__);

# ╔═╡ b41420e8-1710-415c-b7d3-8042a19f660f
results = data[:results];

# ╔═╡ e3288d1f-4b66-49d5-9270-57827151a361
nodes = data[:meta].node_names

# ╔═╡ b4f91614-ada2-4961-8913-96855f7ca81b
md"""
## Load substation data
"""

# ╔═╡ 34b6c866-0fc6-4f7d-91b5-15e27277ce9d
wecc_states = uppercase.(["ca", "or", "wa", "nv", "mt", "id", "wy", "ut", "co", "az", "nm"])

# ╔═╡ b4f18908-f62f-479e-b808-4847c03dfd5d
substations = let
	df = DataFrame(CSV.File("/Users/degleris/Downloads/substations.csv"))
	filter!(r -> !ismissing(r.STATE) && r.STATE in wecc_states, df)
end

# ╔═╡ ec338b2e-4106-454d-8687-237602636cf1
k = 2

# ╔═╡ d987edd6-421d-43f9-b54d-d63429db2a21
nodes[k]

# ╔═╡ e4b7ffa3-9ee9-4e1b-8ec7-cf9965659c0c
findnearest(nodes[k][1:end-7], coalesce.(substations.NAME, ""), DamerauLevenshtein())

# ╔═╡ a0591816-88a4-4144-9d43-09f1205614af
function get_coords(node)
	 _, i = findnearest(node[1:end-7], coalesce.(substations.NAME, ""), DamerauLevenshtein())
	 return substations.LATITUDE[i], substations.LONGITUDE[i]
end

# ╔═╡ 1190a679-bd6c-434f-962c-7e8a26b06213
x = [c[2] for c in get_coords.(nodes)]

# ╔═╡ e6e6a795-13f2-4cc0-97ba-83de95011696
y = [c[1] for c in get_coords.(nodes)]

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

# ╔═╡ cbc71e2e-0bd1-441c-bf17-c60053a60795
md"""
## Plot!
"""

# ╔═╡ 5e79de08-58a1-4ddc-ac01-5a8f48a1b02d
# let
# 	fig = Figure()
# 	ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false)
# 	sct = scatter!(ax, x, y, markersize=7, color=λ, colormap=:jet1, colorrange=(0.0, 1000.0))

# 	A = data[:case].A
# 	for j in 1:size(A, 2)
# 		fr = findfirst(==(-1), A[:, j])
# 		to = findfirst(==(1), A[:, j])

# 		lines!(ax, [x[fr]; x[to]], [y[fr]; y[to]], color=(:gray, 0.1))
# 	end
# 	Colorbar(fig[1, 2], sct, label="Marginal Emissions Rate [kg CO2 / MWh]")

# 	ax.xticks = [-150]
# 	ax.yticks = [20]
# 	ax.title = "WECC 2004: Average Nodal MEFs at Hour $hour"

# 	fig
# end;
"Old plot"

# ╔═╡ 59316c15-a94c-4c56-a30a-0e6c23629de7
hour = 6

# ╔═╡ 161fcf67-b65b-4661-bc9e-ff714268b444
λ = mean(skipmissing(mefs[hour, :]))[:, 1]

# ╔═╡ 2c0c2056-01a0-48eb-852d-92a172703975
all_λs = reduce(vcat, skipmissing(mefs[hour, :]))[:, 1]

# ╔═╡ 7ffbe1bc-8cc6-4033-a70b-880209164199
let
	# Everthing in === is from https://lazarusa.github.io/BeautifulMakie/GeoPlots/geoCoastlinesStatesUS/


	# ===========
	states_url = "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/us-states.json"
    states = Downloads.download(states_url)
    states_geo = GeoJSON.read(read(states, String))
    n = length(GeoInterface.features(states_geo))

    trans = Proj4.Transformation("+proj=longlat +datum=WGS84", "+proj=latlong", 
        always_xy=true) 
    # see https://proj.org/operations/projections/index.html for more options 
    ptrans = Makie.PointTrans{2}(trans)

    fig = Figure(resolution = (1250,700), fontsize = 22)
    ax = Axis(fig[1,1], aspect = DataAspect(), 
        title = "Projection: Winkel Tripel, US States")
    # all input data coordinates are projected using this function
    ax.scene.transformation.transform_func[] = ptrans
	
    xlims!(ax, -125, -100)
	ylims!(ax, 32, 51)
    
	# now the plot 
    lines!(ax, GeoMakie.coastlines(), color = :black)
    poly!(ax, states_geo, color=(:lightgray, 0.2), 
		colormap = :plasma, strokecolor = :black, 
        strokewidth = 2, overdraw = true)
	# ============

	
	# and my part
	
	# Edges
	A = data[:case].A
	for j in 1:size(A, 2)
		fr = findfirst(==(-1), A[:, j])
		to = findfirst(==(1), A[:, j])

		lines!(ax, [x[fr]; x[to]], [y[fr]; y[to]], 
			color=(:black, 0.1))
	end

	# Nodes
	sct = scatter!(ax, x, y, markersize=10, color=λ, 
		colormap=:jet1, colorrange=(0.0, 1000.0))

	# Colorbar
	Colorbar(fig[1, 2], sct, label="Marginal Emissions Rate [kg CO2 / MWh]")

	ax.xticks = [-150]
	ax.yticks = [20]
	ax.title = "WECC 2004: Average Nodal MEFs at Hour $hour"

	
    fig
end

# ╔═╡ b53cc8dd-c36e-4cf8-9f1d-473a0e985234
let
	fig, ax = hist(all_λs, bins=500, normalization=:probability)
	xlims!(ax, -100, 1500)
	ax.xlabel = "MEF"
	ax.ylabel = "Probability"

	fig
end

# ╔═╡ Cell order:
# ╠═0ae79723-a5cf-4508-b41d-9622948185a9
# ╠═79ba1c40-99ab-11ec-161e-836401b0041f
# ╠═5303b439-2bbb-4a04-b17e-7df6f2983493
# ╠═668445dc-2437-421f-9251-b4044e5849f6
# ╠═32e5f26a-9b2f-4fc0-a0cd-1a5f101f0db9
# ╠═7a42f00e-193c-45ea-951f-dcd4e1c1975f
# ╟─6db70f24-e8ba-461e-8d86-00e9a37b44d3
# ╠═6de86962-a420-4885-ae7a-18748549c4c2
# ╠═2757231c-ef30-417a-87dd-7d155049ba47
# ╠═b41420e8-1710-415c-b7d3-8042a19f660f
# ╠═e3288d1f-4b66-49d5-9270-57827151a361
# ╟─b4f91614-ada2-4961-8913-96855f7ca81b
# ╠═34b6c866-0fc6-4f7d-91b5-15e27277ce9d
# ╠═b4f18908-f62f-479e-b808-4847c03dfd5d
# ╠═0f7e3ce7-9cf2-46ea-926b-43b7601246f7
# ╠═ec338b2e-4106-454d-8687-237602636cf1
# ╠═d987edd6-421d-43f9-b54d-d63429db2a21
# ╠═e4b7ffa3-9ee9-4e1b-8ec7-cf9965659c0c
# ╠═a0591816-88a4-4144-9d43-09f1205614af
# ╠═1190a679-bd6c-434f-962c-7e8a26b06213
# ╠═e6e6a795-13f2-4cc0-97ba-83de95011696
# ╟─d2bacf4a-af37-4ff9-bebb-3dc3d06edd8a
# ╠═3c9f279c-ca03-4a8d-b9bc-3ad3668b76f7
# ╠═29d3d70d-03d2-4883-9efd-ff382afef4af
# ╠═cfcc5416-2038-48c4-a2b0-bcd92b574441
# ╠═161fcf67-b65b-4661-bc9e-ff714268b444
# ╠═2c0c2056-01a0-48eb-852d-92a172703975
# ╟─cbc71e2e-0bd1-441c-bf17-c60053a60795
# ╠═2d3cf797-4cc2-4aad-bc3e-94f5474e99f9
# ╟─5e79de08-58a1-4ddc-ac01-5a8f48a1b02d
# ╠═5cb1709a-eda0-41b3-8bff-f58c19608be5
# ╠═59316c15-a94c-4c56-a30a-0e6c23629de7
# ╟─7ffbe1bc-8cc6-4033-a70b-880209164199
# ╟─b53cc8dd-c36e-4cf8-9f1d-473a0e985234
