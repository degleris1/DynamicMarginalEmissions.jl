### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 85394d0a-6f17-11ec-04e7-89c834098e0c
begin
	using Pkg; Pkg.activate(joinpath(@__DIR__, "../../dev"))
	
	using HDF5
	using CairoMakie
	using ColorSchemes
	using Colors
	using TOML
end

# ╔═╡ 00c874fc-f087-4c6a-b653-feb5080ff89d
using ColorSchemeTools

# ╔═╡ 5aab77a4-ef85-4676-aeb6-eb2ce0aaeb22
using LinearAlgebra

# ╔═╡ a643ae79-d173-423a-9b00-0e8bdaded832
config = TOML.parsefile(joinpath(@__DIR__, "../../config.toml"))

# ╔═╡ ff86eb8b-7619-40cf-b414-240ba3e44c7f
config["data"]

# ╔═╡ 177d72ea-5a88-416e-99ee-4ed694b0384e
begin
	#Anthony
	DATA_DIR = config["data"]["DATA_DIR"]
	SAVE_PATH = config["data"]["SAVE_DIR"]
	
	fnm1 = "mef_illustration.pdf"
	fnm2 = "mef_heatmaps.pdf"
end;

# ╔═╡ c041f8ff-df9a-4b6a-a73b-b125c3de4495
md"""
### Load data
"""

# ╔═╡ 6986e19a-c052-4e0f-b765-46496e2c8a36
function load_h5(path)
	return h5open(path) do f
		r = Dict()

		for k in keys(f)
			r[k] = read(f[k])
		end

		r
	end
end

# ╔═╡ 9874d30f-c591-48a4-8118-a523bc758509
f = load_h5(joinpath(DATA_DIR, "fig2_data_1"));

# ╔═╡ be52012c-4a74-4b3d-bdaf-ce380c8b08b0
curve_nodes = [3, 21, 24]

# ╔═╡ 49af4bc2-be22-48b7-a08f-e5db775a64dc
curve_datasets = let
	r = Dict()

	for k in curve_nodes
		r[k] = load_h5(joinpath(DATA_DIR, "fig2_data_sensitivity_node_$k"))
	end

	r
end

# ╔═╡ 48d56784-63e1-408e-8b0b-0dc1d858b836
diff_pts = collect(95:5:105)

# ╔═╡ c33a5eb8-ae22-412c-90ef-9cb1ab385ecb
md"""
### Plot A: MEFs across the network
"""

# ╔═╡ 672984c0-eff8-4156-9ce8-0d2f919a4d4c
figure_nodal_mefs = let
	nodal_mefs = f["MEFs"] / 1000
	
	# Setup
	size_inches = (3, 2)
	size_pt = 72 .* size_inches
	fig = Figure(resolution=size_pt ./ 0.75, fontsize=10)

	# Top panels
	ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false)
	node = 21
	
	data = curve_datasets[node]
	x = data["x"]
		
	lines!(ax, x, data["y_exp"], label="True", linewidth=3)
	for d in diff_pts
		lines!(ax, x, data["y_$d"], label="Estimates", 
			linestyle=:dash, linewidth=2, color=:red)
		scatter!(ax, [d/100], [data["y_$d"][findfirst(==(d/100), x)]],
			color=(:red, 0.75), markersize=8)
	end
	
	xlims!(ax, 0.9, 1.1)
	ylims!(ax, 0.98, 1.09)
	ax.ylabel = "ΔE [%]"
	ax.xlabel = "ΔD [%]"
	#ax.title = "Node $node"
	axislegend(ax, position=:lt, padding=(4.0, 4.0, 2.0, 2.0), rowgap=0, unique=true)

	# # Setup
	# size_inches = (6.5, 4)
	# size_pt = 72 .* size_inches
	# fig = Figure(resolution=size_pt ./ 0.75, fontsize=10)

	# # Top panels
	# axes_top = [Axis(fig[1, k], xgridvisible=false, ygridvisible=false) for k in 1:3]
	# for (ind_node, node) in enumerate(curve_nodes)
	# 	ax = axes_top[ind_node]
	# 	data = curve_datasets[node]
	# 	x = data["x"]
		
	# 	lines!(ax, x, data["y_exp"], label="True", linewidth=3)
	# 	for d in diff_pts
	# 		lines!(ax, x, data["y_$d"], label="Estimates", 
	# 			linestyle=:dash, linewidth=2, color=:red)
	# 		scatter!(ax, [d/100], [data["y_$d"][findfirst(==(d/100), x)]],
	# 			color=(:red, 0.75), markersize=8)
	# 	end
		
	# 	xlims!(ax, 0.9, 1.1)
	# 	ylims!(ax, 0.98, 1.09)
	# 	ax.ylabel = "ΔE [%]"
	# 	ax.xlabel = "ΔD [%]"
	# 	ax.title = "Node $node"
	# 	axislegend(ax, position=:lt, padding=(4.0, 4.0, 2.0, 2.0), rowgap=0, unique=true)
	# end


	# # Bottom panel
	# ax_bot = Axis(fig[2, 1:3], xgridvisible=false, ygridvisible=false)
	# barplot!(ax_bot, nodal_mefs, strokewidth=0.5)
	# xlims!(ax_bot, 0, 31)
	# ax_bot.ylabel = "MEF [mTCO2 / MWh]"
	# ax_bot.xlabel = "Node"

	# # Add panel labels
	# for (label, layout) in zip(["A", "B"], [fig[1, 1], fig[2, 1]])
	#     Label(
	# 		layout[1, 1, TopLeft()], 
	# 		label,
	#         textsize=18,
	# 		font="Noto Sans Bold",
	#         padding=(0, 5, 5, 0),
	#         halign=:right
	# 	)
	# end

	fname = joinpath(SAVE_PATH, fnm1)
	save(fname, fig, pt_per_unit=0.75)
	fig
end

# ╔═╡ b6cf83e6-043e-452a-ba72-2e5911c0cffa
md"""
### Plot B: MEFs over time
"""

# ╔═╡ 274013ad-c939-467f-8dfd-cc199340bd2c
storage_pcs = ["0% Storage", "5% Storage", "10% Storage"]

# ╔═╡ 2c443306-d607-42ff-b36a-4822b61e0846
# http://juliagraphics.github.io/Colors.jl/stable/namedcolors/
# useful link to choose appropriate colors

# ╔═╡ b92163fd-7118-4a9b-b0b0-27bbca819b7f
let
	hours = f["MEF_vs_t_x"]
	mef_grids = [f["hm_23_$(k)"] for k in 1:3] * 1.5

	mef_zero_obs = diag(mef_grids[1])
	mef_zero_true = sum(mef_grids[1], dims=2)[:, 1]

	mef_stor_obs = diag(mef_grids[2])
	mef_stor_true = sum(mef_grids[2], dims=2)[:, 1]

	fig = Figure(resolution=(600, 300), fontsize=12)
	ax = Axis(fig[1, 1])
	hidedecorations!(ax, ticks=false, ticklabels=false, label=false)

	lines!(ax, mef_zero_obs, linewidth=4, color=(:black, 0.5), label="0% Storage (Observed)")
	lines!(ax, mef_stor_obs, linewidth=4, color=(:green, 0.5), label="5% Storage (Observed)")
	
	lines!(ax, mef_stor_true, linewidth=4, color=:green, linestyle=:dash, label="5% Storage - True")

	ax.xticks = [0, 6, 12, 18, 24]
	xlims!(ax, 1, 24)
	ax.xlabel = "Hour"
	
	ax.ylabel = "MEF [kg CO2 / MWh]"
	ylims!(ax, 0, 5000)



	Legend(fig[2, 1], ax)
	save("/Users/degleris/Downloads/mef_storage_C.pdf", fig, pt_per_unit=0.75)

	fig
end

# ╔═╡ 4e579820-e404-4a44-900f-cbd588dd931e
figure_temporal = let
	hours = f["MEF_vs_t_x"]
	total_mefs = [f["MEF_vs_t_y1"], f["MEF_vs_t_y2"]] / 1000
	nodes = Dict(1 => 21, 2 => 23)
	mef_grids = [f["hm_$(nodes[j])_$(k)"] for j in 1:2, k in 1:3] ./ 1000
	cmax = maximum(maximum.(mef_grids))
	cmin = minimum(minimum.(mef_grids))

	lim = max(abs(cmin), abs(cmax))
	cmax = lim
	cmin = -lim
	
	colormap = cgrad(
		[:firebrick4, :orangered2,  :grey96, :skyblue1, :navyblue], 
		[.3, .49, (0-cmin)/(cmax-cmin), .51, .8],
		rev=true
	)
	
	# Setup
	lw = 3
	size_inches = (6.5, 4)
	size_pt = 72 .* size_inches
	
	fig = Figure(resolution=size_pt ./ 0.75, fontsize=10)
	
	# Right panel
	subplts = [
		fig[j, k] = GridLayout()
		for j in 1:2, k in 1:3
	]
	for (j, k) in Iterators.product(1:2, 1:3)
		M = f["hm_$(nodes[j])_$(k)"]./1000
		m = sum(M, dims=1)[1, :]

		ax_top = Axis(subplts[j, k][1, 1], ygridvisible=false, xgridvisible=false)
		ax_main = Axis(subplts[j, k][2, 1])
		linkxaxes!(ax_top, ax_main)

		xlims!(ax_top, 1, 24)
		if j == 1
			ylims!(ax_top, -1, 1)
			ax_top.yticks = [0, 1]
		else
			ylims!(ax_top, -1, 3.5)
			ax_top.yticks = [0, 3]
		end
		
		band!(ax_top, hours, zeros(24), m, linewidth=lw, color=(m .> 0) .- 0.5, colormap=colormap, colorrange=(-1, 1))
		ax_top.xticks = [-100]
		if j == 1
			ax_top.title = storage_pcs[k]
		end

		heatmap!(
			ax_main, M, 
			colorrange=(cmin, cmax), colormap=colormap,
		)
		ax_main.xticks = 0:6:24
		ax_main.yticks = 0:6:24
		

		if k != 1
			ax_main.ytickformat = xs -> ["" for x in xs]
		end
		if j != 2
			ax_main.xtickformat = xs -> ["" for x in xs]
		end

		rowgap!(subplts[j, k], 5)
		rowsize!(subplts[j, k], 1, Auto(1))
		rowsize!(subplts[j, k], 2, Auto(4))

		if k == 1
			Label(subplts[j, k][1, 1, Left()], "Total", 
				rotation=pi/2, padding=(0, 30, 0, 0))
			Label(subplts[j, k][2, 1, Left()], "Emissions Time", 
				rotation=pi/2, padding=(0, 30, 0, 0))
		end

	end
	cb = Colorbar(
		fig[1:2, 4], 
		colormap = colormap, 
		limits=(cmin, cmax), 
		label="MEF [TCO2 / MWh]"
	)
	
	for j in 1:2
		Label(fig[j, 1, Left()], "Node $(nodes[j])", 
			rotation=pi/2, padding=(0, 50, 0, 0))
	end
	
	Label(fig[2, 1:3, Bottom()], "Consumption Time", padding=(0, 0, 0, 24))


	
	
	fname = joinpath(SAVE_PATH, fnm2)
	save(fname, fig, pt_per_unit=0.75)
	fig
end

# ╔═╡ 0b004661-a463-4f00-bad9-3532245000da
maximum(maximum.([f["hm_$(Dict(1=>21, 2=>23)[j])_$(k)"] for j in 1:2, k in 1:3]))

# ╔═╡ Cell order:
# ╠═85394d0a-6f17-11ec-04e7-89c834098e0c
# ╠═a643ae79-d173-423a-9b00-0e8bdaded832
# ╠═ff86eb8b-7619-40cf-b414-240ba3e44c7f
# ╠═177d72ea-5a88-416e-99ee-4ed694b0384e
# ╟─c041f8ff-df9a-4b6a-a73b-b125c3de4495
# ╟─6986e19a-c052-4e0f-b765-46496e2c8a36
# ╠═9874d30f-c591-48a4-8118-a523bc758509
# ╠═be52012c-4a74-4b3d-bdaf-ce380c8b08b0
# ╠═49af4bc2-be22-48b7-a08f-e5db775a64dc
# ╠═48d56784-63e1-408e-8b0b-0dc1d858b836
# ╟─c33a5eb8-ae22-412c-90ef-9cb1ab385ecb
# ╠═672984c0-eff8-4156-9ce8-0d2f919a4d4c
# ╟─b6cf83e6-043e-452a-ba72-2e5911c0cffa
# ╠═00c874fc-f087-4c6a-b653-feb5080ff89d
# ╠═274013ad-c939-467f-8dfd-cc199340bd2c
# ╠═2c443306-d607-42ff-b36a-4822b61e0846
# ╠═5aab77a4-ef85-4676-aeb6-eb2ce0aaeb22
# ╠═b92163fd-7118-4a9b-b0b0-27bbca819b7f
# ╠═4e579820-e404-4a44-900f-cbd588dd931e
# ╠═0b004661-a463-4f00-bad9-3532245000da
