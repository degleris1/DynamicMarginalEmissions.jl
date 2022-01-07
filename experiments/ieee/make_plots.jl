### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 85394d0a-6f17-11ec-04e7-89c834098e0c
begin
	using Pkg; Pkg.activate("")
	
	using HDF5
	using CairoMakie
	using CarbonNetworks
end

# ╔═╡ c041f8ff-df9a-4b6a-a73b-b125c3de4495
md"""
### Load data
"""

# ╔═╡ 9874d30f-c591-48a4-8118-a523bc758509
f = h5open("/Users/degleris/Data/carbon_networks/fig2_data.h5") do f
	r = Dict()

	for k in keys(f)
		r[k] = read(f[k])
	end

	r
end;

# ╔═╡ c33a5eb8-ae22-412c-90ef-9cb1ab385ecb
md"""
### Plot A: MEFs across the network
"""

# ╔═╡ 672984c0-eff8-4156-9ce8-0d2f919a4d4c
figure_nodal_mefs = let
	nodal_mefs = f["MEFs"] / 1000
	Δd = f["exp_vs_th_x"]
	ΔE = f["exp_vs_th_y_exp"]
	ΔE_est = f["exp_vs_th_y_th"]

	# Setup
	size_inches = (6.5, 2)
	size_pt = 72 .* size_inches
	fig = Figure(resolution=size_pt ./ 0.75, fontsize=10)

	# Left Panel
	ax_left = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false)
	barplot!(ax_left, nodal_mefs, strokewidth=0.5)
	xlims!(ax_left, 0, 31)
	ax_left.ylabel = "MEF [mTCO2 / MWh]"
	ax_left.xlabel = "Node"

	# Right Panel
	ax_right = Axis(fig[1, 2], xgridvisible=false, ygridvisible=false)
	lines!(ax_right, Δd, ΔE, label="True", linewidth=3)
	lines!(ax_right, Δd, ΔE_est, label="Estimate", linestyle=:dot, linewidth=3)
	xlims!(ax_right, 0.9, 1.1)
	ylims!(ax_right, 0.98, 1.06)
	ax_right.ylabel = "ΔE [%]"
	ax_right.xlabel = "ΔD [%]"
	axislegend(ax_right, position=:lt, padding=(4.0, 4.0, 2.0, 2.0), rowgap=0)

	# Add panel labels
	for (label, layout) in zip(["A", "B"], [fig[1, 1], fig[1, 2]])
	    Label(
			layout[1, 1, TopLeft()], 
			label,
	        textsize=18,
			font="Noto Sans Bold",
	        padding=(0, 5, 5, 0),
	        halign=:right
		)
	end

	# Decrease size of second column
	colsize!(fig.layout, 2, Auto(0.5))

	fname = "/Users/degleris/Documents/Research/CarbonNetworks/img/network_mefs.pdf"
	save(fname, fig, pt_per_unit=0.75)
	fig
end

# ╔═╡ b6cf83e6-043e-452a-ba72-2e5911c0cffa
md"""
### Plot B: MEFs over time
"""

# ╔═╡ 4e579820-e404-4a44-900f-cbd588dd931e
figure_temporal = let
	hours = f["MEF_vs_t_x"]
	total_mefs = [f["MEF_vs_t_y1"], f["MEF_vs_t_y2"]] / 1000
	nodes = Dict(1 => 21, 2 => 23)
	mef_grids = [f["hm_$(nodes[j])_$(k)"] for j in 1:2, k in 1:3] ./ 1000
	cmax = maximum(maximum.(mef_grids))
	cmin = minimum(minimum.(mef_grids))
	
	# Setup
	lw = 3
	size_inches = (6.5, 3)
	size_pt = 72 .* size_inches
	
	fig = Figure(resolution=size_pt ./ 0.75, fontsize=10)
	#grid_left = fig[1, 1] = GridLayout()
	#grid_right = fig[1, 2] = GridLayout()

	
	# Left panel
	gl_axes = [
		Axis(fig[k, 1], xgridvisible=false, ygridvisible=false) 
		for k in 1:2
	]
	for i in 1:2
		xlims!(gl_axes[i], 1, 24)
		gl_axes[i].xticks = 0:6:24
		
		lines!(gl_axes[i], total_mefs[i][:, 1], linewidth=lw)
		lines!(gl_axes[i], total_mefs[i][:, 2], linewidth=lw)
	end
	Label(fig[1:2, 1, Left()], "MEF [mTCO2 / MWh]", 
		rotation=pi/2, 
		padding=(0, 35, 5, 0)
	)
	gl_axes[2].xlabel = "Hour"

	
	# Right panel
	gr_axes = [
		Axis(fig[j, k+1]) 
		for j in 1:2, k in 1:3
	]
	for (j, k) in Iterators.product(1:2, 1:3)
		
		heatmap!(gr_axes[j, k], f["hm_$(nodes[j])_$(k)"], colorrange=(cmin, cmax))
		gr_axes[j, k].xticks = 0:6:24
		gr_axes[j, k].yticks = 0:6:24

		if k != 1
			gr_axes[j, k].ytickformat = xs -> ["" for x in xs]
		end
		if j != 2
			gr_axes[j, k].xtickformat = xs -> ["" for x in xs]
		end
	end
	cb = Colorbar(
		fig[1:2, 5], 
		colormap = :viridis, 
		limits=(cmin, cmax), 
		label="MEF [mTCO2 / MWh]"
	)
	
	Label(fig[1:2, 2, Left()], "Emissions Time", 
		rotation=pi/2, 
		padding=(0, 30, 0, 0)
	)
	Label(fig[2, 2:4, Bottom()], "Consumption Time", 
		padding=(0, 0, 0, 24)
	)

	
	# Tweaks
	# Grow left figure column
	colsize!(fig.layout, 1, Auto(3))
	
	# Make heatmaps square
	for j in 2:4
		colsize!(fig.layout, j, Aspect(1, 1))
	end

	# Add panel labels
	for (label, layout) in zip(["A", "B"], [fig[1, 1], fig[1, 2]])
	    Label(
			layout[1, 1, TopLeft()], 
			label,
	        textsize=18,
			font="Noto Sans Bold",
	        padding=(0, 5, 5, 0),
	        halign=:right
		)
	end

	
	
	fname = "/Users/degleris/Documents/Research/CarbonNetworks/img/mef_heatmaps.pdf"
	save(fname, fig, pt_per_unit=0.75)
	fig
end

# ╔═╡ 0b004661-a463-4f00-bad9-3532245000da
maximum(maximum.([f["hm_$(Dict(1=>21, 2=>23)[j])_$(k)"] for j in 1:2, k in 1:3]))

# ╔═╡ Cell order:
# ╠═85394d0a-6f17-11ec-04e7-89c834098e0c
# ╟─c041f8ff-df9a-4b6a-a73b-b125c3de4495
# ╠═9874d30f-c591-48a4-8118-a523bc758509
# ╟─c33a5eb8-ae22-412c-90ef-9cb1ab385ecb
# ╠═672984c0-eff8-4156-9ce8-0d2f919a4d4c
# ╟─b6cf83e6-043e-452a-ba72-2e5911c0cffa
# ╠═4e579820-e404-4a44-900f-cbd588dd931e
# ╠═0b004661-a463-4f00-bad9-3532245000da
