### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ f83c0c3a-200d-11ec-2135-1fa17a0df82a
begin
	using Pkg; Pkg.activate("")
	using Convex
	using Statistics
end

# ╔═╡ 5d212716-9661-4d67-9382-4b0c775d6d85
using Plots

# ╔═╡ 8da59504-ad52-466f-8258-bc774e107faf
using StatsPlots

# ╔═╡ a6af42c6-1ce3-4092-a5a5-fde864a55043
using Random

# ╔═╡ 58555d65-d842-423c-8380-f6d38149f7d5
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
end

# ╔═╡ dcb739ab-7837-4807-8df1-b471492f17a2
util = ingredients("util.jl")

# ╔═╡ 3065f8e6-eebc-41ab-8516-e3f019bcfba2
dpi = 100

# ╔═╡ 3b69d7a4-2ca3-484d-b86f-4801414fce53
theme(:default,
	label=nothing,
	tickfont=(:Times, 8),
	guidefont=(:Times, 8),
	legendfont=(:Times, 8),
	titlefont=(:Times, 8, :bold),
	dpi=dpi,
	frame=:box,
	grid=false,
)

# ╔═╡ 49454659-4260-48ed-828b-d5411c9f5ca3
md"""
## Main figure

**TODO**

- [ ] Make heatmaps
- [ ] _Maybe_ Replace prediction error with generation mix after storage
- [ ] _Maybe_ Show dollars/mTCO2 reduced for storage investment, renewabale investment
"""

# ╔═╡ 408baa62-4953-4ac2-8cdc-7da617bf6bea
begin
	main_plt = plot(
		plta, pltb,
		size = dpi .* (8, 3),
		layout = (1, 2),
		bottom_margin = 10Plots.pt,
	)
	
	savefig(main_plt, "../../img/rodm_full.pdf")
	savefig(main_plt, "../../img/rodm_full.png")
	main_plt
end

# ╔═╡ 6dcd2b42-ecfc-4d7c-afc9-13fda3807827
md"""
## Load data

Comparing RODM, RODM + MDT, and RODM + ramping
"""

# ╔═╡ 98c7c6d5-8c08-42ef-9652-798f9c26137d
ed_config, df, ed_data, ed_results = util.load_results_static("base");

# ╔═╡ a6c15810-1286-4442-915f-a0665add81a6
uc_config, _, _, uc_data, uc_results = util.load_results_dynamic("unit");

# ╔═╡ 5ace0589-6462-46f9-8411-02a8bdf71a71
md"""
#### Compute MEFs
"""

# ╔═╡ 04e905b2-325a-404d-a71b-3c2a348889a7
chems = util.CHEMICALS

# ╔═╡ ee089f14-3eea-4907-8540-bf80af0e488d
num_chem = length(chems)

# ╔═╡ 7a457bfe-51d5-44f6-aa9a-f60c00fb1f3f
me_ed = Dict(
	chem => [r.dq[c] for r in ed_results] 
	for (c, chem) in enumerate(keys(ed_data[1].q))
)

# ╔═╡ 213e5fe9-6d37-4913-9b51-132b2d1cec3a
me_uc = Dict(
	chem => vcat([sum.(r.dq[c]) for r in uc_results]...) 
	for (c, chem) in enumerate(chems)
)

# ╔═╡ b3b350cf-3977-4d12-a610-e35a26505112
me_da = Dict(
	chem => df.carbon[util.TIMES, Symbol(chem, :_marg)]
	for chem in chems
)

# ╔═╡ 3a7f292a-0cd3-472a-b968-342841a22d06
me = Dict(
	:ed => me_ed,
	:uc => me_uc,
	:da => me_da,
)

# ╔═╡ a359e195-6121-4438-aae8-0a3b2a1f1aa5
md"""
#### Compute true (historical) MEFs
"""

# ╔═╡ 9c886bfe-7931-4cb7-8ee5-e3ef18b3de8f
ΔD = diff(df.carbon.demand)[util.TIMES]

# ╔═╡ d2d812ff-a7e6-4cf2-9dae-8144ed9ccb6d
ΔE = Dict(
	chem => diff(df.carbon[:, Symbol(chem, :_tot)])[util.TIMES]
	for chem in chems
)

# ╔═╡ be43c4c3-c222-4792-b6ab-3f0efeef117b
ΔE_est = Dict(alg => Dict(chem => me[alg][chem] .* ΔD for chem in chems) for alg in keys(me))

# ╔═╡ 74427ed9-1daf-4606-9bcb-ebcd4e67e42b
ΔE_err = Dict(alg => Dict(chem => ΔE_est[alg][chem] - ΔE[chem] for chem in chems) for alg in keys(me))

# ╔═╡ a64d1047-3174-4f05-b7c7-3885c003c47d
Dict(alg => Dict(chem => mean(ΔE_err[alg][chem]) / std(ΔE_err[alg][chem]) for chem in chems) for alg in keys(me))

# ╔═╡ ef472b30-75cf-4e5d-a851-9f70b129686a
let
	plot(
		[dotplot(ΔE_err[alg][:co2] / mean(abs, ΔE[:co2]), ms=1.1) for alg in keys(me)]...
	)
end

# ╔═╡ 1ae2b093-bed5-4e9c-8a10-e066e449d8e8
ΔE_err_abs = Dict(alg => Dict(chem => abs.(ΔE_err[alg][chem]) for chem in chems) for alg in keys(me))

# ╔═╡ 18237910-8535-466f-a213-ae4a16e8bd80
l1_errs = Dict(alg => Dict(chem => mean(ΔE_err_abs[alg][chem]) for chem in chems) for alg in keys(me))

# ╔═╡ 768363dd-3942-49c4-9a12-ebcd8356a36b
md"""
## Panel ?: Violin plot
"""

# ╔═╡ d01beaff-a7c6-4f3e-9dfe-1713929a4d92
let
	plts = []
	for c in chems
		Random.seed!(1234)
		num_full = length(ΔE_err_abs[:ed][c])
		num_dots = 500
		dots = rand(1:num_full, num_dots)

		_xc = repeat([1, 2, 3], inner=num_dots)
		_yc = [ΔE_err_abs[:ed][c][dots]; ΔE_err_abs[:da][c][dots]; ΔE_err_abs[:uc][c][dots]] / mean(abs, ΔE[c])

		plt = violin(_xc, _yc)
		boxplot!(_xc, _yc, alpha=0.75, lw=2, ms=0)
		dotplot!(_xc, _yc, c=:Black, ms=0.8)

		make_pc = x -> round(100x, digits=1)
		xticklabels = (
			"ED\n$(make_pc(l1_errs[:ed][c] / mean(abs, ΔE[c])))%", 
			"MDT\n$(make_pc(l1_errs[:da][c] / mean(abs, ΔE[c])))%", 
			"UC\n$(make_pc(l1_errs[:uc][c] / mean(abs, ΔE[c])))%"
		)
		plot!(
			xticks=((1, 2, 3), xticklabels),
			ylabel="Normalized Error",
			ylim=(-0.1, quantile(_yc, 0.93)),
			size=(300, 300),
			title=string(c),
		)
		
		push!(plts, plt)
	end
	
	plot(plts..., layout=(1, 3), size=(600, 200),
		bottom_margin=10Plots.pt,
	)
end

# ╔═╡ Cell order:
# ╠═f83c0c3a-200d-11ec-2135-1fa17a0df82a
# ╟─58555d65-d842-423c-8380-f6d38149f7d5
# ╠═dcb739ab-7837-4807-8df1-b471492f17a2
# ╠═5d212716-9661-4d67-9382-4b0c775d6d85
# ╠═3065f8e6-eebc-41ab-8516-e3f019bcfba2
# ╠═3b69d7a4-2ca3-484d-b86f-4801414fce53
# ╠═49454659-4260-48ed-828b-d5411c9f5ca3
# ╠═408baa62-4953-4ac2-8cdc-7da617bf6bea
# ╟─6dcd2b42-ecfc-4d7c-afc9-13fda3807827
# ╠═98c7c6d5-8c08-42ef-9652-798f9c26137d
# ╠═a6c15810-1286-4442-915f-a0665add81a6
# ╟─5ace0589-6462-46f9-8411-02a8bdf71a71
# ╠═04e905b2-325a-404d-a71b-3c2a348889a7
# ╠═ee089f14-3eea-4907-8540-bf80af0e488d
# ╠═7a457bfe-51d5-44f6-aa9a-f60c00fb1f3f
# ╠═213e5fe9-6d37-4913-9b51-132b2d1cec3a
# ╠═b3b350cf-3977-4d12-a610-e35a26505112
# ╠═3a7f292a-0cd3-472a-b968-342841a22d06
# ╟─a359e195-6121-4438-aae8-0a3b2a1f1aa5
# ╠═9c886bfe-7931-4cb7-8ee5-e3ef18b3de8f
# ╠═d2d812ff-a7e6-4cf2-9dae-8144ed9ccb6d
# ╠═be43c4c3-c222-4792-b6ab-3f0efeef117b
# ╠═74427ed9-1daf-4606-9bcb-ebcd4e67e42b
# ╠═a64d1047-3174-4f05-b7c7-3885c003c47d
# ╠═ef472b30-75cf-4e5d-a851-9f70b129686a
# ╠═1ae2b093-bed5-4e9c-8a10-e066e449d8e8
# ╠═18237910-8535-466f-a213-ae4a16e8bd80
# ╠═768363dd-3942-49c4-9a12-ebcd8356a36b
# ╠═8da59504-ad52-466f-8258-bc774e107faf
# ╠═a6af42c6-1ce3-4092-a5a5-fde864a55043
# ╠═d01beaff-a7c6-4f3e-9dfe-1713929a4d92
