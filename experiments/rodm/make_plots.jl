### A Pluto.jl notebook ###
# v0.17.7

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

# ╔═╡ f83c0c3a-200d-11ec-2135-1fa17a0df82a
begin
	using Pkg; Pkg.activate("")
	using Convex
	using Statistics
end

# ╔═╡ 5d212716-9661-4d67-9382-4b0c775d6d85
using CairoMakie

# ╔═╡ 56c26e7e-0ddb-45f8-b589-ee082f4b5c19
using StatsBase: trim

# ╔═╡ a6af42c6-1ce3-4092-a5a5-fde864a55043
using Random

# ╔═╡ d7267646-dfb8-417d-8d08-61d3b0308072
using StatsBase: harmmean

# ╔═╡ 6c69224a-5c5c-4c62-a26a-e5b3d24fc16b
using PlutoUI

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
# dpi = 100

# ╔═╡ 3b69d7a4-2ca3-484d-b86f-4801414fce53
# theme(:default,
# 	label=nothing,
# 	tickfont=(:Times, 8),
# 	guidefont=(:Times, 8),
# 	legendfont=(:Times, 8),
# 	titlefont=(:Times, 8, :bold),
# 	dpi=dpi,
# 	frame=:box,
# 	grid=false,
# )

# ╔═╡ 49454659-4260-48ed-828b-d5411c9f5ca3
md"""
## Main figure

**TODO**

- [ ] Make heatmaps
- [ ] _Maybe_ Replace prediction error with generation mix after storage
- [ ] _Maybe_ Show dollars/mTCO2 reduced for storage investment, renewabale investment
"""

# ╔═╡ 90a254ce-86ec-4340-8761-2c42471a0ae5
figwidth = 469.75502 / 72

# ╔═╡ 96bcea19-a6ec-4505-9f8f-44f3ff15c82a
colors = [:Yellow, :Brown, :Blue]

# ╔═╡ 408baa62-4953-4ac2-8cdc-7da617bf6bea
main_plt = let
	l = @layout [a{0.3w} b]
	main_plt = plot(
		plt_violin, plt_demand_mefs,
		size = dpi .* (figwidth, 2),
		layout = l,
		bottom_margin = 10Plots.pt,
		left_margin = 10Plots.pt
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
# let
# 	plot(
# 		[dotplot(ΔE_err[alg][:co2] / mean(abs, ΔE[:co2]), ms=1.1) for alg in keys(me)]...
# 	)
# end

# ╔═╡ b44cca3a-04ee-43e6-9bc6-9826eb2c683a
μ_data = Dict(c => ΔE[c] ./ (ΔD .+ 0.9) for c in chems)

# ╔═╡ c7487d0d-8f4f-4ac4-b551-9921bd20d13b
num_days = length(ΔD) ÷ 24

# ╔═╡ 768363dd-3942-49c4-9a12-ebcd8356a36b
md"""
## Panel ?: Violin plot
"""

# ╔═╡ a0a8604b-8c8f-44c9-9dcc-b9eb4601757b
let
	σ = norm(ΔE[:co2] .- mean(ΔE[:co2])) / length(ΔE[:co2])
	
	err = mean(abs, ΔE[:co2] .- mean(ΔE[:co2])) / mean(abs, ΔE[:co2])
	
	mean(abs, ΔE[:co2])
end

# ╔═╡ d20f4766-74b2-4302-8cfd-e0e4c370e502
metric(e) = abs(e)

# ╔═╡ 1ae2b093-bed5-4e9c-8a10-e066e449d8e8
ΔE_err_abs = Dict(alg => Dict(chem => metric.(ΔE_err[alg][chem]) for chem in chems) for alg in keys(me))

# ╔═╡ 18237910-8535-466f-a213-ae4a16e8bd80
l1_errs = Dict(alg => Dict(chem => mean(ΔE_err_abs[alg][chem]) for chem in chems) for alg in keys(me))

# ╔═╡ 69453214-9fdb-4b97-9a16-eb2a73342079
ΔE_err_hourly = Dict(alg => Dict(chem => 
		mean([ΔE_err_abs[alg][chem][(d-1)*24+1 : d*24] for d in 1:num_days])
		for chem in chems) for alg in keys(me))

# ╔═╡ d01beaff-a7c6-4f3e-9dfe-1713929a4d92
# plt_violin = let
# 	plts = []
# 	for c in chems
# 		Random.seed!(1234)
# 		num_full = length(ΔE_err_abs[:ed][c])
# 		num_dots = 500
# 		dots = rand(1:num_full, num_dots)

# 		_xc = repeat([1, 2, 3], inner=num_dots)
# 		_yc = 100 * [ΔE_err_abs[:da][c][dots]; ΔE_err_abs[:ed][c][dots]; ΔE_err_abs[:uc][c][dots]] / mean(metric, ΔE[c])

# 		plt = violin(_xc, _yc)
# 		boxplot!(_xc, _yc, alpha=0.75, lw=2, ms=0)
# 		dotplot!(_xc, _yc, c=:Black, ms=2)

# 		make_pc = x -> round(100x, digits=1)
# 		xticklabels = (
# 			"DA\n$(make_pc(l1_errs[:da][c] / mean(metric, ΔE[c])))%", 
# 			"ED\n$(make_pc(l1_errs[:ed][c] / mean(metric, ΔE[c])))%", 
# 			"UC\n$(make_pc(l1_errs[:uc][c] / mean(metric, ΔE[c])))%"
# 		)
# 		plot!(
# 			xticks=((1, 2, 3), xticklabels),
# 			#ylabel="Normalized Error",
# 			ylim=(-0.00, quantile(_yc, 0.9)),
# 			size=(300, 300),
# 			title=string(c),
# 			yticks=(0:25:75, map(x -> string(x) * "%", 0:25:75))
# 		)
		
# 		push!(plts, plt)
# 	end
	
# 	plot!(plts[1], title="Normalized Error")
# end

# ╔═╡ 7f4c72b9-96b7-4cfe-b7aa-552ac20d68d9
md"""
## Error per hour
"""

# ╔═╡ c01e3eda-f4a1-43ac-962e-c7aae2bb2d4e
# plt_eph = let
# 	c = :co2
	
# 	plot(size=(600, 200))
	
# 	normalization = hourly_avg_abs_ΔE
# 	for (i, alg) in enumerate(keys(me))
# 		bar!((1:24) .+ 0.3*(i-1), ΔE_err_hourly[alg][c] ./ normalization,
# 			lw=1, label=uppercase.(string(alg)), bar_width=0.25, c=colors[i])
# 	end
	
# 	plot!(
# 		legend=:outertopright, 
# 		xlim=(0.5, 25),
# 		xticks=6:6:18,
# 		ylim=(0, Inf),
# 		xlabel="Hour",
# 		ylabel="Absolute Error",
# 		bottom_margin=10Plots.pt,
# 	)
# end

# ╔═╡ 8da21a6f-e639-49d8-be1f-c7e594d9f6b9
robust_avg(xs) = [mean(trim(map(x -> x[i], xs), prop=0.05)) for i in 1:length(xs[1])]

# ╔═╡ 4fde0a54-36db-433c-9bce-70c163d6a4f4
hourly_avg_abs_ΔE = robust_avg([abs.(ΔE[:co2][(d-1)*24+1:d*24]) for d in 1:num_days])

# ╔═╡ f90a7b6d-2591-4b6a-b4b8-10f866368759
μ_hourly = Dict(alg => Dict(chem =>
		robust_avg([me[alg][chem][(d-1)*24+1 : d*24] for d in 1:num_days])
		for chem in chems) for alg in keys(me))

# ╔═╡ 681aaea2-1232-42ed-8945-cd019c0b40bb
μ_data_hourly = Dict(c => robust_avg([μ_data[c][(d-1)*24+1:d*24] for d in 1:num_days]) for c in chems)

# ╔═╡ c43f9b7c-ace1-4873-ae9a-cc0d9d8d254d
# plt_mef_time = let
# 	c = :co2
	
# 	plot(size=(600, 200))
	
# 	plot!(1:24, μ_data_hourly[c], lw=2, ls=:dash, label="Data", c=:Black)
# 	for (i, alg) in enumerate(keys(me))
# 		plot!(1:24, μ_hourly[alg][c], 
# 			lw=2, 
# 			label=uppercase.(string(alg)),
# 			c = colors[i],
# 		)
# 	end
	
# 	plot!(
# 		legend=:outertopright, 
# 		xlim=(1, 24),
# 		xticks=6:6:18,
# 		xlabel="Hour",
# 		ylabel="Average MEF",
# 		bottom_margin=10Plots.pt,
# 	)
# end

# ╔═╡ 5bb43074-7ea6-4da6-bf3e-5f26ae520b65
md"""
## Plot distribution of error per hour
"""

# ╔═╡ 7ccac729-a819-4bcc-94bf-14b3c72cfe0a
ΔE_hourly_vec = Dict(alg => Dict(chem => 
		hcat([ΔE_est[alg][chem][(d-1)*24+1 : d*24] for d in 1:num_days]...)
		for chem in chems) for alg in keys(me));

# ╔═╡ ab33c9c8-076a-4be6-9b1f-44b187f8b638
ΔE_err_hourly_vec = Dict(alg => Dict(chem => 
		hcat([ΔE_err[alg][chem][(d-1)*24+1 : d*24] for d in 1:num_days]...)
		for chem in chems) for alg in keys(me));

# ╔═╡ cb1e76cc-8a4b-4004-9527-9ed202725237
μ_hourly_vec = Dict(alg => Dict(chem =>
		hcat([me[alg][chem][(d-1)*24+1 : d*24] for d in 1:num_days]...)
		for chem in chems) for alg in keys(me));

# ╔═╡ 4e3af057-c4f2-463f-ba1a-3c236ff630db
μ_data_hourly_vec = Dict(c => hcat([μ_data[c][(d-1)*24+1:d*24] for d in 1:num_days]...) for c in chems);

# ╔═╡ 8f7f84cf-dfa2-433a-8506-28fc8793b94b
ΔE_data_hourly_vec = Dict(c => hcat([ΔE[c][(d-1)*24+1:d*24] for d in 1:num_days]...) for c in chems);

# ╔═╡ 8ecb77fe-bb5c-48a5-bd79-6af182001123
# let
# 	hour = density_hour
# 	c = :co2
	
# 	plt1 = plot(
# 		##xlim=(100, 1200),
# 		title="PDF of MEFs at hour $hour",
# 	)
# 	density!(ΔE_data_hourly_vec[c][hour, :], lw=2, label="data", c=:black, ls=:dash)
# 	for alg in keys(me)
# 		density!(ΔE_hourly_vec[alg][c][hour, :], lw=2, label=string(alg))
# 	end
# 	plot!()
	
# 	ax_lims = (minimum(ΔE_data_hourly_vec[c][hour, :] / 1e6) + 0.5, maximum(ΔE_data_hourly_vec[c][hour, :] / 1e6) + 0.5)
# 	plot(
# 		[marginalhist(
# 				ΔE_data_hourly_vec[c][hour, :] / 1e6, ΔE_hourly_vec[alg][c][hour, :] / 1e6,
# 				title=string(alg), nbin=30) 
# 		for alg in keys(me)]..., 
# 		layout=(1, 3), size=(650, 200), ylim=ax_lims, xlim=ax_lims,
# 	)
	
# 	plt1
# end

# ╔═╡ 731b2899-d1c9-4d2b-b841-2fe369a0461d
@bind density_hour Slider(1:24)

# ╔═╡ 60fc957d-21ca-43e1-b21b-0b8efd8931d4
density_hour

# ╔═╡ dc4b0662-3a99-4e1a-ae66-f723fcdc0204
# let hour = density_hour, c = :co2
	
# 	plot()
# 	for alg in keys(me)
# 		density!(
# 			ΔE_err_hourly_vec[alg][c][hour, :] / mean(abs, ΔE[c]), 
# 			label=uppercase(string(alg)), 
# 			lw=4, alpha=0.75
# 		)
# 	end
# 	vline!([0.0], lw=4, c=:Black) 
	
# 	plot!(
# 		#xlim=(-2, 2),
# 		title="PDF of Error for Hour $hour",
# 		size=(600, 200),
# 	)
# end

# ╔═╡ 59023006-6a80-4ee9-abc9-8e8ec8f4cfda
# let hour = density_hour
# 	plts = [
# 		scatter( 
# 			ΔE_data_hourly_vec[:co2][hour, :] ./ 1e6, 
# 			ΔE_hourly_vec[alg][:co2][hour, :] ./ 1e6,
# 			xlim=(-5, 5),
# 			ylim=(-5, 5),
# 			title=uppercase(string(alg)),
# 			ms=2,
# 			xlabel="Data",
# 			#smooth=true,
# 			#lw=4,
# 			#linealpha=0.75,
# 			#qqline=:R
# 		)
# 		for alg in keys(me)
# 	]
	
# 	plot!(plts[1], ylabel="Model")
# 	[plot!(plts[k], [-1000, 1000], [-1000, 1000], lw=4, alpha=0.3, c=:Black, ) for k in 1:3]
	
# 	plot(plts..., layout=(1, 3), size=(600, 200), link=:y, bottom_margin=10Plots.pt)
# end

# ╔═╡ 5b7aa658-e793-497c-a3ed-55c78599d1ae
ΔD_hourly = hcat([ΔD[(d-1)*24+1:d*24] for d in 1:num_days]...)

# ╔═╡ 02724af2-8024-4270-bbab-46f57ed86403
md"""
## LAST PLOT

Demand vs MEF
"""

# ╔═╡ 21754fbb-7ecb-4acd-9ff9-6bce91a1c14d
demands = df.carbon.demand[util.TIMES]

# ╔═╡ 2c33db45-9be8-4251-9755-5dc58682a855
demand_order = sortperm(demands)

# ╔═╡ 34aebe17-60a5-47e9-9f97-b06e32a4cc64
dd = demands[demand_order]

# ╔═╡ e9288b4a-4718-483c-900d-62890ce63514
c = :co2

# ╔═╡ d8911b94-d0f4-4acd-9df1-6c69df79b256
μd = μ_data[c][demand_order]

# ╔═╡ be480e0f-d8eb-4524-8b64-2dd15eb85a20
μd_alg = Dict(alg => me[alg][c][demand_order] for alg in keys(me))

# ╔═╡ 66a543b6-e030-45a6-909a-5e286ac82685
smooth(x, w, δ=0.25) = let
	n = length(x)
	vecs = [trim(x[i-w:i+w], prop=0.01) for i in (w+1):(n-w)]
	
	(
		mean.(vecs), 
		quantile.(vecs, 0.5 + δ), 
		quantile.(vecs, 0.5 - δ)
	)
end

# ╔═╡ 62abef96-709d-4413-abc1-0f0ff7bd0d55
# plt_demand_mefs = let
# 	w = floor(Int, 0.05 * length(μd))
# 	smth = x -> smooth(x, w)	
# 	_x = smth(dd)[1] / 1e3
	
# 	μds, μdsu, μdsl = smth(μd)
# 	plot(_x, μds, c=:black, ls=:dash, ribbon=(μdsu, μdsl), fillalpha=0.5, label="Data")
	
# 	for (i, alg) in enumerate(keys(me))
# 		μdsa, μdsau, μdsal = smth(μd_alg[alg])
# 		plot!(_x, μdsa, ribbon=(μdsau, μdsal), fillalpha=0.3, 
# 			label=uppercase(string(alg)), c=i
# 		)
# 	end
	
# 	plot!(
# 		xlabel="Demand [GWh]",
# 		title="Marginal Emissions [kg CO2 / MWh]",
# 		size=(450, 250),
# 		legend=:bottomleft,
# 	)
# end

# ╔═╡ d2040316-ee23-44df-9865-5bf9a64260ad
begin
	paper_theme = Theme()
	set_theme!(paper_theme)
end

# ╔═╡ 5dd2d8cb-042d-4719-b33f-1971a7c0bb6d
function get_color(f, i)
	f.scene.attributes.attributes[:palette].val.attributes[:color].val[i]
end

# ╔═╡ fd57e810-fc54-47d7-8210-defffcfd1591
set_alpha(c, a) = CairoMakie.ColorTypes.RGBA(c.r, c.g, c.b, a)

# ╔═╡ 5b80b8c6-36be-4fbd-970f-eb1273ce9b0e
main_figure = let
	f = Figure(resolution=(650, 250), fontsize=10)
	
	gl = f[1, 1] =  GridLayout()
	gr = f[1, 2] = GridLayout()
	
	
	# ===
	# Left side
	# ===
	c = :co2
	Random.seed!(1234)
	num_full = length(ΔE_err_abs[:ed][c])
	num_dots = 500
	dots = rand(1:num_full, num_dots)
	make_pc = x -> round(100x, digits=1)

	_xc = repeat([1, 2, 3], inner=num_dots)
	_yc = [ΔE_err_abs[:da][c][dots]; ΔE_err_abs[:ed][c][dots]; ΔE_err_abs[:uc][c][dots]] / mean(metric, ΔE[c])
	_cc = repeat([set_alpha(get_color(f, i), 0.25) for i in 1:3], inner=num_dots)
	xtick_labels = [
			"DA\n$(make_pc(l1_errs[:da][c] / mean(metric, ΔE[c])))%", 
			"ED\n$(make_pc(l1_errs[:ed][c] / mean(metric, ΔE[c])))%", 
			"UC\n$(make_pc(l1_errs[:uc][c] / mean(metric, ΔE[c])))%"
	]
	
	
	axl = Axis(
		gl[1, 1], 
		xgridvisible=false, 
		ygridvisible=false,
		xticks=(1:3, xtick_labels),
		yticks=(0:0.25:0.75, [string(make_pc(i)) * "%" for i in 0:0.25:0.75]),
	)
	
	violin!(_xc, _yc, color=(:black, 0.25))
	boxplot!(
		_xc, _yc, 
		color=_cc, 
		show_outliers=false,
		strokewidth=2,
	)
	ylims!(axl, 0, 0.75)
	axl.ylabel = "Normalized Error"
	
	
	
	
	
	# ===
	# Right side
	# ===
	
	w = floor(Int, 0.05 * length(μd))
	smth = x -> smooth(x, w)	
	_x = smth(dd)[1] / 1e3
	μds, μdsu, μdsl = smth(μd)
	
	
	axr = []
	for (i, alg) in enumerate([:da, :ed, :uc])
		μdsa, μdsau, μdsal = smth(μd_alg[alg])
		
		# Plot
		push!(axr, Axis(
				gr[i, 1], 
				xgridvisible=false, 
				ygridvisible=false,
				xticks=([], []),
			))
		ylims!(axr[i], 350 / 1e3, 900 / 1e3)
		xlims!(axr[i], minimum(_x), maximum(_x))
		
		band!(_x, μdsl / 1e3, μdsu / 1e3, color=(:black, 0.25))
	 	lines!(_x, μds / 1e3, color=:black)
		
		band!(_x, μdsal / 1e3, μdsau / 1e3, color=set_alpha(get_color(f, i), 0.25))
		lines!(_x, μdsa / 1e3, color=get_color(f, i), label=uppercase(string(alg)))
	end
	
	axr[3].xlabel = "Demand [GWh]"
	axr[2].ylabel = "MEF [t / MWh]"
	axr[3].xticks = collect(20:10:60)
	rowgap!(gr, 6)
	
	
	
	
	
	# ===
	# Final tweaks
	# ===
	for (label, layout) in zip(["A", "B"], [gl, gr])
		Label(
			layout[1, 1, TopLeft()], 
			label,
			textsize = 18,
			font = "DejaVu Sans",
			padding = (0, -18, 2, 0),
			halign = :right
		)
	end
	
	colsize!(f.layout, 1, Auto(0.5))
	
	f
end

# ╔═╡ 12ae050a-8f3d-435a-b8e9-1c70255596aa
save("makie_mefs_rodm_figure.pdf", main_figure, pt_per_unit=1)

# ╔═╡ Cell order:
# ╠═f83c0c3a-200d-11ec-2135-1fa17a0df82a
# ╟─58555d65-d842-423c-8380-f6d38149f7d5
# ╠═dcb739ab-7837-4807-8df1-b471492f17a2
# ╠═5d212716-9661-4d67-9382-4b0c775d6d85
# ╠═3065f8e6-eebc-41ab-8516-e3f019bcfba2
# ╠═3b69d7a4-2ca3-484d-b86f-4801414fce53
# ╠═49454659-4260-48ed-828b-d5411c9f5ca3
# ╠═90a254ce-86ec-4340-8761-2c42471a0ae5
# ╠═96bcea19-a6ec-4505-9f8f-44f3ff15c82a
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
# ╠═69453214-9fdb-4b97-9a16-eb2a73342079
# ╠═b44cca3a-04ee-43e6-9bc6-9826eb2c683a
# ╠═c7487d0d-8f4f-4ac4-b551-9921bd20d13b
# ╠═56c26e7e-0ddb-45f8-b589-ee082f4b5c19
# ╠═4fde0a54-36db-433c-9bce-70c163d6a4f4
# ╠═f90a7b6d-2591-4b6a-b4b8-10f866368759
# ╠═681aaea2-1232-42ed-8945-cd019c0b40bb
# ╟─768363dd-3942-49c4-9a12-ebcd8356a36b
# ╠═a6af42c6-1ce3-4092-a5a5-fde864a55043
# ╠═a0a8604b-8c8f-44c9-9dcc-b9eb4601757b
# ╠═d20f4766-74b2-4302-8cfd-e0e4c370e502
# ╟─d01beaff-a7c6-4f3e-9dfe-1713929a4d92
# ╟─7f4c72b9-96b7-4cfe-b7aa-552ac20d68d9
# ╟─c01e3eda-f4a1-43ac-962e-c7aae2bb2d4e
# ╠═d7267646-dfb8-417d-8d08-61d3b0308072
# ╠═8da21a6f-e639-49d8-be1f-c7e594d9f6b9
# ╟─c43f9b7c-ace1-4873-ae9a-cc0d9d8d254d
# ╠═5bb43074-7ea6-4da6-bf3e-5f26ae520b65
# ╠═7ccac729-a819-4bcc-94bf-14b3c72cfe0a
# ╠═ab33c9c8-076a-4be6-9b1f-44b187f8b638
# ╠═cb1e76cc-8a4b-4004-9527-9ed202725237
# ╠═4e3af057-c4f2-463f-ba1a-3c236ff630db
# ╠═8f7f84cf-dfa2-433a-8506-28fc8793b94b
# ╠═6c69224a-5c5c-4c62-a26a-e5b3d24fc16b
# ╟─60fc957d-21ca-43e1-b21b-0b8efd8931d4
# ╟─8ecb77fe-bb5c-48a5-bd79-6af182001123
# ╠═731b2899-d1c9-4d2b-b841-2fe369a0461d
# ╟─dc4b0662-3a99-4e1a-ae66-f723fcdc0204
# ╟─59023006-6a80-4ee9-abc9-8e8ec8f4cfda
# ╠═5b7aa658-e793-497c-a3ed-55c78599d1ae
# ╠═02724af2-8024-4270-bbab-46f57ed86403
# ╠═21754fbb-7ecb-4acd-9ff9-6bce91a1c14d
# ╠═2c33db45-9be8-4251-9755-5dc58682a855
# ╠═34aebe17-60a5-47e9-9f97-b06e32a4cc64
# ╠═e9288b4a-4718-483c-900d-62890ce63514
# ╠═d8911b94-d0f4-4acd-9df1-6c69df79b256
# ╠═be480e0f-d8eb-4524-8b64-2dd15eb85a20
# ╠═66a543b6-e030-45a6-909a-5e286ac82685
# ╟─62abef96-709d-4413-abc1-0f0ff7bd0d55
# ╠═d2040316-ee23-44df-9865-5bf9a64260ad
# ╟─5dd2d8cb-042d-4719-b33f-1971a7c0bb6d
# ╟─fd57e810-fc54-47d7-8210-defffcfd1591
# ╠═5b80b8c6-36be-4fbd-970f-eb1273ce9b0e
# ╠═12ae050a-8f3d-435a-b8e9-1c70255596aa
