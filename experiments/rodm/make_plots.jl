### A Pluto.jl notebook ###
# v0.15.1

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

# ╔═╡ 148f5d6f-dc6d-4664-beac-6a0098422cb4
using Dates

# ╔═╡ e2a7f5a6-3a9b-44d1-b6ef-726b9870aa59
using Printf

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

# ╔═╡ 6dcd2b42-ecfc-4d7c-afc9-13fda3807827
md"""
## Load data

Comparing RODM, RODM + MDT, and RODM + ramping
"""

# ╔═╡ d0a30c4e-b675-439b-8bc0-adec6126360e
dyn_config, df, data, dynamic_data, results = util.load_results_dynamic("ramp4");

# ╔═╡ 98c7c6d5-8c08-42ef-9652-798f9c26137d
static_config, _, _, static_results = util.load_results_static("base");

# ╔═╡ 63ac3134-90e6-4351-90de-54ae6c21dee1
begin
	storage_results = []
	for i in 1:length(util.STORAGE_CONFIG)
		cfg, _, _, _, r = util.load_results_dynamic("storage$i")
		push!(storage_results, (cfg=cfg, r=r))
	end
end

# ╔═╡ a359e195-6121-4438-aae8-0a3b2a1f1aa5
md"""
#### Compute true (historical) MEFs
"""

# ╔═╡ 118dc77a-b4d7-44c5-8eda-f9a4b29df23a
df.gen.is_coal

# ╔═╡ 9c78666f-f03c-4b1b-9a8c-8b784eeb0b65
bar(
	results[25].g[34] ./ data[25*7*24].θ.gmax, 
	size=(600, 200), group=df.gen.is_coal .+ 1, lw=0, legend=false
)

# ╔═╡ 9c886bfe-7931-4cb7-8ee5-e3ef18b3de8f
demand = df.carbon.demand

# ╔═╡ 8c22b4e7-6ecf-4188-a34b-2b8e4d5f9b58
co2 = df.carbon.co2_tot

# ╔═╡ 93a455cc-c620-463e-90da-654bb96dca8b
mef_hist = (diff(co2) ./ diff(demand))[util.TIMES]

# ╔═╡ 33d7aace-51f7-411c-bef3-79e6c672152d
valid_times = findall(x -> 0 < x < Inf, mef_hist)

# ╔═╡ 0c58311b-3001-4307-b8d3-e9365a0829c1
get_err = mef_est -> (
	mean(abs, (mef_hist - mef_est)[valid_times]) 
	/ mean(mef_hist[valid_times])
)

# ╔═╡ 3330f77a-73a3-4a8e-90b3-27824d356514
md"""
#### Load MEF estimates
"""

# ╔═╡ 002f2a27-8489-4521-8ed9-660ec6cd5760
mefs_ed = vcat(getfield.(static_results, :dq)...)[:, 1]

# ╔═╡ 94ccf5d2-fa1f-4b08-aee5-51340bad7e56
mefs_mdt = df.carbon.co2_marg[util.TIMES]

# ╔═╡ 54db5151-a50c-43a0-aa2e-f873a1de30b1
mefs_ramp = vcat([sum.(r.dq[1]) for r in results]...)

# ╔═╡ bdac9b46-6773-4546-8f03-4a4e1e3c4b88
md"""
#### Compute errors
"""

# ╔═╡ dc592388-31ab-4428-9bac-9671f12b0c5c
err_ed = get_err(mefs_ed)

# ╔═╡ e9623b71-b6de-455d-8a26-c41be6d58f8c
err_mdt = get_err(mefs_mdt)

# ╔═╡ 9a65cc16-734c-4101-b65c-cda75aa9b75e
err_ramp = get_err(mefs_ramp)

# ╔═╡ 74dc860b-0ff0-4cf4-90e1-0505a9d952a4
md"""
## Panel A1: Load profile
"""

# ╔═╡ 36c864b0-abf3-48e2-8c99-48669f142896
load_profile_start = Date(2000, 6, 1)

# ╔═╡ e3aeecd6-0a14-484b-aa28-54a360214dcb
load_profile_stop = Date(2000, 6, 7)

# ╔═╡ 0502409b-2a9f-4fc5-ac75-057de60736f6
load_profile_start + Day(1)

# ╔═╡ 3d83859f-391f-45cd-bf8a-99c9e1a9cc80
begin
	x1 = (dayofyear(load_profile_start)-1)*24+1
	x2 = dayofyear(load_profile_stop)*24
	
	
	plta1 = plot(size = dpi .* (4, 1.5))
	
	plot!(x1:x2, df.carbon.demand[x1:x2] / 1e3, lw=2)
	
	plot!(ylabel="Demand [GW]", xlabel="Day")
	plot!(xticks=(x1:24:x2, 1:7))
	plot!(xlim=(x1, x2), ylim=(0, Inf))
	plot!(bottom_margin=10Plots.pt)
end

# ╔═╡ 1ad7f79b-b1ae-4f2a-9ed4-060a5e190c9c
md"""
## Panel A2: Average generation mix
"""

# ╔═╡ 7fe1d7c4-7c0a-4875-8c1f-72316554cb98
avg_gen = mean(x -> x.g, static_results)

# ╔═╡ 7af91ce1-3c0d-409c-b4af-6a7aee3e889d
resources = unique(df.gen.fuel_type)[1:2]  # Drop biomass

# ╔═╡ 8e4213fe-ec26-4c5f-b3a2-ae11b9bf51bd
gen_per_res = [sum(avg_gen[df.gen.fuel_type .== r]) for r in resources]

# ╔═╡ c9635b9e-5c9b-4437-b8af-3ebd27e0c85d
begin
	plta2 = plot(size = dpi .* (1.5, 1.5))
	
	bar!(resources, gen_per_res/1e3)
	
	plot!(ylabel="Generation [GW]")
end

# ╔═╡ 000fd4ee-c2a7-4924-ac34-b76f0ca80fe7
md"""
## Panel A3: Prediction error of each model
"""

# ╔═╡ e467b463-0228-4eae-9833-7b0afcb14bf3
begin
	plta3 = plot(size = dpi .* (2, 1.5))
	
	bar!(["ROD", "MDT", "RMP"], 100*[err_ed, err_mdt, err_ramp]) 
	plot!(ylabel="Relative Absolute Error", yticks=([0, 25, 50], ["0%", "25%", "50%"]))
	
end

# ╔═╡ d18bfc20-9646-4d90-a0e3-bb9bdbdff04e
md"""
## Panel A
"""

# ╔═╡ 864a14fb-61d7-4570-a1d5-0ae48af91090
begin
	l = @layout [a{0.5h}; [b{0.5w} c]]
	
	plta = plot(plta1, plta2, plta3, layout=l, size= dpi .* (4, 3))
end

# ╔═╡ 5bac312b-4a72-43fe-8dd2-2ac058b6c0f9
md"""
## Panel B1: Storage MEFs
"""

# ╔═╡ d462c03e-4c18-435a-84ec-789a21178c21
mefs_storage_full = [
	vcat(sr.r[t].dq[1]...) 
	for t in util.BLOCKS, sr in storage_results
];

# ╔═╡ 39132631-11df-4c07-a90d-fc0c962accd5
mefs_storage_total = [
	hcat([vcat(sum.(sr.r[t].dq[1])...) for t in util.BLOCKS]...)
	for sr in storage_results
];

# ╔═╡ e9044c8d-88e2-432e-958d-2219c096475f
mefs_storage_block_avg = [mean(M, dims=2)[:, 1] for M in mefs_storage_total]

# ╔═╡ a0de4556-4ccf-40b2-8629-b16d142d7f58
hpb = util.HOURS_PER_BLOCK

# ╔═╡ 04a6c4f0-9529-4f07-b788-1fc7266f665e
hpd = util.HOURS_PER_DAY

# ╔═╡ 4e80cd11-6d78-4ed2-bcce-4feace8859b7
dpb = util.HOURS_PER_BLOCK ÷ util.HOURS_PER_DAY

# ╔═╡ 6ab6dffe-d000-450d-989a-664eececb807
mefs_storage_daily_avg = [
	mean(hcat([b[(d-1)*hpd+1 : d*hpd] for d in 1:dpb]...), dims=2)[:, 1]
	for b in mefs_storage_block_avg
]

# ╔═╡ 7774426b-2380-4ca1-941f-8847da4e70ce
storage_pcs = collect(adjoint([sr.cfg[:storage_percentage] for sr in storage_results]))

# ╔═╡ 85a8b914-52e3-49ac-b35d-972745d520cb
storage_labels = [@sprintf("%.2f", pc) for pc in storage_pcs]

# ╔═╡ de47d9b0-9417-4ec8-8295-7e89c86cefd1
begin
	pltb1 = plot(size = dpi .* (4, 1.5))
	
	plot!(hcat(mefs_storage_daily_avg...), label=storage_labels, lw=2)
	
	plot!(xticks=nothing) #plot!(xlabel="Hour", xticks=6:6:24)
	plot!(ylabel="kgCO2/MWh", title="Emissions Rate")
	plot!(legend=:outertopright)
end

# ╔═╡ cd9addf7-3edc-4f97-997e-91697a42d8a6
md"""
## Panel B2: Storage total emissions
"""

# ╔═╡ 7ee21ef1-bb83-4e45-af5e-4c731928f8cb
ξs = [dd[1].q.co2 for dd in dynamic_data]

# ╔═╡ ae600073-62d2-4312-b1a1-bc4e81032f36
storage_total_emissions = [
	vcat([[sr.r[d].g[t]'ξs[d] for t in 1:hpb] for d in util.BLOCKS]...)
	for sr in storage_results
]

# ╔═╡ 5d34a1a6-700a-4769-bf5a-5628aa23136d
wpy = util.DAYS_PER_WEEK * maximum(util.WEEKS)

# ╔═╡ 36efc389-2de2-4c95-9f4c-b6aefbc36010
storage_total_emissions_daily_avg = [
	mean(hcat([v[(d-1)*hpd+1 : d*hpd] for d in 1:wpy]...), dims=2)[:, 1]
	for v in storage_total_emissions
]

# ╔═╡ 4cf4bd4e-2b94-46d2-8f86-6eead23a9798
begin
	pltb2 = plot(size = dpi .* (4, 1.5))
	
	plot!(hcat(storage_total_emissions_daily_avg...) / 1e6, 
		label=storage_labels, lw=2, alpha=0.8)
	
	plot!(xlabel="Hour", xticks=6:6:24)
	plot!(ylabel="kTCO2", title="Total Emissions")
	plot!(legend=:outertopright)
end

# ╔═╡ 6e357cae-6401-4ab9-9a88-afd1d6161507
md"""
## Panel B
"""

# ╔═╡ 0c9aedaa-3525-48eb-83cc-a282a016a165
begin
	pltb = plot(
		pltb1, pltb2, 
		size = dpi .* (4, 3),
		layout=(2, 1)
	)
end

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

# ╔═╡ 7d3094d0-eb6f-4257-a123-c4dd71fd2644
md"""
## Panel C: Heatmaps
"""

# ╔═╡ 51121966-1028-41c3-9b2d-a9ae86005076


# ╔═╡ deae71cb-3b2e-43ca-bce7-9484a93bc9fe
md"""
## Panel D: Investment (\$/mTCO2)
"""

# ╔═╡ ed1f5880-13fa-4fda-9587-8f7ab33d8143


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
# ╠═d0a30c4e-b675-439b-8bc0-adec6126360e
# ╠═98c7c6d5-8c08-42ef-9652-798f9c26137d
# ╠═63ac3134-90e6-4351-90de-54ae6c21dee1
# ╟─a359e195-6121-4438-aae8-0a3b2a1f1aa5
# ╠═118dc77a-b4d7-44c5-8eda-f9a4b29df23a
# ╠═9c78666f-f03c-4b1b-9a8c-8b784eeb0b65
# ╠═9c886bfe-7931-4cb7-8ee5-e3ef18b3de8f
# ╠═8c22b4e7-6ecf-4188-a34b-2b8e4d5f9b58
# ╠═93a455cc-c620-463e-90da-654bb96dca8b
# ╠═33d7aace-51f7-411c-bef3-79e6c672152d
# ╠═0c58311b-3001-4307-b8d3-e9365a0829c1
# ╟─3330f77a-73a3-4a8e-90b3-27824d356514
# ╠═002f2a27-8489-4521-8ed9-660ec6cd5760
# ╠═94ccf5d2-fa1f-4b08-aee5-51340bad7e56
# ╠═54db5151-a50c-43a0-aa2e-f873a1de30b1
# ╟─bdac9b46-6773-4546-8f03-4a4e1e3c4b88
# ╠═dc592388-31ab-4428-9bac-9671f12b0c5c
# ╠═e9623b71-b6de-455d-8a26-c41be6d58f8c
# ╠═9a65cc16-734c-4101-b65c-cda75aa9b75e
# ╠═74dc860b-0ff0-4cf4-90e1-0505a9d952a4
# ╠═148f5d6f-dc6d-4664-beac-6a0098422cb4
# ╠═36c864b0-abf3-48e2-8c99-48669f142896
# ╠═e3aeecd6-0a14-484b-aa28-54a360214dcb
# ╠═0502409b-2a9f-4fc5-ac75-057de60736f6
# ╠═3d83859f-391f-45cd-bf8a-99c9e1a9cc80
# ╟─1ad7f79b-b1ae-4f2a-9ed4-060a5e190c9c
# ╠═7fe1d7c4-7c0a-4875-8c1f-72316554cb98
# ╠═7af91ce1-3c0d-409c-b4af-6a7aee3e889d
# ╠═8e4213fe-ec26-4c5f-b3a2-ae11b9bf51bd
# ╠═c9635b9e-5c9b-4437-b8af-3ebd27e0c85d
# ╟─000fd4ee-c2a7-4924-ac34-b76f0ca80fe7
# ╠═e467b463-0228-4eae-9833-7b0afcb14bf3
# ╟─d18bfc20-9646-4d90-a0e3-bb9bdbdff04e
# ╠═864a14fb-61d7-4570-a1d5-0ae48af91090
# ╠═5bac312b-4a72-43fe-8dd2-2ac058b6c0f9
# ╠═d462c03e-4c18-435a-84ec-789a21178c21
# ╠═39132631-11df-4c07-a90d-fc0c962accd5
# ╠═e9044c8d-88e2-432e-958d-2219c096475f
# ╠═a0de4556-4ccf-40b2-8629-b16d142d7f58
# ╠═04a6c4f0-9529-4f07-b788-1fc7266f665e
# ╠═4e80cd11-6d78-4ed2-bcce-4feace8859b7
# ╠═6ab6dffe-d000-450d-989a-664eececb807
# ╠═7774426b-2380-4ca1-941f-8847da4e70ce
# ╠═e2a7f5a6-3a9b-44d1-b6ef-726b9870aa59
# ╠═85a8b914-52e3-49ac-b35d-972745d520cb
# ╠═de47d9b0-9417-4ec8-8295-7e89c86cefd1
# ╠═cd9addf7-3edc-4f97-997e-91697a42d8a6
# ╠═7ee21ef1-bb83-4e45-af5e-4c731928f8cb
# ╠═ae600073-62d2-4312-b1a1-bc4e81032f36
# ╠═5d34a1a6-700a-4769-bf5a-5628aa23136d
# ╠═36efc389-2de2-4c95-9f4c-b6aefbc36010
# ╠═4cf4bd4e-2b94-46d2-8f86-6eead23a9798
# ╠═6e357cae-6401-4ab9-9a88-afd1d6161507
# ╠═0c9aedaa-3525-48eb-83cc-a282a016a165
# ╠═7d3094d0-eb6f-4257-a123-c4dd71fd2644
# ╠═51121966-1028-41c3-9b2d-a9ae86005076
# ╠═deae71cb-3b2e-43ca-bce7-9484a93bc9fe
# ╠═ed1f5880-13fa-4fda-9587-8f7ab33d8143
