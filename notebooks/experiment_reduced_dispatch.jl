### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ f6d57d0c-ff67-11eb-23fd-05b201176fd1
begin
	using Pkg; Pkg.activate()
	using CSV
	using DataFrames
	using PlutoUI
	using Dates
	using SparseArrays
	using Convex, ECOS
end

# ╔═╡ 20df5713-9530-4100-90ed-b949d005c71b
begin
	using Revise
	using CarbonNetworks
end

# ╔═╡ 883eca64-7507-4c17-bf30-4c7198e39091
using Plots

# ╔═╡ 5c4ef9bc-2f95-400f-b4ce-45796369eff7
using Statistics

# ╔═╡ 03270f7c-0653-4b2e-a0fe-e0945adeceae
theme(:default, label=nothing)

# ╔═╡ caf2109d-9ccd-4c5e-8b71-4215d2bfc2b0
OPT = () -> ECOS.Optimizer(verbose=false)

# ╔═╡ de200dab-4758-4338-aca9-b2778fabd417
@bind datapath TextField((50, 1); default="~/Data/carbon_networks/")

# ╔═╡ 1eaa75e5-6664-4abf-8aba-42ccdd053492
datapath

# ╔═╡ fadc80fd-8917-4ef4-ae00-170d669c9cb1
PATH_CARBON = "NoEVs_year2019_dpdf_20210816.csv"

# ╔═╡ 669602ef-85e8-44dd-b306-62656dbf0aea
PATH_GENERATOR = "NoEVs_year2019_bsdf_20210816.csv"

# ╔═╡ 08aa41ac-5bb1-4113-8080-5802999ca3b4
df_carbon = DataFrame(CSV.File(expanduser(joinpath(datapath, PATH_CARBON))))

# ╔═╡ 70dfc33c-566c-417f-b09e-d529ac8a3dc1
begin
	df_gen = DataFrame(CSV.File(expanduser(joinpath(datapath, PATH_GENERATOR))))
	dropmissing!(df_gen)
	df_gen = df_gen[3:end, :]  # <-- First two rows are all zero...
end

# ╔═╡ 0672086d-ea46-49c5-8bdc-e0c0abe92a94
md"""
## Construct network

We need to create
- The `n` nodes of the network
- The `l` generators of the network
- The node-generator map `B`
- The generation capacities `gmax_t` for each week
- The generator prices `fl_t` (we assume `fq = 0`) for each week
- The generator carbon emission rates for each week `co2_t, so2_t`
"""

# ╔═╡ c56a6e96-5483-4619-99dc-f47db24a8545
node_labels = unique(df_gen.nerc)

# ╔═╡ f2b1f653-1ad0-4d81-a17a-3d56c7edbbc7
n = length(node_labels)

# ╔═╡ e9ed83b1-779e-4381-830d-c269e6fca29a
node_label_ids = Dict(label => i for (i, label) in enumerate(node_labels))

# ╔═╡ 097bab7d-5e94-43a8-93e1-5b693cc42f6d
node_label_to_id = label -> node_label_ids[label]

# ╔═╡ da0bf30f-1869-4415-a907-5141effe0134
nrow(df_gen)

# ╔═╡ b64137cf-3c29-472b-b2b2-f01705f4c9ce
l = length(df_gen.nerc)

# ╔═╡ 6f1609c5-e35d-4d29-8f6b-0d49f7485c68
generator_nodes = node_label_to_id.(df_gen.nerc)

# ╔═╡ f96577db-ad64-4b3d-8e3b-79d647d6dd7b
begin
	B = zeros(n, l)
	for (gen_ind, node_ind) in enumerate(generator_nodes)
		B[node_ind, gen_ind] = 1
	end
	B
end

# ╔═╡ de89fae0-0d64-432e-aa2b-685a533225c5
weeks = 1:52

# ╔═╡ 32144008-61a4-4ce8-8fd0-be131912427f
gmax_week = [df_gen[:, "mw"*t] for t in string.(weeks)]

# ╔═╡ 5a844b5b-9eb8-4743-ab37-0de66ea100d5
fl_week = [df_gen[:, "fuel_price"*t] .* df_gen[:, "heat_rate"*t] + df_gen.vom for t in string.(weeks)]

# ╔═╡ 2403bbf7-900e-45af-b576-783036591876
fq = zeros(l);

# ╔═╡ 61de6d84-8482-4d4a-ae08-86a5eb1c433b
co2_week = [df_gen[:, "co2"*t] for t in string.(weeks)]

# ╔═╡ f2f5446f-bce7-451b-8014-abaf87bf2915
so2_week = [df_gen[:, "so2"*t] for t in string.(weeks)]

# ╔═╡ 558b4e6a-0f4b-474e-a84b-fed466d8be4c
A = zeros(n, 1)

# ╔═╡ f37ab50c-2307-4baf-8952-61b2dcc4453c
m = size(A, 2)

# ╔═╡ 65a03949-8950-4874-b944-d1eb3f11804b
pmax = 0.1 * ones(m)

# ╔═╡ fe8fe8c7-2a51-4515-a7e3-f67ed3abdf62
md"""
## Process carbon data

The demand data is at the hourly resolution (which we denote with `τ`), in contrast with the generator data, which changes weekly (which we denote by `t`).

We need to create
- A demand time series `d_τ`
- The week of each demand time `week_τ`
- The calculated marginal emission rates from (Deetjen, Azevedo 2019) `Δco2_ref_τ, Δso2_ref_τ`
"""

# ╔═╡ b3cb2534-d21b-420f-a87a-ae63312f7050
date_format = DateFormat("yyyy-mm-dd HH:MM:SS")

# ╔═╡ 3c634cbe-b324-4915-913f-37059fc45bce
times = 1:(maximum(weeks)*7*24)

# ╔═╡ 709319da-bbea-49d7-9706-5e2be0dbafb9
demand = df_carbon.demand[times]

# ╔═╡ 205c588e-bc6d-4a07-b984-08ced604f5fe
datetimes = DateTime.(df_carbon.datetime[times], date_format)

# ╔═╡ a2a70891-3cb7-4c73-b9b9-42f1ee7e1954
normalized_week(dt) = ((dayofyear(dt) - 1) ÷ 7) + 1

# ╔═╡ cb7b7aea-b150-4706-a387-99fdfb3c54dc
week_time = normalized_week.(datetimes)

# ╔═╡ ed3ea202-2ec4-40d2-afc7-4d0b65534f76
Δco2_ref = df_carbon.co2_marg[times]

# ╔═╡ a388ae1d-3cbe-45e1-b339-cfbe1efa1f61
Δso2_ref = df_carbon.so2_marg[times]

# ╔═╡ 345b2bf9-3b24-4048-8c3e-ee269e4ead09
md"""
## Set up dispatch problems with `CarbonNetworks`
"""

# ╔═╡ 7a3132f1-99ce-47dc-b384-e2f0f9dccbc4
net_week = [PowerNetwork(fq, fl_week[w], pmax, gmax_week[w], A, B) for w in weeks]

# ╔═╡ c9f6b5ed-3e13-4e6c-8613-6c41f6b43926
opfs = [PowerManagementProblem(net_week[week_time[t]], demand[t]) for t in times];

# ╔═╡ 7c3bdd70-5879-4bd6-9fad-8e3dee00d265
solve_ecos! = x -> solve!(x, OPT)

# ╔═╡ 4c97c2af-be2a-4010-b26a-5c1a66cb8246
solve_ecos!.(opfs);

# ╔═╡ 3fe5d64d-2a9d-441a-8568-bfa80fa3521b
get_Δco2 = τ -> begin
	w = week_time[τ]
	return compute_mefs(opfs[τ], net_week[w], [demand[τ]], co2_week[w])[1]
end

# ╔═╡ 01e04ca5-03e4-4deb-9ae9-964ef3642bb6
[co2_week[1] so2_week[1]]

# ╔═╡ b74cd811-9880-4685-8e9f-35db4b1c603e
compute_mefs(opfs[1], net_week[1], [demand[1]], [co2_week[1] so2_week[1]])

# ╔═╡ da42dca1-257d-4624-9bcc-9885480bd808
Δco2_diff = get_Δco2.(times)

# ╔═╡ ba1c65bd-f4c3-438f-8697-2f58a579b134
md"""
## How do the methods compare?
"""

# ╔═╡ 798df6bc-553f-454c-878f-a007d8b93928
md"""
### First, to one another
"""

# ╔═╡ 42d6329a-ed86-45ad-9832-bda5dc96745a
let
	bw = 0.25
	selected = 1:24
	
	plt_mefs = plot(ylabel="Δco2 (kg / mwh)")
	bar!(selected, Δco2_ref[selected], bar_width=bw)
	bar!(selected .+ bw, Δco2_diff[selected], bar_width=bw)
	
	selected = 1:50
	plt_errs = plot(ylabel="error")
	bar!(selected, (Δco2_diff - Δco2_ref)[selected])
	
	plot(plt_mefs, plt_errs, layout=(2, 1))
end

# ╔═╡ ebe751b0-04f7-4f64-ab00-3454d28fa204
md"""
### Second, to the true data
"""

# ╔═╡ 3a6f7a93-b116-4974-b992-f843b9ad7cac
begin
	Δco2 = (diff(df_carbon.co2_tot) ./ diff(df_carbon.demand))[times]
	Δco2[isnan.(Δco2)] .= -1
	Δco2[Δco2 .> 2000] .= -1
	Δco2[Δco2 .< 0] .= -1
	
	Δco2
end

# ╔═╡ 735e95e9-6a3d-41ef-813c-b83a0cb2ab7d
is_valid = (Δco2 .!= -1)

# ╔═╡ 9b40aec8-fecd-4089-99b1-c61dcb74d774
get_err(Δco2_est) = mean(abs, (Δco2_est - Δco2)[is_valid]) / mean(Δco2[is_valid])

# ╔═╡ 1040aa20-dff5-453f-991a-f00ce9098b5f
err_ref = get_err(Δco2_ref)

# ╔═╡ 2f8194e1-73c5-4342-9b18-613df39f91f6
err_diff = get_err(Δco2_diff)

# ╔═╡ 200d7d2f-7979-4751-8d55-8095555ac09c
md"""
## Incorporate the minimum downtime heuristic

The details describing this procedure are listed [here](https://pubs.acs.org/doi/suppl/10.1021/acs.est.9b02500/suppl_file/es9b02500_si_001.pdf), and the code used to implement this procedure is [here](https://github.com/tdeetjen/simple_dispatch/blob/master/simple_dispatch.py#L723).


**Update.** We were unable to implement this heuristic as implemented in DA19. As a result, we have just chosen to ignore it, and to compute errors between DA19's results and the proposed approach using hours outside the MDT events.
"""

# ╔═╡ b48488b7-7e8f-493b-8e29-6b755a197386
function find_mdt_events(d; Δt=12, ϵ=0.05)
	@assert iseven(Δt)
	
	T = length(d)
	times = 1:(T-Δt)
	
	∫dτ_τ = [sum(d[t:t+Δt]) for t in times]
	∫dt_τ = d[times]*(Δt+1)
	∫min_dtdτ_τ = [sum(min.(d[t:t+Δt], d[t])) for t in times]
	
	# Compute the intergral of the convex regions
	At = ∫dt_τ - ∫min_dtdτ_τ
	
	# Filter to ensure these are true valleys (not just actual decreases)
	At = At .* (d[times] .<= (1+ϵ)*d[times .+ Δt])
		
	# Find local maximums
	max_times = times[(Δt÷2)+1 : end-(Δt÷2)]
	window = (At, t, Δt) -> At[t-(Δt÷2) : t+(Δt÷2)]
	
	is_max_t = [
		(At[t] == maximum(window(At, t, Δt))) && 
		(At[t] != 0) #&& 
		#(∫dt_τ[t] >= ∫dτ_τ[t])
		for t in max_times
	]
	
	#@show 309 .+ findall(x->x==1, is_max_t[310:340])
	#@show max_times[310]
	@show [x for x in zip(max_times[300:330], At[max_times[300:330]], is_max_t[300:330], ∫dt_τ[max_times[300:330]], ∫dτ_τ[max_times[300:330]])]
	
	events = [
		(start=t, stop=t+12, demand=d[t]) 
		for (ind, t) in enumerate(max_times) if is_max_t[ind]
	]
	
	# @assert all([
	# 		!((events[k].stop - 1) in [e.start:e.stop for e in events[k+1:end]]) 
	# 		for k in 1:length(events)
	# ])
	
	return events
end

# ╔═╡ aea874b7-828d-4928-8012-f3b70381d729
mdt_events = find_mdt_events(demand)

# ╔═╡ 02f2017a-53c6-482b-8feb-78997d263527
is_in_mdt_event(τ) = any([start <= τ <= stop for (start, stop, _) in mdt_events])

# ╔═╡ 01c95de8-8dec-4e94-9a50-028e1f387aee
not_mdt_times = .! is_in_mdt_event.(times)

# ╔═╡ caa4ed02-da02-494d-a985-e3dc57a48846
num_non_mdt_times = sum(not_mdt_times)

# ╔═╡ 609bb8c5-0ea7-41eb-a562-e3ca0e1948a6
num_non_mdt_errs = sum(
	e -> abs(e) > 0.1, 
	(Δco2_diff - Δco2_ref)[not_mdt_times]
)

# ╔═╡ 421dcc77-888a-49b0-b9f8-d815f8f4746f
"Percent error: $(num_non_mdt_errs / num_non_mdt_times)"

# ╔═╡ 662d5a0e-94a4-47f9-b276-4ca3c5c07119
let
	x = 310:340
	plot(x, demand[x])
	for ev in mdt_events[1:100]
		vline!([ev.start], color=:Green, lw=3, ls=:dash)
		vline!([ev.stop], color=:Red, lw=3, ls=:dot)
	end
	plot!(xlim=(minimum(x), maximum(x)))
	plot!(size=(600, 200))
end

# ╔═╡ d18ee136-f8d2-42a5-993a-d14d98cf6a4c
let
	x = times
	bar(x, ((Δco2_diff - Δco2_ref) .* not_mdt_times)[x], size=(600, 200))
end

# ╔═╡ a257eb36-b584-4620-bba2-9e5c4b859d48
# function reshape_demand_and_capacity(mdt_events, d, opfs, is_coal, gmin_pc, c; tol=1e-1)
# 	times = 1:length(d)
	
# 	# d_new = deepcopy(d)
# 	gmax_new = [deepcopy(opf.params.gmax) for opf in opfs]
# 	fl_new = [deepcopy(opf.params.fl) for opf in opfs]
# 	c_new = [deepcopy(cτ) for cτ in c]
	
# 	is_active = t -> evaluate(opfs[t].g) .> tol
	
# 	for (start, stop, dt) in mdt_events
# 		# Find coal generators active during start
# 		active_coal_gens = is_coal .& is_active(start) 
				
# 		#@show start, stop
# 		#@show sum(is_coal), sum(active_coal_gens)
		
# 		for t in start:(stop-1)
# 			gmin = opfs[t].params.gmax .* gmin_pc
# 			g_mdt = gmin .* active_coal_gens
			
			
			
# 			# Create baseload generator
# 			#weights = active_coal_gens .* opfs[t].params.gmax	
# 			g_base = sum(g_mdt)
# 			c_base = sum(c[t] .* g_mdt) / sum(g_mdt)
			
			
# 			# Set active coal generators new capacity to (1-gmin_pc) 
# 			# of their original capacity			
			
# 			# And add baseload generator
# 			gmax_new[t] = [opfs[t].params.gmax - g_mdt; g_base]
# 			fl_new[t] = [opfs[t].params.fl; 0.0]
# 			c_new[t] = [c[t]; c_base]
				
				
# 			@assert all(g -> g > 0, gmax_new[t])
# 		end
# 	end
	
# 	return gmax_new, fl_new, c_new
# end
"Disabled"

# ╔═╡ 16b192e3-f182-4be4-b4c3-1daf974d08e7
# gmax_new, fl_new, c_new = reshape_demand_and_capacity(
# 	mdt_events, 
# 	d_τ[cases], 
# 	opf_τ, 
# 	Bool.(df_gen.is_coal), 
# 	df_gen.min_out_multiplier,
# 	[co2_t[week_τ[τ]] for τ in cases]
# );
"Disabled"

# ╔═╡ 84d2bd49-dab0-43c3-91c4-569bdc550b88
# B_new = [
# 	(length(fl_new[τ]) == length(fl_t[1])) ? B : [B 1]
# 	for τ in cases
# ]
"Disabled"

# ╔═╡ 05b223e4-1db4-450c-8c14-ce09a5c87f74
# fq_new = [
# 	(length(fl_new[τ]) == length(fl_t[1])) ? fq : [fq; 0]
# 	for τ in cases
# ]
"Disabled"

# ╔═╡ 9a84321b-8d92-410e-97e8-c6d058b6b5d9
"Disabled"
# nets_new = [PowerNetwork(fq_new[τ], fl_new[τ], pmax, gmax_new[τ], A, B_new[τ]) for τ in cases]

# ╔═╡ e04b2c1b-5866-4ac3-abf7-306bc5fe47dc
# opfs_new = [PowerManagementProblem(nets_new[τ], d_τ[τ]) for τ in cases];
"Disabled"

# ╔═╡ 2649fe7a-8e45-4bd3-8a28-33a718efd722
# solve_ecos!.(opfs_new[mdt_selected]);
"Disabled"

# ╔═╡ 521250c7-2f35-4846-9d53-a4db4223e76c
# get_Δco2_mdt = τ -> begin
# 	return compute_mefs(opfs_new[τ], nets_new[τ], [d_τ[τ]], c_new[τ])[1]
# end
"Disabled"

# ╔═╡ 5f3101b0-b506-45e5-9d9f-e642caa23cfa
# Δco2_diff_mdt = get_Δco2_mdt.(cases[mdt_selected])
"Disabled"

# ╔═╡ e83dbdbd-0739-491f-96a5-8feb51ff3ebe
md"""
## Part 2: Incorporating storage
"""

# ╔═╡ d4e4109d-2113-4607-9d78-6a37f6980e37
md"""
#### Add a carbon tax
"""

# ╔═╡ 6b4d9b5c-eef9-4e39-b141-216ac5a642a1
carbon_tax = 0.03  # cents / kg

# ╔═╡ 0fed2542-b66f-40e7-9c02-afd289aa742f
fl_tax = [fl_week[w] + carbon_tax * co2_week[w] for w in weeks]

# ╔═╡ 76162299-5e8e-46ea-ad4f-d5e83fe1098f
plot(demand[1:24])

# ╔═╡ e77aee2e-dd93-435d-8add-badc926f54aa
let
	w = 1
	selected = 1:550
	
	perma = sortperm(fl_week[w])[selected]
	pa1 = bar(fl_week[w][perma], lw=0, title="current", ylabel="\$/mwh")
	pa2 = bar(co2_week[w][perma], lw=0, ylabel="kg co2 / mwh", ylim=(0, 1200))
	
	permb = sortperm(fl_tax[w])[selected]
	pb1 = bar(fl_tax[w][permb], lw=0, title="carbon tax")
	pb2 = bar(co2_week[w][permb], lw=0, ylim=(0, 1200))
	
	plot(pa1, pb1, pa2, pb2, layout=(2, 2))
end

# ╔═╡ aef64b05-737d-4607-80c4-83b7e0917ad2
md"""
#### Make dynamic problems
"""

# ╔═╡ ae6886a3-c1ee-462c-99bf-8d8f193f9ef4
c_rate = 0.25

# ╔═╡ 4718e5e9-b04c-4bde-9216-6b955b33f23d
storage_pc = 0.05

# ╔═╡ 764fc0b6-9aee-49dc-8a0b-c436b36e226d
hpd = 24

# ╔═╡ a4870995-b696-45d6-8795-7e9ed96bdcbb
days = 1 : hpd : (times[end] - hpd)

# ╔═╡ b4e6ec61-a1e0-4628-b2d3-90d8ef09684d
TEST_DAYS = 1:length(days)

# ╔═╡ 5f152c7a-9e3a-42d1-96fd-ea0b25399e1c
net_day = [
	DynamicPowerNetwork(
		[zeros(l) for h in 1:hpd],
		[fl_tax[week_time[td+h-1]] for h in 1:hpd],
		[pmax for h in 1:hpd],
		[gmax_week[week_time[td+h-1]] for h in 1:hpd],
		A,
		B,
		[c_rate * storage_pc * sum(demand[td:td+hpd-1])],
		[storage_pc * sum(demand[td:td+hpd-1])],
		hpd,
		η_c=0.8
	)
	for td in days
];

# ╔═╡ 548e3481-80ad-4248-9823-650254f5218f
co2_day = co2_week[week_time[days]];

# ╔═╡ 903995a5-143c-40b4-a9bf-7a999775cbfd
begin
	dyn_opf_day = [
		DynamicPowerManagementProblem(
			net_day[td_ind], 
			[[demand[h]] for h in td:td+hpd-1]
		)
		for (td_ind, td) in enumerate(days)
	]
	solve_ecos!.(dyn_opf_day[TEST_DAYS])
end;

# ╔═╡ eee98608-292a-4828-b647-b48d8792a8ab
function reduce_problem(opf, net, co2; tol=1e-4)
	gmax = net.gmax
	
	T = length(net.gmax)
	l = length(net.gmax[1])
	
	# Find generators that are empty and at capacity at time t
	rel_gen = t -> evaluate(opf.g[t]) ./ gmax[t]
	empty_gens = [rel_gen(t) .< tol for t in 1:T]
	capped_gens = [rel_gen(t) .> (1-tol) for t in 1:T]
	
	# Find generators that are always active / always at capacity
	never_active = [prod(x -> x[i], empty_gens) for i in 1:l]
	always_capped = [prod(x -> x[i], capped_gens) for i in 1:l]
	
	# Keep relevant generators (that aren't always bound to the same constraint)
	relevant_gens = .! (always_capped .| never_active)
	l_red = sum(relevant_gens)
	
	net_red = DynamicPowerNetwork(
		[net.fq[t][relevant_gens] for t in 1:T],
		[net.fl[t][relevant_gens] for t in 1:T],
		[net.pmax[t] for t in 1:T],
		[net.gmax[t][relevant_gens] for t in 1:T],
		A,
		B[:, relevant_gens],
		net.C,
		net.P,
		T,
		η_c=0.8
	)

	demand_red = [
		dt .- sum(net.gmax[t][always_capped]) 
		for (t, dt) in enumerate(opf.params.d)
	]
	
	
	opf_red = DynamicPowerManagementProblem(net_red, demand_red)
	solve_ecos!(opf_red)
	
	return opf_red, net_red, co2[relevant_gens]
end

# ╔═╡ 480b023a-3aac-4dba-8157-663cb1eae539
dyn_opf_red, net_red, co2_red = zip(reduce_problem.(
	dyn_opf_day[TEST_DAYS], net_day[TEST_DAYS], co2_day[TEST_DAYS]
)...);

# ╔═╡ 2a47df56-d56b-4ea6-b9b9-fb501b2f8a30
get_Δco2_dyn = day -> begin
	println(day)
	return compute_mefs(
		dyn_opf_red[day], 
		net_red[day], 
		dyn_opf_red[day].params.d, 
		co2_red[day]
	)
end

# ╔═╡ 698b6698-3837-4305-b16e-6af70b26ed1b
Δco2_dyn = get_Δco2_dyn.(TEST_DAYS)

# ╔═╡ 756b5a25-f661-402d-a0a1-ed9b1b25aea5
plt_day = dayofyear(Date(2020, 7, 2))

# ╔═╡ 491585a4-0f5d-4d78-9e7f-aacc68b72123
let
	d = plt_day
	
	
	plot(vcat([sum(Δco2_dyn[dp])' for dp in d:(d+6)]...), lw=4, size=(600, 200), label="storage")
	plot!(Δco2_diff[days[d]:days[d+6]+hpd-1], lw=4, ylabel="mef")
end

# ╔═╡ ba628523-7b90-4fac-bf6a-f5d637ed7317
let
	d = plt_day
	
	
	heatmap(abs.(vcat(Δco2_dyn[d]...)), c=:tempo)
end

# ╔═╡ Cell order:
# ╠═f6d57d0c-ff67-11eb-23fd-05b201176fd1
# ╠═20df5713-9530-4100-90ed-b949d005c71b
# ╠═883eca64-7507-4c17-bf30-4c7198e39091
# ╠═03270f7c-0653-4b2e-a0fe-e0945adeceae
# ╠═caf2109d-9ccd-4c5e-8b71-4215d2bfc2b0
# ╠═de200dab-4758-4338-aca9-b2778fabd417
# ╠═1eaa75e5-6664-4abf-8aba-42ccdd053492
# ╠═fadc80fd-8917-4ef4-ae00-170d669c9cb1
# ╠═669602ef-85e8-44dd-b306-62656dbf0aea
# ╠═08aa41ac-5bb1-4113-8080-5802999ca3b4
# ╠═70dfc33c-566c-417f-b09e-d529ac8a3dc1
# ╟─0672086d-ea46-49c5-8bdc-e0c0abe92a94
# ╠═f2b1f653-1ad0-4d81-a17a-3d56c7edbbc7
# ╠═c56a6e96-5483-4619-99dc-f47db24a8545
# ╠═e9ed83b1-779e-4381-830d-c269e6fca29a
# ╠═097bab7d-5e94-43a8-93e1-5b693cc42f6d
# ╠═da0bf30f-1869-4415-a907-5141effe0134
# ╠═b64137cf-3c29-472b-b2b2-f01705f4c9ce
# ╠═6f1609c5-e35d-4d29-8f6b-0d49f7485c68
# ╠═f96577db-ad64-4b3d-8e3b-79d647d6dd7b
# ╠═de89fae0-0d64-432e-aa2b-685a533225c5
# ╠═32144008-61a4-4ce8-8fd0-be131912427f
# ╠═5a844b5b-9eb8-4743-ab37-0de66ea100d5
# ╠═2403bbf7-900e-45af-b576-783036591876
# ╠═61de6d84-8482-4d4a-ae08-86a5eb1c433b
# ╠═f2f5446f-bce7-451b-8014-abaf87bf2915
# ╠═558b4e6a-0f4b-474e-a84b-fed466d8be4c
# ╠═f37ab50c-2307-4baf-8952-61b2dcc4453c
# ╠═65a03949-8950-4874-b944-d1eb3f11804b
# ╟─fe8fe8c7-2a51-4515-a7e3-f67ed3abdf62
# ╠═b3cb2534-d21b-420f-a87a-ae63312f7050
# ╠═3c634cbe-b324-4915-913f-37059fc45bce
# ╠═709319da-bbea-49d7-9706-5e2be0dbafb9
# ╠═205c588e-bc6d-4a07-b984-08ced604f5fe
# ╠═a2a70891-3cb7-4c73-b9b9-42f1ee7e1954
# ╠═cb7b7aea-b150-4706-a387-99fdfb3c54dc
# ╠═ed3ea202-2ec4-40d2-afc7-4d0b65534f76
# ╠═a388ae1d-3cbe-45e1-b339-cfbe1efa1f61
# ╟─345b2bf9-3b24-4048-8c3e-ee269e4ead09
# ╠═7a3132f1-99ce-47dc-b384-e2f0f9dccbc4
# ╠═c9f6b5ed-3e13-4e6c-8613-6c41f6b43926
# ╠═7c3bdd70-5879-4bd6-9fad-8e3dee00d265
# ╠═4c97c2af-be2a-4010-b26a-5c1a66cb8246
# ╠═3fe5d64d-2a9d-441a-8568-bfa80fa3521b
# ╠═01e04ca5-03e4-4deb-9ae9-964ef3642bb6
# ╠═b74cd811-9880-4685-8e9f-35db4b1c603e
# ╠═da42dca1-257d-4624-9bcc-9885480bd808
# ╟─ba1c65bd-f4c3-438f-8697-2f58a579b134
# ╟─798df6bc-553f-454c-878f-a007d8b93928
# ╠═42d6329a-ed86-45ad-9832-bda5dc96745a
# ╟─ebe751b0-04f7-4f64-ab00-3454d28fa204
# ╠═3a6f7a93-b116-4974-b992-f843b9ad7cac
# ╠═735e95e9-6a3d-41ef-813c-b83a0cb2ab7d
# ╠═5c4ef9bc-2f95-400f-b4ce-45796369eff7
# ╠═9b40aec8-fecd-4089-99b1-c61dcb74d774
# ╠═1040aa20-dff5-453f-991a-f00ce9098b5f
# ╠═2f8194e1-73c5-4342-9b18-613df39f91f6
# ╟─200d7d2f-7979-4751-8d55-8095555ac09c
# ╠═b48488b7-7e8f-493b-8e29-6b755a197386
# ╠═aea874b7-828d-4928-8012-f3b70381d729
# ╠═02f2017a-53c6-482b-8feb-78997d263527
# ╠═01c95de8-8dec-4e94-9a50-028e1f387aee
# ╠═caa4ed02-da02-494d-a985-e3dc57a48846
# ╠═609bb8c5-0ea7-41eb-a562-e3ca0e1948a6
# ╟─421dcc77-888a-49b0-b9f8-d815f8f4746f
# ╠═662d5a0e-94a4-47f9-b276-4ca3c5c07119
# ╠═d18ee136-f8d2-42a5-993a-d14d98cf6a4c
# ╟─a257eb36-b584-4620-bba2-9e5c4b859d48
# ╟─16b192e3-f182-4be4-b4c3-1daf974d08e7
# ╟─84d2bd49-dab0-43c3-91c4-569bdc550b88
# ╟─05b223e4-1db4-450c-8c14-ce09a5c87f74
# ╟─9a84321b-8d92-410e-97e8-c6d058b6b5d9
# ╟─e04b2c1b-5866-4ac3-abf7-306bc5fe47dc
# ╟─2649fe7a-8e45-4bd3-8a28-33a718efd722
# ╟─521250c7-2f35-4846-9d53-a4db4223e76c
# ╟─5f3101b0-b506-45e5-9d9f-e642caa23cfa
# ╟─e83dbdbd-0739-491f-96a5-8feb51ff3ebe
# ╟─d4e4109d-2113-4607-9d78-6a37f6980e37
# ╠═6b4d9b5c-eef9-4e39-b141-216ac5a642a1
# ╠═0fed2542-b66f-40e7-9c02-afd289aa742f
# ╠═76162299-5e8e-46ea-ad4f-d5e83fe1098f
# ╠═e77aee2e-dd93-435d-8add-badc926f54aa
# ╟─aef64b05-737d-4607-80c4-83b7e0917ad2
# ╠═ae6886a3-c1ee-462c-99bf-8d8f193f9ef4
# ╠═4718e5e9-b04c-4bde-9216-6b955b33f23d
# ╠═764fc0b6-9aee-49dc-8a0b-c436b36e226d
# ╠═a4870995-b696-45d6-8795-7e9ed96bdcbb
# ╠═b4e6ec61-a1e0-4628-b2d3-90d8ef09684d
# ╠═5f152c7a-9e3a-42d1-96fd-ea0b25399e1c
# ╠═548e3481-80ad-4248-9823-650254f5218f
# ╠═903995a5-143c-40b4-a9bf-7a999775cbfd
# ╠═eee98608-292a-4828-b647-b48d8792a8ab
# ╠═480b023a-3aac-4dba-8157-663cb1eae539
# ╠═2a47df56-d56b-4ea6-b9b9-fb501b2f8a30
# ╠═698b6698-3837-4305-b16e-6af70b26ed1b
# ╠═756b5a25-f661-402d-a0a1-ed9b1b25aea5
# ╠═491585a4-0f5d-4d78-9e7f-aacc68b72123
# ╠═ba628523-7b90-4fac-bf6a-f5d637ed7317
