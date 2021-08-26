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

# ╔═╡ 819cff42-f4c4-4d63-9272-427321d4ec93
using Profile

# ╔═╡ 883eca64-7507-4c17-bf30-4c7198e39091
using Plots

# ╔═╡ caf2109d-9ccd-4c5e-8b71-4215d2bfc2b0
OPT = () -> ECOS.Optimizer(verbose=false)

# ╔═╡ de200dab-4758-4338-aca9-b2778fabd417
@bind datapath TextField((50, 1))

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
gmax_t = [df_gen[:, "mw"*t] for t in string.(weeks)]

# ╔═╡ 5a844b5b-9eb8-4743-ab37-0de66ea100d5
fl_t = [df_gen[:, "fuel_price"*t] .* df_gen[:, "heat_rate"*t] + df_gen.vom for t in string.(weeks)]

# ╔═╡ 2403bbf7-900e-45af-b576-783036591876
fq = zeros(l);

# ╔═╡ 61de6d84-8482-4d4a-ae08-86a5eb1c433b
co2_t = [df_gen[:, "co2"*t] for t in string.(weeks)]

# ╔═╡ f2f5446f-bce7-451b-8014-abaf87bf2915
so2_t = [df_gen[:, "so2"*t] for t in string.(weeks)]

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

# ╔═╡ 8613c9b7-4df7-4d23-8245-6eaeba072fac
52*7*24

# ╔═╡ b3cb2534-d21b-420f-a87a-ae63312f7050
date_format = DateFormat("yyyy-mm-dd HH:MM:SS")

# ╔═╡ 3c634cbe-b324-4915-913f-37059fc45bce
times = 1:(maximum(weeks)*7*24)

# ╔═╡ 709319da-bbea-49d7-9706-5e2be0dbafb9
d_τ = df_carbon.demand[times]

# ╔═╡ 205c588e-bc6d-4a07-b984-08ced604f5fe
datetimes = DateTime.(df_carbon.datetime[times], date_format)

# ╔═╡ a2a70891-3cb7-4c73-b9b9-42f1ee7e1954
normalized_week(dt) = ((dayofyear(dt) - 1) ÷ 7) + 1

# ╔═╡ cb7b7aea-b150-4706-a387-99fdfb3c54dc
week_τ = normalized_week.(datetimes)

# ╔═╡ ed3ea202-2ec4-40d2-afc7-4d0b65534f76
Δco2_ref_τ = df_carbon.co2_marg[times]

# ╔═╡ a388ae1d-3cbe-45e1-b339-cfbe1efa1f61
Δso2_ref_τ = df_carbon.so2_marg[times]

# ╔═╡ 345b2bf9-3b24-4048-8c3e-ee269e4ead09
md"""
## Set up dispatch problems with `CarbonNetworks`
"""

# ╔═╡ 1ab326cd-c719-4081-b7a5-db72c341ca34
Profile.clear()

# ╔═╡ dcc6aee1-152f-433d-aa09-4d1968abb0fa
Profile.print(mincount=10)

# ╔═╡ 7a3132f1-99ce-47dc-b384-e2f0f9dccbc4
net_t = [PowerNetwork(fq, fl_t[t], pmax, gmax_t[t], A, B) for t in weeks]

# ╔═╡ 2fd1e4d9-4b82-4d74-9fde-2eb75429a9c3
cases = copy(times)

# ╔═╡ c9f6b5ed-3e13-4e6c-8613-6c41f6b43926
opf_τ = [PowerManagementProblem(net_t[week_τ[τ]], d_τ[τ]) for τ in cases];

# ╔═╡ 7c3bdd70-5879-4bd6-9fad-8e3dee00d265
solve_ecos! = x -> solve!(x, OPT)

# ╔═╡ 4c97c2af-be2a-4010-b26a-5c1a66cb8246
solve_ecos!.(opf_τ);

# ╔═╡ 3fe5d64d-2a9d-441a-8568-bfa80fa3521b
get_Δco2 = τ -> begin
	w = week_τ[τ]
	return compute_mefs(opf_τ[τ], net_t[w], [d_τ[τ]], co2_t[w])[1]
end

# ╔═╡ da42dca1-257d-4624-9bcc-9885480bd808
Δco2_diff = get_Δco2.(cases)

# ╔═╡ ba1c65bd-f4c3-438f-8697-2f58a579b134
md"""
## How do the methods compare?
"""

# ╔═╡ 03270f7c-0653-4b2e-a0fe-e0945adeceae
theme(:default, label=nothing)

# ╔═╡ 200d7d2f-7979-4751-8d55-8095555ac09c
md"""
## Incorporate the minimum downtime heuristic

The details describing this procedure are listed [here](https://pubs.acs.org/doi/suppl/10.1021/acs.est.9b02500/suppl_file/es9b02500_si_001.pdf), and the code used to implement this procedure is [here](https://github.com/tdeetjen/simple_dispatch/blob/master/simple_dispatch.py#L723).
"""

# ╔═╡ b48488b7-7e8f-493b-8e29-6b755a197386
function find_mdt_events(d; Δt=12, ϵ=0.05)
	@assert iseven(Δt)
	
	T = length(d)
	times = 1:(T-Δt)
	
	∫dτ_τ = [sum(d[t:t+Δt]) for t in times]
	∫dt_τ = [d[t]*(Δt+1) for t in times]
	∫min_dtdτ_τ = [sum(min.(d[t:t+Δt], d[t])) for t in times]
	
	# Compute the intergral of the convex regions
	At = ∫dt_τ - ∫min_dtdτ_τ
	
	# Filter to ensure these are true valleys (not just actual decreases)
	At = At .* (d[times] .<= (1+ϵ)*d[times .+ Δt])
		
	# Find local maximums
	max_times = times[(Δt÷2)+2 : end-(Δt÷2)-1]
	window = (At, t, Δt) -> At[t-(Δt÷2)-1 : t+(Δt÷2)+1]
	is_max_t = [
		(At[t] == maximum(window(At, t, Δt))) && 
		(At[t] != 0) && 
		(∫dt_τ[t] >= ∫dτ_τ[t])
		for t in max_times
	]
	
	@show At[1:24]
	@show is_max_t[1:24]
	@show max_times[1:24]
	
	events = [
		(start=t, stop=t+12, demand=d[t]) 
		for (ind, t) in enumerate(max_times) if is_max_t[ind]
	]
	
	@assert all([
			!((events[k].stop - 1) in [e.start:e.stop for e in events[k+1:end]]) 
			for k in 1:length(events)
	])
	
	return events
end

# ╔═╡ aea874b7-828d-4928-8012-f3b70381d729
mdt_events = find_mdt_events(d_τ)

# ╔═╡ 42d6329a-ed86-45ad-9832-bda5dc96745a
let
	bw = 0.25
	selected = 1:48
	
	plt_mefs = plot(ylabel="Δco2 (kg / mwh)")
	bar!(cases[selected], Δco2_ref_τ[cases][selected], bar_width=bw)
	bar!(cases[selected] .+ bw, Δco2_diff[selected], bar_width=bw)
	for ev in mdt_events[1:5]
		vline!([ev.start], c=:Green, lw=4)
		vline!([ev.stop], c=:Red, lw=4)
	end
	plot!(xlim=(minimum(cases[selected])-1, maximum(cases[selected])+1))
	
	selected = 1:100
	plt_errs = plot(ylabel="error")
	bar!(cases[selected], Δco2_ref_τ[cases][selected] - Δco2_diff[selected])
	
	plot(plt_mefs, plt_errs, layout=(2, 1))
end

# ╔═╡ 662d5a0e-94a4-47f9-b276-4ca3c5c07119
let
	x = 330:400
	plot(x, d_τ[x])
	for ev in mdt_events[1:100]
		vline!([ev.start], color=:Black, lw=3)
		vline!([ev.stop], color=:Red, lw=3)
	end
	plot!(xlim=(minimum(x), maximum(x)))
	plot!()
end

# ╔═╡ a257eb36-b584-4620-bba2-9e5c4b859d48
function reshape_demand_and_capacity(mdt_events, d, opfs, is_coal, gmin_pc, c; tol=1e-1)
	times = 1:length(d)
	
	# d_new = deepcopy(d)
	gmax_new = [deepcopy(opf.params.gmax) for opf in opfs]
	fl_new = [deepcopy(opf.params.fl) for opf in opfs]
	c_new = [deepcopy(cτ) for cτ in c]
	
	is_active = t -> evaluate(opfs[t].g) .> tol
	
	for (start, stop, dt) in mdt_events
		# Find coal generators active during start
		active_coal_gens = is_coal .& is_active(start) 
				
		#@show start, stop
		#@show sum(is_coal), sum(active_coal_gens)
		
		for t in start:(stop-1)
			gmin = opfs[t].params.gmax .* gmin_pc
			g_mdt = gmin .* active_coal_gens
			
			
			
			# Create baseload generator
			#weights = active_coal_gens .* opfs[t].params.gmax	
			g_base = sum(g_mdt)
			c_base = sum(c[t] .* g_mdt) / sum(g_mdt)
			
			
			# Set active coal generators new capacity to (1-gmin_pc) 
			# of their original capacity			
			
			# And add baseload generator
			gmax_new[t] = [opfs[t].params.gmax - g_mdt; g_base]
			fl_new[t] = [opfs[t].params.fl; 0.0]
			c_new[t] = [c[t]; c_base]
				
				
			@assert all(g -> g > 0, gmax_new[t])
		end
	end
	
	return gmax_new, fl_new, c_new
end

# ╔═╡ c1b9729b-6463-485c-8b00-3116d0dbf88e
d_τ[5] - sum(df_gen.min_out)

# ╔═╡ 16b192e3-f182-4be4-b4c3-1daf974d08e7
gmax_new, fl_new, c_new = reshape_demand_and_capacity(
	mdt_events, 
	d_τ, 
	opf_τ, 
	Bool.(df_gen.is_coal), 
	df_gen.min_out_multiplier,
	[co2_t[week_τ[τ]] for τ in cases]
);

# ╔═╡ 38ba43ad-b057-4125-bbc6-baff41b7ec3e
let
	_t = 46
	_rng = 300:400
	bar((evaluate(opf_τ[_t].g) ./ opf_τ[_t].params.gmax)[_rng])
	# bar!(evaluate(opf_τ[10].g)[rng])
end

# ╔═╡ 84d2bd49-dab0-43c3-91c4-569bdc550b88
B_new = [
	(length(fl_new[τ]) == length(fl_t[1])) ? B : [B 1]
	for τ in cases
]

# ╔═╡ 05b223e4-1db4-450c-8c14-ce09a5c87f74
fq_new = [
	(length(fl_new[τ]) == length(fl_t[1])) ? fq : [fq; 0]
	for τ in cases
]

# ╔═╡ 1f6f530b-cc10-444d-831c-bc719a4a8c58
length.(B_new)

# ╔═╡ 12aedbdd-d3bc-41e9-bf1c-d3e6f6e3cf30
length.(fl_new)

# ╔═╡ 7165bd84-3d32-4127-bedf-4a4e3d82feb5
all(length.(gmax_new) .== length.(fl_new))

# ╔═╡ bac0645a-1780-4c09-912c-fc7c5b2e88d1
all(length.(gmax_new) .== length.(B_new))

# ╔═╡ 9a84321b-8d92-410e-97e8-c6d058b6b5d9
nets_new = [PowerNetwork(fq_new[τ], fl_new[τ], pmax, gmax_new[τ], A, B_new[τ]) for τ in cases]

# ╔═╡ e04b2c1b-5866-4ac3-abf7-306bc5fe47dc
opfs_new = [PowerManagementProblem(nets_new[τ], d_τ[τ]) for τ in cases];

# ╔═╡ ae394cfc-a925-436f-99ad-6caf68b551e9
bar(evaluate(opfs_new[10].g) ./ gmax_new[10])

# ╔═╡ 269d4461-cafb-4d1a-996b-0fc7a19a55a3
mdt_selected = cases

# ╔═╡ 2649fe7a-8e45-4bd3-8a28-33a718efd722
solve_ecos!.(opfs_new[mdt_selected]);

# ╔═╡ 521250c7-2f35-4846-9d53-a4db4223e76c
get_Δco2_mdt = τ -> begin
	return compute_mefs(opfs_new[τ], nets_new[τ], [d_τ[τ]], c_new[τ])[1]
end

# ╔═╡ 5f3101b0-b506-45e5-9d9f-e642caa23cfa
Δco2_diff_mdt = get_Δco2_mdt.(cases[mdt_selected])

# ╔═╡ 265fea1a-aa56-4004-a66d-3d9671aeafcc
is_below_thresh = [
	(length(fl_new[τ]) != length(fl_t[1])) && # MDT event
		(evaluate(opfs_new[τ].g)[end] < 0.99*gmax_new[τ][end])  # Baseload active
	for τ in mdt_selected
]

# ╔═╡ 6dfd8108-d3c6-4c11-a320-b241f80396fe
Δco2_heur = 0.5*Δco2_diff + [
	is_below_thresh[τ] ?
	0.5*c_new[τ][end] :
	0.5*Δco2_diff[τ]
	for τ in mdt_selected
]

# ╔═╡ 78a0f4cd-e4c6-4d60-b4e7-41147f079ad3
agreement_diff = sum(abs.(Δco2_diff - Δco2_ref_τ[cases]) .< 40) / length(cases)

# ╔═╡ 68836b97-3b90-48a5-899c-dbe03cfe262a
agreement_heur = sum(abs.(Δco2_diff_mdt - Δco2_ref_τ[cases]) .< 40) / length(cases)

# ╔═╡ ad52bf07-66e3-401e-8582-85825cab58ef
let
	cost_marg = [opf_τ[τ].problem.constraints[end].dual for τ in cases]
	cost_marg = [opfs_new[τ].problem.constraints[end].dual for τ in cases]
	
	bw = 0.2
	_x = 1:24

	
	bar(_x, cost_marg[_x], bar_width=bw)
	bar!(_x .+ bw, df_carbon.gen_cost_marg[cases][_x], bar_width=bw)
	bar!(_x .+ 2bw, df_carbon.gen_cost_marg[cases][_x], bar_width=bw)
	plot!(size=(600, 200), ylim=(28, Inf))
end

# ╔═╡ 1f00958c-beb3-4575-bf38-332b07f58285
(1:length(gmax_new[60]))[0.999 .> evaluate(opfs_new[60].g) ./ gmax_new[60] .> 0.001]

# ╔═╡ 6555a9f4-1086-43ac-87b9-2c8336306b9e
c_new[60][237]

# ╔═╡ 880e04a2-43ba-42a4-893f-a88bc0317730
co2_t[week_τ[60]][237]

# ╔═╡ 3a03867a-d66c-446e-9b33-34cb705a713a
let
	bw = 0.2
	day = 5
	selected = (12*(day-1)+1):(12*day)
	
	
	
	plt_mefs = plot(ylabel="Δco2 (kg / mwh)")
	bar!(cases[selected], Δco2_ref_τ[cases][selected], bar_width=bw)
	bar!(cases[selected] .+ bw, Δco2_diff[selected], bar_width=bw)
	bar!(cases[selected] .+ 2bw, Δco2_diff_mdt[selected], bar_width=bw)
	#bar!(cases[selected] .+ 3bw, Δco2_heur[selected], bar_width=bw)
	plot!(ylim=(300, Inf))
	
	for ev in mdt_events[1:10]
		vline!([ev.start], c=:Black, lw=4, alpha=0.5)
		vline!([ev.stop], c=:Red, lw=4)
	end
	plot!(xlim=(minimum(cases[selected])-1, maximum(cases[selected])+1))
	
	
	selected = 1:100
	plt_errs = plot(ylabel="error")
	bar!(cases[selected], Δco2_diff_mdt[selected] - Δco2_ref_τ[cases][selected])
	
	plot(plt_mefs, plt_errs, layout=(2, 1))
end

# ╔═╡ 3592b9f6-fd94-4def-b080-50850584acc0
argmax(abs.(Δco2_diff[mdt_selected] - Δco2_diff_mdt))

# ╔═╡ cc69f3c2-6fb1-43b8-bac4-a160f26a57a1
df_carbon

# ╔═╡ d4730d8e-5ee4-40c9-ae06-058d85a21c8d
names(df_gen)

# ╔═╡ 4991f1d3-8a6d-4e56-8b61-2b551df296b4
df_gen.orispl

# ╔═╡ f0e5b8fc-ea94-44e4-aeef-f66baa491b34
names(df_carbon)

# ╔═╡ 3a6f7a93-b116-4974-b992-f843b9ad7cac
begin
	co2_marg_reg = diff(df_carbon.co2_tot) ./ diff(df_carbon.demand)
	co2_marg_reg[isnan.(co2_marg_reg)] .= 0
	co2_marg_reg[co2_marg_reg .> 2000] .= 0
	co2_marg_reg[co2_marg_reg .< 0] .= 0
end

# ╔═╡ d0519350-52f1-4f75-92f2-c3b249ac6643
mef_true_mav = sum(abs, co2_marg_reg[cases]) / length(cases)

# ╔═╡ 8a8fe652-589b-4a7c-8e1a-ab52d8f772af
rodm_mad_rel = (sum(abs, Δco2_ref_τ[cases] - co2_marg_reg[cases]) / length(cases)) / mef_true_mav

# ╔═╡ c79dc171-7c98-4510-8b12-06b8be8313a5
diff_mad_rel = (sum(abs, Δco2_diff - co2_marg_reg[cases]) / length(cases)) / mef_true_mav

# ╔═╡ cd40dbc9-c1e9-40aa-be46-b93b46b5afee
diff_mdt_mad_rel = (sum(abs, Δco2_diff_mdt - co2_marg_reg[cases]) / length(cases[mdt_selected])) / mef_true_mav

# ╔═╡ 9b33484a-2c0b-4eee-b512-d46396748235
heur_mad_rel = (sum(abs, Δco2_diff_mdt[cases[mdt_selected]] - co2_marg_reg[cases[mdt_selected]]) / length(cases[mdt_selected])) / mef_true_mav

# ╔═╡ 3245d81a-4ec6-43b8-be7b-25e259a15e90
plot(co2_marg_reg)

# ╔═╡ dd38f967-71cc-41b4-a2d0-2c669c8f546b
maximum(co2_marg_reg)

# ╔═╡ 7df343e6-f0c8-4525-8538-559ea799d1e3
df_carbon

# ╔═╡ 5be2e515-067d-4e29-9613-110d4b40fca5
bar(df_gen.mw52)

# ╔═╡ Cell order:
# ╠═f6d57d0c-ff67-11eb-23fd-05b201176fd1
# ╠═20df5713-9530-4100-90ed-b949d005c71b
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
# ╠═8613c9b7-4df7-4d23-8245-6eaeba072fac
# ╠═b3cb2534-d21b-420f-a87a-ae63312f7050
# ╠═3c634cbe-b324-4915-913f-37059fc45bce
# ╠═709319da-bbea-49d7-9706-5e2be0dbafb9
# ╠═205c588e-bc6d-4a07-b984-08ced604f5fe
# ╠═a2a70891-3cb7-4c73-b9b9-42f1ee7e1954
# ╠═cb7b7aea-b150-4706-a387-99fdfb3c54dc
# ╠═ed3ea202-2ec4-40d2-afc7-4d0b65534f76
# ╠═a388ae1d-3cbe-45e1-b339-cfbe1efa1f61
# ╟─345b2bf9-3b24-4048-8c3e-ee269e4ead09
# ╠═819cff42-f4c4-4d63-9272-427321d4ec93
# ╠═1ab326cd-c719-4081-b7a5-db72c341ca34
# ╠═dcc6aee1-152f-433d-aa09-4d1968abb0fa
# ╠═7a3132f1-99ce-47dc-b384-e2f0f9dccbc4
# ╠═2fd1e4d9-4b82-4d74-9fde-2eb75429a9c3
# ╠═c9f6b5ed-3e13-4e6c-8613-6c41f6b43926
# ╠═7c3bdd70-5879-4bd6-9fad-8e3dee00d265
# ╠═4c97c2af-be2a-4010-b26a-5c1a66cb8246
# ╠═3fe5d64d-2a9d-441a-8568-bfa80fa3521b
# ╠═da42dca1-257d-4624-9bcc-9885480bd808
# ╟─ba1c65bd-f4c3-438f-8697-2f58a579b134
# ╠═883eca64-7507-4c17-bf30-4c7198e39091
# ╠═03270f7c-0653-4b2e-a0fe-e0945adeceae
# ╠═42d6329a-ed86-45ad-9832-bda5dc96745a
# ╟─200d7d2f-7979-4751-8d55-8095555ac09c
# ╠═662d5a0e-94a4-47f9-b276-4ca3c5c07119
# ╠═b48488b7-7e8f-493b-8e29-6b755a197386
# ╟─aea874b7-828d-4928-8012-f3b70381d729
# ╠═a257eb36-b584-4620-bba2-9e5c4b859d48
# ╠═c1b9729b-6463-485c-8b00-3116d0dbf88e
# ╠═16b192e3-f182-4be4-b4c3-1daf974d08e7
# ╠═38ba43ad-b057-4125-bbc6-baff41b7ec3e
# ╠═84d2bd49-dab0-43c3-91c4-569bdc550b88
# ╠═05b223e4-1db4-450c-8c14-ce09a5c87f74
# ╠═ae394cfc-a925-436f-99ad-6caf68b551e9
# ╠═1f6f530b-cc10-444d-831c-bc719a4a8c58
# ╠═12aedbdd-d3bc-41e9-bf1c-d3e6f6e3cf30
# ╠═7165bd84-3d32-4127-bedf-4a4e3d82feb5
# ╠═bac0645a-1780-4c09-912c-fc7c5b2e88d1
# ╠═9a84321b-8d92-410e-97e8-c6d058b6b5d9
# ╠═e04b2c1b-5866-4ac3-abf7-306bc5fe47dc
# ╠═269d4461-cafb-4d1a-996b-0fc7a19a55a3
# ╠═2649fe7a-8e45-4bd3-8a28-33a718efd722
# ╠═521250c7-2f35-4846-9d53-a4db4223e76c
# ╠═5f3101b0-b506-45e5-9d9f-e642caa23cfa
# ╠═d0519350-52f1-4f75-92f2-c3b249ac6643
# ╠═8a8fe652-589b-4a7c-8e1a-ab52d8f772af
# ╠═c79dc171-7c98-4510-8b12-06b8be8313a5
# ╠═cd40dbc9-c1e9-40aa-be46-b93b46b5afee
# ╠═9b33484a-2c0b-4eee-b512-d46396748235
# ╠═265fea1a-aa56-4004-a66d-3d9671aeafcc
# ╠═6dfd8108-d3c6-4c11-a320-b241f80396fe
# ╠═78a0f4cd-e4c6-4d60-b4e7-41147f079ad3
# ╠═68836b97-3b90-48a5-899c-dbe03cfe262a
# ╟─ad52bf07-66e3-401e-8582-85825cab58ef
# ╠═1f00958c-beb3-4575-bf38-332b07f58285
# ╠═6555a9f4-1086-43ac-87b9-2c8336306b9e
# ╠═880e04a2-43ba-42a4-893f-a88bc0317730
# ╠═3a03867a-d66c-446e-9b33-34cb705a713a
# ╠═3592b9f6-fd94-4def-b080-50850584acc0
# ╠═cc69f3c2-6fb1-43b8-bac4-a160f26a57a1
# ╠═d4730d8e-5ee4-40c9-ae06-058d85a21c8d
# ╠═4991f1d3-8a6d-4e56-8b61-2b551df296b4
# ╠═f0e5b8fc-ea94-44e4-aeef-f66baa491b34
# ╠═3a6f7a93-b116-4974-b992-f843b9ad7cac
# ╠═3245d81a-4ec6-43b8-be7b-25e259a15e90
# ╠═dd38f967-71cc-41b4-a2d0-2c669c8f546b
# ╠═7df343e6-f0c8-4525-8538-559ea799d1e3
# ╠═5be2e515-067d-4e29-9613-110d4b40fca5
