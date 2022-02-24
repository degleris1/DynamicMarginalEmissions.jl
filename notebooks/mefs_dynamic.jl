### A Pluto.jl notebook ###
# v0.16.1

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

# ╔═╡ 26a4178a-de77-11eb-1fa3-5fc861d32904
begin
	using Pkg; Pkg.activate()
	using Distributions, LightGraphs, Random
	using SparseArrays
	using Convex, ECOS
	using Plots
	using PlutoUI
	
	using Revise
	using CarbonNetworks
	
	OPT = () -> ECOS.Optimizer(verbose=false)
	theme(:default, label=nothing, tickfont=(:Times, 8), guidefont=(:Times, 8))
end;

# ╔═╡ c8406984-bd2c-40fd-92b9-397e07d9f213
md"
## Generate data
"

# ╔═╡ ad0d7c81-c88e-4d21-84e0-1b8e8f44bf5c
md"
### Network
"

# ╔═╡ b3ce8622-b8dd-4f9b-8593-4bd21468285d
begin
	Random.seed!(6)
	
	n = 10
	l = 30

	# Make graph
	G = watts_strogatz(n, 3, 0.2)

	# Convert to incidence matrix
	A = incidence_matrix(G, oriented=true)
	m = size(A, 2)

	# Create generator-node mapping
	node_map = vcat(collect(1:n), rand(1:n, l-n))
	B = spzeros(n, l)
	for (gen, node) in enumerate(node_map)
		B[node, gen] = 1
	end

	# Generate carbon costs
	cq = zeros(l)
	cl = rand(Exponential(2), l)
	
	# Generate regular costs
	fq = rand(Exponential(2), l) .+ 2
	fl = rand(Exponential(2), l) .+ 2

	# Generate generation and flow capacities
	gmax = rand(Gamma(4.0, 3.0), l)
	pmax = rand(Gamma(1.0, 0.1), m)
end;

# ╔═╡ 99af5f08-2ad3-4ee0-840f-474e20a05ea3
md"
### Generate cases
"

# ╔═╡ edd937f1-c5f3-45dc-94e9-84234858db2f
T = 12;

# ╔═╡ bcb8825f-afb0-4f5b-a25d-4891579f9d39
ds = [rand(Uniform(0.25, 1.0), n) .* (B*gmax) for t in 1:T];

# ╔═╡ eea80101-05d7-4f01-bab8-b89d0aef7911
md"
## Solve static problem
"

# ╔═╡ ea6cb153-e45c-480b-a36d-4cf7711def62
net = PowerNetwork(fq, fl, pmax, gmax, A, B);

# ╔═╡ f4ab0602-3954-4a5d-b770-3d9e455b500a
begin
	opf_static = PowerManagementProblem(net, ds[1])
	solve!(opf_static, OPT, verbose=true)
end;

# ╔═╡ c6e83723-ab3d-41e0-b323-07283ce72281
begin
	g_static = evaluate(opf_static.g);
	bar(g_static, size=(800, 200))
end

# ╔═╡ ae595fbf-f501-49e1-93f4-217646c29396
md"
## Solve dynamic problem - no storage
"

# ╔═╡ 88490287-4f63-4c45-bf58-18df9fe711b6
begin
	fqs = [fq for _ in 1:T]
	fls = [fl for _ in 1:T]
	pmaxs = [pmax for _ in 1:T]
	gmaxs = [gmax for _ in 1:T]
end;

# ╔═╡ af67cc5a-ca58-448b-912a-f438172deba8
net_no_cap = DynamicPowerNetwork(fqs, fls, pmaxs, gmaxs, A, B, zeros(n), zeros(n));

# ╔═╡ 7321565d-9eab-4871-8636-4bc6a57bc10f
opf_no_cap = DynamicPowerManagementProblem(net_no_cap, ds);

# ╔═╡ dbb14945-5fea-479d-be31-e9c14e1b8672
solve!(opf_no_cap, OPT, verbose=true)

# ╔═╡ a0380595-bdaa-4e88-b2af-e2800ab74fb6
begin
	g_no_cap = evaluate(opf_no_cap.g[1]);
	plt1 = bar(g_no_cap, size=(800, 200))
	plt2 = bar(g_static - g_no_cap, size=(800, 200))
	plot(plt1, plt2)
end

# ╔═╡ bfbbdaa0-b75a-416c-99be-35daf3492cb7
@show norm(g_no_cap - g_static) / sqrt(l)

# ╔═╡ 1e29eeb1-ff5d-41c2-8e10-081636569b78
md"
## Solve dynamic problem, with storage
"

# ╔═╡ e79d155e-df93-4854-9701-bffbf4fbe2cd
C = rand(Uniform(0.25, 0.75), n) .* (B*gmax);

# ╔═╡ fb2c65d0-62aa-4b14-92a7-951c1959030d
P = rand(Uniform(0.2, 0.3), n) #.* C

# ╔═╡ bd45770e-6f4e-4420-9818-628c1e521dd4
net_dynamic = DynamicPowerNetwork(fqs, fls, pmaxs, gmaxs, A, B, P, C);

# ╔═╡ bc9e7ae7-41b7-4eba-8e40-eedeb0ea47bd
begin
	opf_dynamic = DynamicPowerManagementProblem(net_dynamic, ds)
	solve!(opf_dynamic, ECOS.Optimizer, verbose=true)
end;

# ╔═╡ f0d4f300-af94-4ed9-8e1d-e27d76f0de69
begin
	Random.seed!(949)
	
	plot(xlabel="T", ylabel="s", size=(750, 250), legend=:topleft)
	for i in [1, 2, 6]
		si = [evaluate(opf_dynamic.s[t][i])[1] for t in 1:T]
		plot!(0:T, [0; si], lw=4, label=i)
	end
	plot!()
end

# ╔═╡ 14cd88e7-fb7d-4844-b12c-1c452984acb5
initial_storage = evaluate(opf_dynamic.s[1]);

# ╔═╡ 8fdf7946-a939-496c-be3c-5d49fe755675
@show mean(initial_storage)

# ╔═╡ 182a657c-ed7f-47d9-8c4e-9f0715e39c84
@show maximum(initial_storage)

# ╔═╡ 4f8ea4df-8d96-4053-965b-73fd4bb0787d
md"
## Compute marginal emission factors
"

# ╔═╡ 85c422fd-f3a7-4e8f-81c5-47e8934e7558
c = rand(Exponential(2), l) .+ 1

# ╔═╡ 7214657b-270c-4387-89a9-3d4994f9212e
# Static MEFs
mef_static = compute_mefs(opf_static, net, ds[1], c);

# ╔═╡ 09a5a187-7ec4-4e91-8dc8-619bc9f7653e
t = 1

# ╔═╡ a2a0f75c-9adc-4a63-9291-205d92036df2
# Dynamic MEFs, without storage	
mef_no_cap = compute_mefs(opf_no_cap, net_no_cap, ds, c, t);

# ╔═╡ f40d58fd-2790-4733-a660-81ec67603b5e
@show norm(mef_static - mef_no_cap)

# ╔═╡ 7edf8ad7-5371-4684-85fe-73df835492b2
# Dynamic MEFs, with storage
mef_dynamic = compute_mefs(opf_dynamic, net_dynamic, ds, c, t);

# ╔═╡ 3739a541-1687-49b8-9c52-e08da51a3014
begin
	plot(ylabel="MEF", xlabel="node", size=(600, 300), legend=:topleft)
	bar!(1:n, mef_no_cap, bar_width=0.2, label="no storage")
	bar!((1:n) .+ 0.25, mef_dynamic, bar_width=0.2, label="storage")
end

# ╔═╡ 698a2598-6a72-4c7d-b556-75e1313baf77
@show norm(mef_static - mef_dynamic)

# ╔═╡ fbc4318c-216d-4408-a3d2-95c25c8d9b99
md"
## How do MEFs vary with storage penetration?
"

# ╔═╡ d941788d-61f0-4927-9a29-5b0c519eed08
storage_penetrations = [0.0, 0.10, 0.25, 0.50, 1.0]

# ╔═╡ 642c826c-6b2b-4626-95f4-cef256046784
begin
	results = []
	for C_rel in storage_penetrations
		# Construct problem
		C = C_rel * (B*gmax)
		P = C / 2
		net_C = DynamicPowerNetwork(fqs, fls, pmaxs, gmaxs, A, B, P, C)
		opf_C = DynamicPowerManagementProblem(net_C, ds)
		
		# Solve
		solve!(opf_C, OPT, verbose=true)
		println("Starting C=$C_rel")
		
		# Compute MEFs
		mefs = [compute_mefs(opf_C, net_C, ds, c, t) for t in 1:T]
		push!(results, mefs)
		println("C=$C_rel done")
	end
end

# ╔═╡ 74419b49-55ce-40d4-bce7-d6a1db40448a
@bind node Slider(1:n)

# ╔═╡ b3b68d17-a83b-4aa7-93e8-b740855e58fd
begin
	plot(size=(600, 250), title="node=$node", legend=:bottomleft)
	plot!(xlabel="t", ylabel="CO2 / MW")
	
	for (idx, C_rel) in enumerate(storage_penetrations)
		# Construct MEF over time
		mef_temporal = [results[idx][t][node] for t in 1:T]
		
		# Plot
		plot!(mef_temporal, lw=4, label="C=$C_rel")
	end
	plot!() #ylim=(0, 6))
end

# ╔═╡ Cell order:
# ╠═26a4178a-de77-11eb-1fa3-5fc861d32904
# ╟─c8406984-bd2c-40fd-92b9-397e07d9f213
# ╟─ad0d7c81-c88e-4d21-84e0-1b8e8f44bf5c
# ╠═b3ce8622-b8dd-4f9b-8593-4bd21468285d
# ╟─99af5f08-2ad3-4ee0-840f-474e20a05ea3
# ╠═edd937f1-c5f3-45dc-94e9-84234858db2f
# ╠═bcb8825f-afb0-4f5b-a25d-4891579f9d39
# ╟─eea80101-05d7-4f01-bab8-b89d0aef7911
# ╠═ea6cb153-e45c-480b-a36d-4cf7711def62
# ╠═f4ab0602-3954-4a5d-b770-3d9e455b500a
# ╠═c6e83723-ab3d-41e0-b323-07283ce72281
# ╟─ae595fbf-f501-49e1-93f4-217646c29396
# ╠═88490287-4f63-4c45-bf58-18df9fe711b6
# ╠═af67cc5a-ca58-448b-912a-f438172deba8
# ╠═7321565d-9eab-4871-8636-4bc6a57bc10f
# ╠═dbb14945-5fea-479d-be31-e9c14e1b8672
# ╠═a0380595-bdaa-4e88-b2af-e2800ab74fb6
# ╠═bfbbdaa0-b75a-416c-99be-35daf3492cb7
# ╟─1e29eeb1-ff5d-41c2-8e10-081636569b78
# ╠═e79d155e-df93-4854-9701-bffbf4fbe2cd
# ╠═fb2c65d0-62aa-4b14-92a7-951c1959030d
# ╠═bd45770e-6f4e-4420-9818-628c1e521dd4
# ╠═bc9e7ae7-41b7-4eba-8e40-eedeb0ea47bd
# ╠═f0d4f300-af94-4ed9-8e1d-e27d76f0de69
# ╠═14cd88e7-fb7d-4844-b12c-1c452984acb5
# ╠═8fdf7946-a939-496c-be3c-5d49fe755675
# ╠═182a657c-ed7f-47d9-8c4e-9f0715e39c84
# ╟─4f8ea4df-8d96-4053-965b-73fd4bb0787d
# ╠═85c422fd-f3a7-4e8f-81c5-47e8934e7558
# ╠═7214657b-270c-4387-89a9-3d4994f9212e
# ╠═09a5a187-7ec4-4e91-8dc8-619bc9f7653e
# ╠═a2a0f75c-9adc-4a63-9291-205d92036df2
# ╠═f40d58fd-2790-4733-a660-81ec67603b5e
# ╠═7edf8ad7-5371-4684-85fe-73df835492b2
# ╠═3739a541-1687-49b8-9c52-e08da51a3014
# ╠═698a2598-6a72-4c7d-b556-75e1313baf77
# ╟─fbc4318c-216d-4408-a3d2-95c25c8d9b99
# ╠═d941788d-61f0-4927-9a29-5b0c519eed08
# ╠═642c826c-6b2b-4626-95f4-cef256046784
# ╠═74419b49-55ce-40d4-bce7-d6a1db40448a
# ╠═b3b68d17-a83b-4aa7-93e8-b740855e58fd
