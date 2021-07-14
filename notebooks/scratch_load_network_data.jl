### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 7b1febd2-e4ca-11eb-1907-d15ce3b2a418
begin
	using Pkg; Pkg.activate()
	using PowerModels
	using SparseArrays
	using Plots
end

# ╔═╡ 83a3ff1b-cb31-4494-ac42-771319888e88
begin
	network_data = PowerModels.parse_file("../data/case118.m");
	net = make_basic_network(network_data)
	make_per_unit!(net)
end

# ╔═╡ 17b115c1-fc14-4dfc-a9ff-32f1119d42c9
net["branch"]["1"]["rate_a"]

# ╔═╡ 72d8fad5-040d-4e2b-befc-14a7db1932b1
md"""
## Functions
"""

# ╔═╡ 70713496-b390-4b05-9d7e-22e974971dc1
"""
	make_B(gen, n, l)
"""
function make_B(gen, n, l)
	nodes = [gen[i]["gen_bus"] for i in string.(1:l)]
	
	B = spzeros(n, l)
	for (ind_g, ind_n) in enumerate(nodes)
		B[ind_n, ind_g] = 1.0
	end
	
	return B
end

# ╔═╡ 706642f6-3f5a-416f-993a-915659acc868
"""
	make_d(load, n)
"""
function make_d(load, n)
	n_load = length(load)
	nodes = [load[i]["load_bus"] for i in string.(1:n_load)]
	#dp = [load[i]["pd"] for i in string.(1:n_load)]
	
	d = zeros(n)
	for (ind_d, ind_n) in enumerate(nodes)
		d[ind_n] = load[string(ind_d)]["pd"]
	end
	
	return d
end

# ╔═╡ 244d1250-7d4d-4cb0-9a89-af29a926b768
begin
	base_mva = net["baseMVA"]
	gen, load, bus, branch = net["gen"], net["load"], net["bus"], net["branch"]
	
	# Dimensions
	l = length(net["gen"])
	n = length(net["bus"])
	m = length(net["branch"])
	
	# Topology
	A = calc_basic_incidence_matrix(net)'
	B = make_B(gen, n, l)
	
	@assert size(A) == (n, m)
	@assert size(B) == (n, l)
	
	# Capacities
	pmax = [branch[i]["rate_a"] for i in string.(1:m)]
	gmax = [gen[i]["mbase"] * gen[i]["pmax"] for i in string.(1:l)] / base_mva
	
	# Generator costs
	fq = [gen[i]["cost"][1] for i in string.(1:l)]
	fl = [gen[i]["cost"][2] for i in string.(1:l)]
	
	# Carbon cost
	c = nothing
	
	# Demand
	d = make_d(load, n)
	
	n, m, l
end

# ╔═╡ 99a14cd9-49ae-4f28-a0ad-1b212728176b
begin
	plot()
	histogram(pmax[pmax .< 10], label="pmax")
	histogram!(gmax, label="gmax", alpha=0.75)
	histogram!(d[d .!= 0], label="d", alpha=0.5)
end

# ╔═╡ Cell order:
# ╠═7b1febd2-e4ca-11eb-1907-d15ce3b2a418
# ╠═83a3ff1b-cb31-4494-ac42-771319888e88
# ╠═17b115c1-fc14-4dfc-a9ff-32f1119d42c9
# ╠═244d1250-7d4d-4cb0-9a89-af29a926b768
# ╠═99a14cd9-49ae-4f28-a0ad-1b212728176b
# ╟─72d8fad5-040d-4e2b-befc-14a7db1932b1
# ╠═70713496-b390-4b05-9d7e-22e974971dc1
# ╠═706642f6-3f5a-416f-993a-915659acc868
