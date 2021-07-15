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
	
	using Revise
	using CarbonNetworks
end

# ╔═╡ 83a3ff1b-cb31-4494-ac42-771319888e88
net, d, _ = load_synthetic_network("case14.m");

# ╔═╡ 99a14cd9-49ae-4f28-a0ad-1b212728176b
begin
	plot()
	histogram(net.pmax, label="pmax")
	histogram!(net.gmax, label="gmax", alpha=0.75)
	histogram!(d[d .!= 0], label="d", alpha=0.5)
end

# ╔═╡ 49d0fa29-4a52-4b07-a110-95876a353041
net.A

# ╔═╡ Cell order:
# ╠═7b1febd2-e4ca-11eb-1907-d15ce3b2a418
# ╠═83a3ff1b-cb31-4494-ac42-771319888e88
# ╠═99a14cd9-49ae-4f28-a0ad-1b212728176b
# ╠═49d0fa29-4a52-4b07-a110-95876a353041
