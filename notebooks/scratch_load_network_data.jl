### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 7b1febd2-e4ca-11eb-1907-d15ce3b2a418
begin
	using Pkg; Pkg.activate()
	using Plots
	
	using Revise
	using CarbonNetworks
	
	theme(:default, label=nothing)
end

# ╔═╡ c21d3111-159e-47b9-8c8d-4dee105fd527
md"""
## Get network data
"""

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

# ╔═╡ c68d38d4-7057-4e69-8e1a-921061a555f8
md"""
## Get load and solar data
"""

# ╔═╡ e57a9392-49d3-483e-8de5-f24a71bbfd00
load_data = load_demand_data("2021_07_01", normalize_rows=true);

# ╔═╡ 947fc22a-f85b-46c5-bc9e-6c7510292040
plot(load_data, lw=3, alpha=0.7, palette=palette([:SkyBlue, :RoyalBlue], 24))

# ╔═╡ eda9ae84-5950-4602-800a-1e271d87195c
renew_data, renew_labels = load_renewable_data("2021_07_01"; normalize_rows=true);

# ╔═╡ 2204b359-8913-4ac5-8f80-a772045d0cb6
begin
	T = size(renew_data, 2)
	
	# Color according to resource type
	fc = l -> l == "solar" ? fill(:Orange, T) : fill(:Blue, T)
	c = hcat([fc(l) for l in renew_labels]...)
	
	plot(renew_data, lw=4, alpha=0.75, c=c)
end

# ╔═╡ Cell order:
# ╠═7b1febd2-e4ca-11eb-1907-d15ce3b2a418
# ╟─c21d3111-159e-47b9-8c8d-4dee105fd527
# ╠═83a3ff1b-cb31-4494-ac42-771319888e88
# ╠═99a14cd9-49ae-4f28-a0ad-1b212728176b
# ╠═49d0fa29-4a52-4b07-a110-95876a353041
# ╟─c68d38d4-7057-4e69-8e1a-921061a555f8
# ╠═e57a9392-49d3-483e-8de5-f24a71bbfd00
# ╠═947fc22a-f85b-46c5-bc9e-6c7510292040
# ╠═eda9ae84-5950-4602-800a-1e271d87195c
# ╠═2204b359-8913-4ac5-8f80-a772045d0cb6
