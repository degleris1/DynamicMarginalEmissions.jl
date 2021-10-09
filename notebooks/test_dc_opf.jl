### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 74d94937-1cc3-49d3-b56a-d857a4a9c595
using Pkg; Pkg.activate("")

# ╔═╡ 093f668a-27b1-11ec-0f49-7d244fb901b7
begin
	using Revise
	using CarbonNetworks
end

# ╔═╡ f49550c4-ad45-41ee-a78f-231158a3162a
using Convex, ECOS

# ╔═╡ 1316c162-4145-473a-b947-8f24329dea8e
using Plots

# ╔═╡ f246ca6e-9aa9-4552-90f7-e299c7bfa02e
using LinearAlgebra

# ╔═╡ 399aeec4-f6b2-4225-bf7b-1cdbf5a8afc6
solve_ecos!(x) = solve!(x, () -> ECOS.Optimizer())

# ╔═╡ 89ed6548-c338-4056-abd3-ad75689f8015
md"""
## Load data, solve problem
"""

# ╔═╡ d4717e90-fde2-4d2a-9b51-ab4cdd45f06a
begin
	net, d, casedata = load_synthetic_network("case14.m")
	
	net.fq .+= 10 * rand(length(net.fq))
	net.fl .+= 100 * rand(length(net.fl))
	net.pmax ./= 27
end;

# ╔═╡ 68c6ca27-c096-4252-be3f-f4348c01fcf4
begin
	opf = PowerManagementProblem(net, d)
	solve_ecos!(opf)
	
	opf.problem.status, opf.problem.optval
end

# ╔═╡ df3e9efd-abb8-4c15-9972-24b3b2f21d45
bar(abs.(evaluate(opf.p)) ./ net.pmax, size=(600, 200), 
	label=nothing, ylabel="p / pmax")

# ╔═╡ 1674446e-503f-4a6b-b715-15f3025893bc
heatmap(opf.params.F, 
	yflip=true, c=cgrad(:bluesreds, [0.0]), clim=(-1, 1), size=(300, 300)
)

# ╔═╡ 46e74004-cbbd-4657-b176-3c352013ccf6
n, m = size(net.A)

# ╔═╡ 3b6380d3-3d13-456c-b7ce-ee183dd5e836
_, l = size(net.B)

# ╔═╡ 09da247a-c4ff-492b-bf9f-4ca4a4ea7ce0
md"""
## Analyze sensitivities
"""

# ╔═╡ 67f95780-d48f-4a8b-b830-2fdc10f5ec25
begin
	c = zeros(kkt_dims(n, m, l))
	c[1:l] .= rand(length(net.fq))
	
	c
end

# ╔═╡ 4b52f0c1-f51e-4645-b36f-f5893eb6c275
x = flatten_variables(opf)

# ╔═╡ e7b8e288-7ac4-4b29-9ec6-513570ecfac9
DK_x = CarbonNetworks.compute_jacobian_kkt(
	net.fq, net.fl, d, net.pmax, net.gmax, net.A, net.B, net.F, 
	x
)

# ╔═╡ aab8766e-2dee-4cf5-9857-e7b0b58239b8
mefs = sensitivity_demand(opf, c, net, d)

# ╔═╡ 55efbe79-01ce-407b-8f92-2ac81b9d811c
lmps = get_lmps(opf)

# ╔═╡ 6e0c7b10-b8ef-4351-bedf-ce95f5c24405
bar(mefs, size=(600, 200), label=nothing, ylabel="mef")

# ╔═╡ ab8949f7-8271-4c8e-8ce3-54ea8feb834a
bar(lmps, size=(600, 200), label=nothing, ylabel="lmp")

# ╔═╡ fcb08184-9f13-44df-ada0-bcb5674011f9
md"""
## Dynamic model

Does it still work?
"""

# ╔═╡ 0c47e14b-24db-44c8-93cd-a2585d46d732


# ╔═╡ Cell order:
# ╠═74d94937-1cc3-49d3-b56a-d857a4a9c595
# ╠═093f668a-27b1-11ec-0f49-7d244fb901b7
# ╠═f49550c4-ad45-41ee-a78f-231158a3162a
# ╠═399aeec4-f6b2-4225-bf7b-1cdbf5a8afc6
# ╠═1316c162-4145-473a-b947-8f24329dea8e
# ╟─89ed6548-c338-4056-abd3-ad75689f8015
# ╠═d4717e90-fde2-4d2a-9b51-ab4cdd45f06a
# ╠═68c6ca27-c096-4252-be3f-f4348c01fcf4
# ╠═df3e9efd-abb8-4c15-9972-24b3b2f21d45
# ╠═1674446e-503f-4a6b-b715-15f3025893bc
# ╠═46e74004-cbbd-4657-b176-3c352013ccf6
# ╠═3b6380d3-3d13-456c-b7ce-ee183dd5e836
# ╟─09da247a-c4ff-492b-bf9f-4ca4a4ea7ce0
# ╠═f246ca6e-9aa9-4552-90f7-e299c7bfa02e
# ╠═67f95780-d48f-4a8b-b830-2fdc10f5ec25
# ╠═4b52f0c1-f51e-4645-b36f-f5893eb6c275
# ╠═e7b8e288-7ac4-4b29-9ec6-513570ecfac9
# ╠═aab8766e-2dee-4cf5-9857-e7b0b58239b8
# ╠═55efbe79-01ce-407b-8f92-2ac81b9d811c
# ╠═6e0c7b10-b8ef-4351-bedf-ce95f5c24405
# ╠═ab8949f7-8271-4c8e-8ce3-54ea8feb834a
# ╟─fcb08184-9f13-44df-ada0-bcb5674011f9
# ╠═0c47e14b-24db-44c8-93cd-a2585d46d732
