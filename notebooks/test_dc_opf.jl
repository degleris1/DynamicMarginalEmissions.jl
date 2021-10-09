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

# ╔═╡ 584da6cc-383a-4fa2-82f9-4cfd72c8a2e8
using Distributions

# ╔═╡ d3bd8561-248b-4e19-b556-c4274cd0e77c
using ForwardDiff


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
	net.pmax ./= 25
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

# ╔═╡ 961403c7-1452-4a2c-86a4-3cc9a442df50
rms(x) = norm(x) / length(x)

# ╔═╡ 67f95780-d48f-4a8b-b830-2fdc10f5ec25
begin
	c = zeros(kkt_dims(n, m, l))
	c[1:l] .= rand(length(net.fq))
	
	c
end

# ╔═╡ 4b52f0c1-f51e-4645-b36f-f5893eb6c275
x = flatten_variables(opf)

# ╔═╡ 9f207ea7-f01c-4519-b586-8f30e70a1aa7
K = kkt(x, net, d)

# ╔═╡ 74af0dc5-9ce5-40f7-b11e-f2a45e814773
rms(K)

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

# ╔═╡ c5dc8c59-a7c4-4bad-95f9-fdc187f28c27
T = 3

# ╔═╡ dd7eae71-4b65-47b5-8d7f-dfb94768538a
d_dyn = [rand(Uniform(0.95, 1.05), n) .* d for t in 1:T]

# ╔═╡ a6b0715e-d836-4f16-9c3f-6725ad176d6f
η_c, η_d = 0.95, 0.95

# ╔═╡ 0d66e4f3-0f83-46ad-b7cd-c34ea6860bce
C = 0.02 * sum(d_dyn)

# ╔═╡ 2283c971-9d13-4ef6-a7e3-b038742284a7
P = 0.25 * C

# ╔═╡ 0c47e14b-24db-44c8-93cd-a2585d46d732
dnet = make_dynamic(net, T, P, C, [net.gmax for _ in 1:T], η_c, η_d);

# ╔═╡ 6282d5d8-24e1-4f2f-98de-4bb7829b3786
begin
	dopf = DynamicPowerManagementProblem(dnet, d_dyn)
	solve_ecos!(dopf)
	
	dopf.problem.status, dopf.problem.optval
end

# ╔═╡ 6d319e6d-6214-4039-b8ca-085ab42ac445
xdyn = flatten_variables_dyn(dopf)[:, 1]

# ╔═╡ 88e3b8ed-c459-40fd-bec2-1e8ac966a364
Kdyn = kkt_dyn(xdyn, dnet, d_dyn)

# ╔═╡ 79a1456e-a778-4dd1-9991-8087500a1d42
rms(Kdyn)

# ╔═╡ 42fe6178-f3ff-4e42-8d61-4a98d96d3101
md"""
## Dynamic jacobian
"""

# ╔═╡ 4e1c81b0-52bb-40a4-a708-ecd18ad3e148
DK_xdyn = compute_jacobian_kkt_dyn(xdyn, dnet, d_dyn)

# ╔═╡ a518b330-cd00-4253-8f5b-8c7443f5874d
DK_xdyn_zyg = ForwardDiff.jacobian(x -> kkt_dyn(x, dnet, d_dyn), xdyn)

# ╔═╡ d3c64336-2df7-4935-9d58-a1fb4bdb0a4e
err = (DK_xdyn_zyg - DK_xdyn)[length(x)*3:end, :];

# ╔═╡ 628cfcc8-fb29-4fd0-a6a2-e76fe1c89c5b
err

# ╔═╡ fe39cfab-9373-473b-a939-e1a64bdd0262
sum(x -> x != 0, err)

# ╔═╡ e8024227-3633-4671-baf9-bc5a19b1794f
md"""
## Could this be the Zygote bug again?
"""

# ╔═╡ 70c92bb8-c95f-40ba-8d2d-89f36a0f1ee4
jac_man(x) = compute_jacobian_kkt_dyn(x, dnet, d_dyn)

# ╔═╡ 7faee8eb-d78d-4576-af8b-e2ea4d0e5ada
jac_zyg(x) = ForwardDiff.jacobian(x -> kkt_dyn(x, dnet, d_dyn), x)

# ╔═╡ c5aaa8de-a2ac-41c8-9474-761d57225f2b
norm(jac_man(xdyn) - jac_man(xdyn .+ 1e-5))

# ╔═╡ 4d609864-4fab-409d-a915-c08bdd8301dc
norm(jac_zyg(xdyn) - jac_zyg(xdyn .+ 1e-5))

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
# ╠═961403c7-1452-4a2c-86a4-3cc9a442df50
# ╠═67f95780-d48f-4a8b-b830-2fdc10f5ec25
# ╠═4b52f0c1-f51e-4645-b36f-f5893eb6c275
# ╠═9f207ea7-f01c-4519-b586-8f30e70a1aa7
# ╠═74af0dc5-9ce5-40f7-b11e-f2a45e814773
# ╠═e7b8e288-7ac4-4b29-9ec6-513570ecfac9
# ╠═aab8766e-2dee-4cf5-9857-e7b0b58239b8
# ╠═55efbe79-01ce-407b-8f92-2ac81b9d811c
# ╠═6e0c7b10-b8ef-4351-bedf-ce95f5c24405
# ╠═ab8949f7-8271-4c8e-8ce3-54ea8feb834a
# ╟─fcb08184-9f13-44df-ada0-bcb5674011f9
# ╠═584da6cc-383a-4fa2-82f9-4cfd72c8a2e8
# ╠═c5dc8c59-a7c4-4bad-95f9-fdc187f28c27
# ╠═dd7eae71-4b65-47b5-8d7f-dfb94768538a
# ╠═a6b0715e-d836-4f16-9c3f-6725ad176d6f
# ╠═0d66e4f3-0f83-46ad-b7cd-c34ea6860bce
# ╠═2283c971-9d13-4ef6-a7e3-b038742284a7
# ╠═0c47e14b-24db-44c8-93cd-a2585d46d732
# ╠═6282d5d8-24e1-4f2f-98de-4bb7829b3786
# ╠═6d319e6d-6214-4039-b8ca-085ab42ac445
# ╠═88e3b8ed-c459-40fd-bec2-1e8ac966a364
# ╠═79a1456e-a778-4dd1-9991-8087500a1d42
# ╠═42fe6178-f3ff-4e42-8d61-4a98d96d3101
# ╠═d3bd8561-248b-4e19-b556-c4274cd0e77c
# ╠═4e1c81b0-52bb-40a4-a708-ecd18ad3e148
# ╠═a518b330-cd00-4253-8f5b-8c7443f5874d
# ╠═d3c64336-2df7-4935-9d58-a1fb4bdb0a4e
# ╠═628cfcc8-fb29-4fd0-a6a2-e76fe1c89c5b
# ╠═fe39cfab-9373-473b-a939-e1a64bdd0262
# ╟─e8024227-3633-4671-baf9-bc5a19b1794f
# ╠═70c92bb8-c95f-40ba-8d2d-89f36a0f1ee4
# ╠═7faee8eb-d78d-4576-af8b-e2ea4d0e5ada
# ╠═c5aaa8de-a2ac-41c8-9474-761d57225f2b
# ╠═4d609864-4fab-409d-a915-c08bdd8301dc
