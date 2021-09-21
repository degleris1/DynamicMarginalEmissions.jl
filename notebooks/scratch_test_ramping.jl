### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ ef1bccb8-1b11-11ec-02b5-af1e0e589de0
begin
	using Pkg; Pkg.activate("")
	using Convex, ECOS
end

# ╔═╡ 589ab2e5-319a-487c-855d-508af37a277f
using Revise; using CarbonNetworks

# ╔═╡ 78dff8f2-9454-414d-8767-02b928d6daa6
using SparseArrays

# ╔═╡ 575745b0-7b13-47b7-8021-63fa29c80b34
using Zygote

# ╔═╡ cab7047e-8b71-4652-a668-4016349b0f3c
solve_ecos!(p) = solve!(p, () -> ECOS.Optimizer(verbose=true))

# ╔═╡ d0f3079b-d3db-48a4-923f-c977c7ce0c34
md"""
## Make random data
"""

# ╔═╡ 0b98c228-3413-49ef-805e-0e0bd46f7373
n, l, T = 5, 20, 3

# ╔═╡ df56e6de-d953-47c4-b846-5338e98cf597
A, B, fq, fl, d, gmax, pmax, P, C = generate_random_data(n, l, T);

# ╔═╡ c62d92b0-b850-47fe-b521-0ad63609e5f0
m = size(A, 2)

# ╔═╡ 26345b70-ae32-423f-9aa3-a61a78857ca0
ρ = 0.05 * gmax[1]

# ╔═╡ c9fd438b-d43a-4a84-be49-12d048ee7e72
size(ρ)

# ╔═╡ 681e1106-e41b-4417-a958-c5a23b734fe3
dnet = DynamicPowerNetwork(fq, fl, pmax, gmax, A, B, P, C, T; 
	η_c=0.95, η_d=0.95, ρ=ρ);

# ╔═╡ 837f4451-f246-436a-bf30-5d5d70791d34
begin
	dmin = DynamicPowerManagementProblem(dnet, d)
	solve_ecos!(dmin)
end;

# ╔═╡ 9edb5562-3fc0-4abe-aea6-d2dd859ecf5c
dmin.problem.status

# ╔═╡ 024d1409-ba5f-401e-b28d-ace0b8c024b4
md"""
## Check KKT conditions
"""

# ╔═╡ fb761e75-9607-4675-b91c-807e4ba4495e
x = flatten_variables_dyn(dmin);

# ╔═╡ 838ade3d-8a86-4683-8f13-49a48292338d
K = kkt_dyn(x, dnet, d);

# ╔═╡ d5876c0e-024d-4fe0-a7d7-c12a1b9d96e8
length(K)

# ╔═╡ 97d40580-6dbe-46dc-912e-e9412d94227d
sum(abs, K)

# ╔═╡ f10c1627-3578-41e0-993f-3f3c313c2372
md"""
## Check Jacobian
"""

# ╔═╡ ddbc75b5-9dc4-4fb4-81e7-be2aafb556cd
y, ∂KT_zyg = Zygote.forward_jacobian(x -> kkt_dyn(x, dnet, d), x);

# ╔═╡ d5021fd3-064d-4d56-ae7c-a2a13b468f05
∂K_zyg = sparse(adjoint(∂KT_zyg))

# ╔═╡ Cell order:
# ╠═ef1bccb8-1b11-11ec-02b5-af1e0e589de0
# ╠═589ab2e5-319a-487c-855d-508af37a277f
# ╠═cab7047e-8b71-4652-a668-4016349b0f3c
# ╠═d0f3079b-d3db-48a4-923f-c977c7ce0c34
# ╠═0b98c228-3413-49ef-805e-0e0bd46f7373
# ╠═df56e6de-d953-47c4-b846-5338e98cf597
# ╠═c62d92b0-b850-47fe-b521-0ad63609e5f0
# ╠═26345b70-ae32-423f-9aa3-a61a78857ca0
# ╠═c9fd438b-d43a-4a84-be49-12d048ee7e72
# ╠═681e1106-e41b-4417-a958-c5a23b734fe3
# ╠═837f4451-f246-436a-bf30-5d5d70791d34
# ╠═9edb5562-3fc0-4abe-aea6-d2dd859ecf5c
# ╠═024d1409-ba5f-401e-b28d-ace0b8c024b4
# ╠═fb761e75-9607-4675-b91c-807e4ba4495e
# ╠═838ade3d-8a86-4683-8f13-49a48292338d
# ╠═d5876c0e-024d-4fe0-a7d7-c12a1b9d96e8
# ╠═97d40580-6dbe-46dc-912e-e9412d94227d
# ╠═f10c1627-3578-41e0-993f-3f3c313c2372
# ╠═78dff8f2-9454-414d-8767-02b928d6daa6
# ╠═575745b0-7b13-47b7-8021-63fa29c80b34
# ╠═ddbc75b5-9dc4-4fb4-81e7-be2aafb556cd
# ╠═d5021fd3-064d-4d56-ae7c-a2a13b468f05
