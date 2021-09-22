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

# ╔═╡ f5afa002-63d8-4b8a-81a2-a867aeb8d612
using Plots

# ╔═╡ 78dff8f2-9454-414d-8767-02b928d6daa6
using SparseArrays

# ╔═╡ 575745b0-7b13-47b7-8021-63fa29c80b34
using Zygote

# ╔═╡ f84b6daf-ade9-4ef0-b0d9-88a93d27ea1a
using ForwardDiff

# ╔═╡ cab7047e-8b71-4652-a668-4016349b0f3c
solve_ecos!(p) = solve!(p, () -> ECOS.Optimizer(verbose=true))

# ╔═╡ d0f3079b-d3db-48a4-923f-c977c7ce0c34
md"""
## Make random data
"""

# ╔═╡ 0b98c228-3413-49ef-805e-0e0bd46f7373
n, l, T = 5, 40, 5

# ╔═╡ df56e6de-d953-47c4-b846-5338e98cf597
A, B, fq, fl, d, gmax, pmax, P, C = generate_random_data(n, l, T);

# ╔═╡ c62d92b0-b850-47fe-b521-0ad63609e5f0
m = size(A, 2)

# ╔═╡ 024d1409-ba5f-401e-b28d-ace0b8c024b4
md"""
## Check KKT conditions
"""

# ╔═╡ f10c1627-3578-41e0-993f-3f3c313c2372
md"""
## Check Jacobian
"""

# ╔═╡ e1519269-475b-409f-9f34-a45769877f5b
md"""
## Look at error
"""

# ╔═╡ 488cbe38-c106-4d19-ad34-363d763a4a15
α = 0.2 #+ eps()

# Weird values...

# n =5 and 
# α = 0.2

# n=20 and
# α = 0.79999999 OR

# Others...

# ╔═╡ 26345b70-ae32-423f-9aa3-a61a78857ca0
ρ = α * gmax[1]

# ╔═╡ 681e1106-e41b-4417-a958-c5a23b734fe3
dnet = DynamicPowerNetwork(fq, fl, pmax, gmax, A, B, P, C, T; 
	η_c=0.95, η_d=0.95, ρ=ρ);

# ╔═╡ 837f4451-f246-436a-bf30-5d5d70791d34
begin
	dmin = DynamicPowerManagementProblem(dnet, d)
	solve_ecos!(dmin)
end;

# ╔═╡ fb761e75-9607-4675-b91c-807e4ba4495e
x = flatten_variables_dyn(dmin);

# ╔═╡ 838ade3d-8a86-4683-8f13-49a48292338d
K = kkt_dyn(x, dnet, d);

# ╔═╡ d5876c0e-024d-4fe0-a7d7-c12a1b9d96e8
length(K)

# ╔═╡ 97d40580-6dbe-46dc-912e-e9412d94227d
sum(abs, K)

# ╔═╡ d55c0cbe-8463-4d97-b75c-614a902728ea
∂K_fwd = ForwardDiff.jacobian(x -> kkt_dyn(x, dnet, d), x);

# ╔═╡ ddbc75b5-9dc4-4fb4-81e7-be2aafb556cd
y, ∂KT_zyg = Zygote.forward_jacobian(x -> kkt_dyn(x, dnet, d), x);

# ╔═╡ d5021fd3-064d-4d56-ae7c-a2a13b468f05
∂K_zyg = sparse(adjoint(∂KT_zyg));

# ╔═╡ b4b13c32-9572-4411-90e7-e7ff4667737b
sum(abs, ∂K_fwd - ∂K_zyg)

# ╔═╡ 3e9e7c0c-643a-4ad9-bbb8-d02a09d0fde8
∂K = compute_jacobian_kkt_dyn(x, dnet, d);

# ╔═╡ 9edb5562-3fc0-4abe-aea6-d2dd859ecf5c
dmin.problem.status

# ╔═╡ f5eb20a0-6468-411f-8569-42969c8cf2be
err = ∂K - ∂K_zyg;

# ╔═╡ 5f91864e-8469-4b78-92c9-80337c23eea5
sum(abs, ∂K)

# ╔═╡ c523ce28-f026-4b67-bde2-d0dac3e9eeec
sum(abs, ∂K_zyg)

# ╔═╡ 79c24ade-65b1-4ac9-8162-4bed95bb85bf
sum(abs, err)

# ╔═╡ be6e4906-e560-4258-98f9-a273fea2bba9
sum(!=(0), err)

# ╔═╡ 835d31ae-03a7-4839-84e2-d5e21e3a15cd
err.nzval

# ╔═╡ 5873ffa5-5962-471f-a579-32470d555f7a
findall(!=(0), err)

# ╔═╡ 1f8a655e-a914-468b-98df-f9a3cfd19171
∂K[705, :].nzval

# ╔═╡ 668a9da8-41b7-4931-9637-c38e299265cd
∂K_zyg[705, :].nzval

# ╔═╡ 0661ea23-a9ce-4530-aa19-9040d4296c25
kkt_dims(n,m,l)*T

# ╔═╡ 3561ab12-574f-4a57-b60c-432532df49e5
kkt_dims(n,m,l)*T + 3n

# ╔═╡ 5d0da605-8cc5-46be-91d3-e004c31d05bc
plot(
	heatmap(Matrix(∂K[701:705, 716:720]), yflip=true),
	heatmap(Matrix(∂K_zyg[701:705, 716:720]), yflip=true),
)

# ╔═╡ c449e275-3353-4296-8a2f-7f8bfa6132c9
b = rand(size(∂K, 2));

# ╔═╡ df67d1b6-71bf-4792-80e3-6eac11312e24
v1 = ∂K \ b;

# ╔═╡ 2bf3d695-ffed-4745-9277-d4668b3d7ae0
sum(abs, ∂K*v1 - b)

# ╔═╡ 5edb6c01-66fe-4c18-a8ee-6078e8255fe8
v2 = ∂K_zyg \ b;

# ╔═╡ Cell order:
# ╠═ef1bccb8-1b11-11ec-02b5-af1e0e589de0
# ╠═589ab2e5-319a-487c-855d-508af37a277f
# ╠═cab7047e-8b71-4652-a668-4016349b0f3c
# ╠═f5afa002-63d8-4b8a-81a2-a867aeb8d612
# ╠═d0f3079b-d3db-48a4-923f-c977c7ce0c34
# ╠═0b98c228-3413-49ef-805e-0e0bd46f7373
# ╠═df56e6de-d953-47c4-b846-5338e98cf597
# ╠═c62d92b0-b850-47fe-b521-0ad63609e5f0
# ╠═26345b70-ae32-423f-9aa3-a61a78857ca0
# ╠═681e1106-e41b-4417-a958-c5a23b734fe3
# ╠═837f4451-f246-436a-bf30-5d5d70791d34
# ╠═024d1409-ba5f-401e-b28d-ace0b8c024b4
# ╠═fb761e75-9607-4675-b91c-807e4ba4495e
# ╠═838ade3d-8a86-4683-8f13-49a48292338d
# ╠═d5876c0e-024d-4fe0-a7d7-c12a1b9d96e8
# ╠═97d40580-6dbe-46dc-912e-e9412d94227d
# ╠═f10c1627-3578-41e0-993f-3f3c313c2372
# ╠═78dff8f2-9454-414d-8767-02b928d6daa6
# ╠═575745b0-7b13-47b7-8021-63fa29c80b34
# ╠═f84b6daf-ade9-4ef0-b0d9-88a93d27ea1a
# ╠═d55c0cbe-8463-4d97-b75c-614a902728ea
# ╠═ddbc75b5-9dc4-4fb4-81e7-be2aafb556cd
# ╠═d5021fd3-064d-4d56-ae7c-a2a13b468f05
# ╠═b4b13c32-9572-4411-90e7-e7ff4667737b
# ╠═3e9e7c0c-643a-4ad9-bbb8-d02a09d0fde8
# ╟─e1519269-475b-409f-9f34-a45769877f5b
# ╠═488cbe38-c106-4d19-ad34-363d763a4a15
# ╠═9edb5562-3fc0-4abe-aea6-d2dd859ecf5c
# ╠═f5eb20a0-6468-411f-8569-42969c8cf2be
# ╠═5f91864e-8469-4b78-92c9-80337c23eea5
# ╠═c523ce28-f026-4b67-bde2-d0dac3e9eeec
# ╠═79c24ade-65b1-4ac9-8162-4bed95bb85bf
# ╠═be6e4906-e560-4258-98f9-a273fea2bba9
# ╠═835d31ae-03a7-4839-84e2-d5e21e3a15cd
# ╠═5873ffa5-5962-471f-a579-32470d555f7a
# ╠═1f8a655e-a914-468b-98df-f9a3cfd19171
# ╠═668a9da8-41b7-4931-9637-c38e299265cd
# ╠═0661ea23-a9ce-4530-aa19-9040d4296c25
# ╠═3561ab12-574f-4a57-b60c-432532df49e5
# ╠═5d0da605-8cc5-46be-91d3-e004c31d05bc
# ╠═c449e275-3353-4296-8a2f-7f8bfa6132c9
# ╠═df67d1b6-71bf-4792-80e3-6eac11312e24
# ╠═2bf3d695-ffed-4745-9277-d4668b3d7ae0
# ╠═5edb6c01-66fe-4c18-a8ee-6078e8255fe8
