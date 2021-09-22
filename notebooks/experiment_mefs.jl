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

# ╔═╡ f3a6e0d0-f3c0-11eb-3aee-4545a96915ba
begin
	using Pkg; Pkg.activate()
	using Random
	using Convex, ECOS
	using Plots
	using PlutoUI
	using JLD
	using LinearAlgebra
	
	OPT = () -> ECOS.Optimizer(verbose=false)
	δ = 0.01
end;

# ╔═╡ b0b551d5-c484-4499-b307-7c809d444a8a
Pkg.add("SparsityDetection")

# ╔═╡ 373f8fdd-2363-4f03-be52-974d1fd781a7
Pkg.add("SparseDiffTools")

# ╔═╡ a16793d7-0725-4bc5-bb5c-b7c44fd59375
begin
	using Revise
	using CarbonNetworks
end

# ╔═╡ 57c6f8cc-5fc3-4738-aa9f-7552d9ebe417
using BenchmarkTools

# ╔═╡ cfc87b04-1bbf-496c-8a71-5edfc9d6399c
using CarbonNetworks: kkt, compute_jacobian_kkt

# ╔═╡ c0174666-9813-4ecf-ba65-281b2f4dc0aa
using Zygote

# ╔═╡ d0d9a342-dc26-4957-961c-e50a0f8f440a
using ForwardDiff

# ╔═╡ 454392a0-8ab7-4ee4-a0fb-fa7e8a20f8de
using SparsityDetection

# ╔═╡ 939f931c-fb72-4cc8-81cf-87ed9941f2ae
using SparseDiffTools

# ╔═╡ ef235d66-30e9-4b88-9bc5-440f88aa948c
using SparseArrays

# ╔═╡ 76b72e95-e8fd-4c49-94fc-72ed94f5ba10
theme(:default, label=nothing, 
		tickfont=(:Times, 8), guidefont=(:Times, 8), legendfont=(:Times, 8))

# ╔═╡ 15dc8f16-49f9-40c5-bfa4-2cbc1d892e52
md"""
## Generate data
"""

# ╔═╡ 50d152c6-937f-41f8-868e-83fbc46a5a94
begin
	net, d_peak = load_synthetic_network("case118.m")
	n = length(d_peak)
	ℓ = size(net.B, 2)
	
	# Eliminate transmission
	net.pmax .= sum(net.gmax)
	net.fq .= 0
	
	"Network loaded."
end

# ╔═╡ 0abc0612-c92c-49aa-8967-48a09c6bcf89
demand_data = load_demand_data("2021_07_01", normalize_rows=true);

# ╔═╡ b015b97c-3fef-47ee-b36a-dab80ff26008
T_day, n_demand = size(demand_data);

# ╔═╡ 78c059d7-2891-438f-ab7a-dac9bd0780b7
begin
	Random.seed!(1)
	demand_profiles = rand(1:n_demand, n)
	
	d_dyn = hcat([
		d_peak .* demand_data[t, demand_profiles]
		for t in 1:T_day
	]...)
	
	d_total = sum(d_dyn, dims=2)
end;

# ╔═╡ 762d6c7e-ef48-4c56-9c07-3237e815bc47
size(d_dyn)

# ╔═╡ 2ccde2fb-489c-4436-8876-87cbe26881b8
num_days = 3*28

# ╔═╡ 392accfa-b8be-483d-9f47-35cf2a45f87c
ϵ_daily = 0.10

# ╔═╡ 6e189371-6f72-4f93-8fb2-179c8fade636
begin
	Random.seed!(2)
	d = hcat([
		(1 - ϵ_daily .+ ϵ_daily*rand(n)) .* d_dyn for _ in 1:num_days
	]...)
end;

# ╔═╡ a2465906-4342-4106-aa90-6ac8a4568489
T = size(d, 2)

# ╔═╡ f272e6f5-9475-44e4-8d67-2df690329cb6
times = (1/T_day) : (1/T_day) : num_days;

# ╔═╡ 584d88a5-4e25-4a48-9e68-bb0596a624a1
# Emissions rates
b = [0.8e3, 0.4e3, 0.9e3, 1.4e3, 1.8e3, 1.8e3];

# ╔═╡ fda9fc17-cf64-4737-b2fc-9e9aa53b5cc3
begin
	g = zeros(ℓ, T)
	opfs = []
	
	for t in 1:T
		# Solve OPF problem
		opf_t = PowerManagementProblem(net, d[:, t])
		solve!(opf_t, OPT)
		@assert opf_t.problem.status in [Convex.MOI.OPTIMAL, Convex.MOI.ALMOST_OPTIMAL] "t=$t, status=$(opf_t.problem.status)"
		
		# Get local emissions
		g[:, t] = evaluate(opf_t.g)
		push!(opfs, opf_t)
	end	
end

# ╔═╡ b29037ad-6f1a-4e7c-8c9e-f6619484aa49
md"""
## Compute MEFs via regression
"""

# ╔═╡ 6723e9df-8d34-4230-a69a-30620b740eff
E = reshape(sum(b .* g, dims=1), :)

# ╔═╡ bf2773a5-8355-4207-81b0-b2099c175490
plot(times, g', lw=2, xlabel="t", ylabel="generation", size=(650, 200))

# ╔═╡ aeb11efe-77f7-4afb-a6a7-64b3b5ddad7c
plot(times, E, lw=2, xlabel="t", ylabel="generation", size=(650, 200))

# ╔═╡ 6bcc5093-3aa3-43e4-b20a-051e0aeb18de
agg_demand = sum(d, dims=1)[1, :]

# ╔═╡ 14638a4b-ed45-48b4-add3-fa2a7406fb0a
ΔD = diff(agg_demand)

# ╔═╡ 95555f94-8d09-42e1-9195-0efb09ba1893
ΔE = diff(E)

# ╔═╡ 9ac37885-2495-4e96-b3b4-8464e1018a56
μ_regression = ΔD \ ΔE

# ╔═╡ 163b893e-cba0-40b6-ad4f-4d25723e46fb
begin
	plt = plot(size=(250, 250), xlabel="Δmw", ylabel="Δ(lbs co2 / h)")
	scatter!(ΔD, ΔE)
	plot!(sort(ΔD), μ_regression * sort(ΔD), lw=3, alpha=0.9)
end

# ╔═╡ c5535139-26ba-4380-bd95-fe1f2f3ef654
md"""
## Compute MEFs via differentiation
"""

# ╔═╡ db589efb-1e3f-45bc-83e5-77443bcddb59
@bind diff_case Slider(1:(T-1))

# ╔═╡ 2128bf24-9be6-43c0-8754-b7cddb844ab2
diff_case

# ╔═╡ 2c63d8ea-5d84-4b2e-94b7-3418ff466d31
ordering = sortperm(ΔD);

# ╔═╡ 44be20a2-3c13-438c-9cd3-5d8bb6f542d7
μ_diff = first(compute_mefs(opfs[ordering .+ 1][diff_case], net, d[:, ordering .+ 1][:, diff_case], b))

# ╔═╡ 84ef4c02-c754-4300-a52b-b212c866be3b
begin
	plt2 = plot(size=(350, 200), xlabel="Δmw", ylabel="Δ(lbs co2 / h)")
	scatter!(ΔD, ΔE)
	scatter!([ΔD[ordering][diff_case]], [ΔE[ordering][diff_case]], ms=6)
	plot!(sort(ΔD), μ_diff * sort(ΔD), lw=3, alpha=0.9)
end

# ╔═╡ 5cde7ef1-f9e5-49a1-84b0-43124ddc3bd9
md"""
## Analyze results
"""

# ╔═╡ 3434dc78-bab4-4882-b5b2-f69ca035c0e0
P = opfs[1];

# ╔═╡ 253d5e18-8b49-4a3a-ace2-e8a87b0f29dc
x = flatten_variables(P);

# ╔═╡ f87a29d0-bd2e-4a1c-9dcb-859e866a5487
d_test = d[:, 1];

# ╔═╡ 0e25643a-a67b-47b3-82a5-47f848e8849e
kkt_op(x) = kkt(x, net, d_test)

# ╔═╡ b490ca51-2acc-471f-8a62-1e125cccdd6c
output_buffer = zeros(size(ForwardDiff.jacobian(kkt_op, x)));

# ╔═╡ 4ebfb393-b630-414b-8e46-cbc77dbc096e
cfg = ForwardDiff.JacobianConfig(kkt_op, x);

# ╔═╡ 37242e45-70f0-4978-93da-43997510d810
md"""
### Timing
"""

# ╔═╡ 789dfcf2-54f9-4316-96fd-580d98249154
md"""
##### Solver
"""

# ╔═╡ 006e5b5a-8446-4804-891a-05f33dbe978b
@benchmark solve!($P, OPT)

# ╔═╡ 382948f7-69e9-454b-89a5-f1ff991115a7
md"""
##### Zygote
"""

# ╔═╡ 15177ca0-9e14-4be5-b845-14979d09c301
@benchmark Zygote.forward_jacobian(kkt_op, $x);

# ╔═╡ f30f3b03-99ea-49fa-ba04-6b160d71b497
md"""
##### ForwardDiff
"""

# ╔═╡ ae944a0a-0ac2-4b6e-a45f-27083675ac86
grad_fd(x) = ForwardDiff.jacobian(kkt_op, x, cfg, Val{false}());

# ╔═╡ f1f28dca-f739-4c4d-b103-186ac19be4fe
@benchmark grad_fd($x)

# ╔═╡ 91b43131-a5a0-44b6-9c69-73e104b9e5df
md"""
##### ForwardDiff (in-place)
"""

# ╔═╡ a34d8280-6673-48da-b9c6-bf7c7f262915
grad_fd!(out, x) = ForwardDiff.jacobian!(out, kkt_op, x)

# ╔═╡ 582becdd-10e7-4e89-8cb3-28235ff4cb1b
@benchmark grad_fd!($output_buffer, $x)

# ╔═╡ e2aa1ab3-ef06-427f-b958-5fc28d8a22b3
md"""
##### Manual Jacobian
"""

# ╔═╡ 209e64c1-c47e-43d8-a63c-21fcc19dfb30
grad_man(x) = compute_jacobian_kkt(net.fq, net.fl, d_test, net.pmax, net.gmax, net.A, net.B, x)

# ╔═╡ 9b2da7b6-a5e6-441c-8488-89a8efdc0064
@benchmark grad_man($x)

# ╔═╡ 8f1e0d9a-016b-408e-a408-7c612a4dd2c4
md"""
##### Sparsity Tools
"""

# ╔═╡ ea52b823-b72a-40ae-be9b-d6151ba70c77
function kkt_op!(y, x)
	y .= kkt_op(x)
	return y
end

# ╔═╡ 938ecac8-ccd9-4f86-b14c-a190f4fe2514
y = kkt_op(x);

# ╔═╡ c81dee2c-37a1-4cd3-b9da-4edd7b3187be
# sparsity_pattern = jacobian_sparsity(kkt_op!, y, x);

# ╔═╡ 18a9a913-8f20-406e-8383-cc4e243f84f7
md"""
##### ForwardDiff for each block
"""

# ╔═╡ d3d1ffc2-dd29-429c-9910-272daa025854
function lagrangian(x, net, d)
	fq, fl, A, B, pmax, gmax = net.fq, net.fl, net.A, net.B, net.pmax, net.gmax
	n, m = size(A)
    n, l = size(B)
	
	g, p, λpl, λpu, λgl, λgu, ν = unflatten_variables(x, n, m, l)
	
	return lagrangian(g, p, λpl, λpu, λgl, λgu, ν, net, d)
end

# ╔═╡ bd1b2a8c-3ada-42b2-be21-09e3207c9cb0
function lagrangian(g, p, λpl, λpu, λgl, λgu, ν, net, d)
	fq, fl, A, B, pmax, gmax = net.fq, net.fl, net.A, net.B, net.pmax, net.gmax

	return (
		(1/2)*g'*(fq .* g) + fl'g
		+ λpl' * (-p - pmax) + λpu' * (p - pmax)
		+ λgl' * (-g) + λgu' * (g - gmax)
		+ ν' * (A*p - B*g + d)
	)
end

# ╔═╡ 3b1ce0cd-2dbb-405d-887c-8cbfeab047ae
function grad_lagrangian(g, p, λpl, λpu, λgl, λgu, ν, net, d)
	fq, fl, A, B, pmax, gmax = net.fq, net.fl, net.A, net.B, net.pmax, net.gmax
	return [fq .* g + fl - B'ν - λgl + λgu; A'ν + λpu - λpl]
end

# ╔═╡ 0f4e9a29-8945-47ce-bc73-56119469ee1a
function inequality(g, p, net, d)
	fq, fl, A, B, pmax, gmax = net.fq, net.fl, net.A, net.B, net.pmax, net.gmax
	return [-p - pmax; p - pmax; -g; g - gmax]
end

# ╔═╡ 8b2cc7b0-6f25-4c63-8735-0c0822a9c9b0
function equality(g, p, net, d)
	fq, fl, A, B, pmax, gmax = net.fq, net.fl, net.A, net.B, net.pmax, net.gmax
	return A*p - B*g + d
end

# ╔═╡ cebfb82b-a3f2-4220-97e4-2ae25e6ee434
function jacobian_semimanual(x, net, d)
	g, p, λpl, λpu, λgl, λgu, ν = unflatten_variables(x, n, m, ℓ)
	
	J_g = ForwardDiff.jacobian(
		g -> grad_lagrangian(g, p, λpl, λpu, λgl, λgu, ν, net, d_test),
		g
	)
	J_p = ForwardDiff.jacobian(
		p -> grad_lagrangian(g, p, λpl, λpu, λgl, λgu, ν, net, d_test),
		p
	)
	J11 = [J_g J_p]
	

# ╔═╡ 2dc6a92b-bb1b-44d3-ac58-841a50309ba8
md"""
TODO: The challenge here is that we need to write things in terms of three variables: x, λ, and ν!
"""

# ╔═╡ 61df8c87-91d3-4fd2-ba5e-7cad2aa1d63b
md"""
## Apendix A: Making sure the Jacobians agree
"""

# ╔═╡ 45dbb744-21b4-4fcb-a74c-01ac5fb30ef3
J_man = compute_jacobian_kkt(net.fq, net.fl, d_test, net.pmax, net.gmax, net.A, net.B, x);

# ╔═╡ cbd745cc-5f1e-4eb5-934f-dd4d5fe95f9a
jac = Float64.(sparse(J_man))  #Float64.(sparse(sparsity_pattern));

# ╔═╡ 62388672-c0cc-4eff-a948-a9ed6f9dd0c1
colors = matrix_colors(jac);

# ╔═╡ 503f44a3-e798-4ccf-bed5-256448f2c3f8
maximum(colors)

# ╔═╡ dfe3f7d1-57f2-4e35-8c39-530606d9ac57
J_sparse = forwarddiff_color_jacobian!(jac, kkt_op!, x, colorvec = colors);

# ╔═╡ fc1305a1-09e1-4e9e-83bc-755e5650a26b
@benchmark forwarddiff_color_jacobian!(jac, kkt_op!, x, colorvec = colors)

# ╔═╡ 56b9af8d-774b-4aa7-9351-8986257d2ced
jac

# ╔═╡ 26624f7f-878b-48a8-a9e7-05d650b5628c
norm(J_sparse - J_man)

# ╔═╡ 585e662e-a1d2-436c-9332-d0845157cbee
sparse(J_man)

# ╔═╡ 02554976-765c-4306-ac09-7bdda637f921
J_fd = ForwardDiff.jacobian(kkt_op, x);

# ╔═╡ 4aafd4e9-bc06-4d1f-b334-aeadc994e959
sparse(J_fd)

# ╔═╡ 4ef65e6e-6cb4-46fe-9077-7eee660c423f
J_zyg = Zygote.forward_jacobian(kkt_op, x)[2];

# ╔═╡ 7725e477-cd0b-4374-a08e-e099aab8953b
norm(J_fd - J_zyg')

# ╔═╡ 9659a65a-6b91-4e0d-b3d6-9bf6c170e8c3
norm(J_fd - J_man)

# ╔═╡ 0bdaf75e-5f94-439f-a60f-c79ee7be2c84
resid = J_fd - Matrix(J_man);

# ╔═╡ 48c8a0a8-a386-4188-8482-9ce9c7e83468
m = size(net.A, 2)

# ╔═╡ f3fddd28-6ae7-42e4-a253-4fcce5c8e670
let
	g, p, λpl, λpu, λgl, λgu, ν = unflatten_variables(x, n, m, ℓ)
	
	
	#grad_lagrangian(g, p, λpl, λpu, λgl, λgu, ν, net, d_test)
	J_g = ForwardDiff.jacobian(
		g -> grad_lagrangian(g, p, λpl, λpu, λgl, λgu, ν, net, d_test),
		g
	)
	
	
end;

# ╔═╡ c7248266-be03-47b9-9afe-680228b3bca8
n, ℓ, m

# ╔═╡ 82624ece-7b63-4ffd-a92f-fd45fb0fb7c7
size(J_fd)

# ╔═╡ 38cde096-558e-446a-bdb1-e7ac8fc65cc0
3m + 3ℓ + n

# ╔═╡ e76f1a9d-ba63-46c4-9818-49d2f1ab5054
block1, block2 = (3m+2ℓ+1):(3m+3ℓ), 1:2ℓ

# ╔═╡ 8e9b9cf7-2f31-42c6-9868-100093ef13f4
plot(
	heatmap(abs.(J_fd)[block1, block2], yflip=true, clim=(0, 0.001), colorbar=false),
	heatmap(abs.(resid)[block1, block2], yflip=true, clim=(0, 0.001), colorbar=false)
);

# ╔═╡ 533027ae-df7a-4c93-822c-df378b87ef99
function my_compute_jacobian_kkt(fq, fl, d, pmax, gmax, A, B, x; τ=0.0)
    n, m = size(A)
    n, l = size(B)

    g, p, λpl, λpu, λgl, λgu, _ = CarbonNetworks.unflatten_variables(x, n, m, l)

    K11 = [
        Diagonal(fq)     spzeros(l, m);
        spzeros(m, l)    τ * I(m)
    ]
   
    K12 = [
        spzeros(l, 2 * m)   -I(l)       I(l);
        -I(m)             I(m)        spzeros(m, 2l)
    ]

    K13 = [-B'; A']

    K21 = [
        spzeros(m, l) -Diagonal(λpl);
        spzeros(m, l) Diagonal(λpu);
        -Diagonal(λgl) spzeros(l, m);
        Diagonal(λgu) spzeros(l, m)
    ]

    K22 = Diagonal([-p - pmax; p - pmax; -g; g - gmax])
    
    return [
        K11 K12 K13;
        K21 K22 spzeros(2 * (m + l), n);
        K13' spzeros(n, 2 * (m + l) + n)
    ]
end

# ╔═╡ Cell order:
# ╠═f3a6e0d0-f3c0-11eb-3aee-4545a96915ba
# ╠═a16793d7-0725-4bc5-bb5c-b7c44fd59375
# ╠═76b72e95-e8fd-4c49-94fc-72ed94f5ba10
# ╟─15dc8f16-49f9-40c5-bfa4-2cbc1d892e52
# ╠═50d152c6-937f-41f8-868e-83fbc46a5a94
# ╠═0abc0612-c92c-49aa-8967-48a09c6bcf89
# ╠═b015b97c-3fef-47ee-b36a-dab80ff26008
# ╠═78c059d7-2891-438f-ab7a-dac9bd0780b7
# ╠═762d6c7e-ef48-4c56-9c07-3237e815bc47
# ╠═2ccde2fb-489c-4436-8876-87cbe26881b8
# ╠═392accfa-b8be-483d-9f47-35cf2a45f87c
# ╠═6e189371-6f72-4f93-8fb2-179c8fade636
# ╠═a2465906-4342-4106-aa90-6ac8a4568489
# ╠═f272e6f5-9475-44e4-8d67-2df690329cb6
# ╠═584d88a5-4e25-4a48-9e68-bb0596a624a1
# ╠═fda9fc17-cf64-4737-b2fc-9e9aa53b5cc3
# ╟─b29037ad-6f1a-4e7c-8c9e-f6619484aa49
# ╠═6723e9df-8d34-4230-a69a-30620b740eff
# ╠═bf2773a5-8355-4207-81b0-b2099c175490
# ╠═aeb11efe-77f7-4afb-a6a7-64b3b5ddad7c
# ╠═6bcc5093-3aa3-43e4-b20a-051e0aeb18de
# ╠═14638a4b-ed45-48b4-add3-fa2a7406fb0a
# ╠═95555f94-8d09-42e1-9195-0efb09ba1893
# ╠═9ac37885-2495-4e96-b3b4-8464e1018a56
# ╠═163b893e-cba0-40b6-ad4f-4d25723e46fb
# ╟─c5535139-26ba-4380-bd95-fe1f2f3ef654
# ╠═db589efb-1e3f-45bc-83e5-77443bcddb59
# ╟─2128bf24-9be6-43c0-8754-b7cddb844ab2
# ╠═2c63d8ea-5d84-4b2e-94b7-3418ff466d31
# ╠═44be20a2-3c13-438c-9cd3-5d8bb6f542d7
# ╠═84ef4c02-c754-4300-a52b-b212c866be3b
# ╟─5cde7ef1-f9e5-49a1-84b0-43124ddc3bd9
# ╠═57c6f8cc-5fc3-4738-aa9f-7552d9ebe417
# ╠═cfc87b04-1bbf-496c-8a71-5edfc9d6399c
# ╠═c0174666-9813-4ecf-ba65-281b2f4dc0aa
# ╠═d0d9a342-dc26-4957-961c-e50a0f8f440a
# ╠═3434dc78-bab4-4882-b5b2-f69ca035c0e0
# ╠═253d5e18-8b49-4a3a-ace2-e8a87b0f29dc
# ╠═f87a29d0-bd2e-4a1c-9dcb-859e866a5487
# ╠═0e25643a-a67b-47b3-82a5-47f848e8849e
# ╠═b490ca51-2acc-471f-8a62-1e125cccdd6c
# ╠═4ebfb393-b630-414b-8e46-cbc77dbc096e
# ╟─37242e45-70f0-4978-93da-43997510d810
# ╟─789dfcf2-54f9-4316-96fd-580d98249154
# ╠═006e5b5a-8446-4804-891a-05f33dbe978b
# ╟─382948f7-69e9-454b-89a5-f1ff991115a7
# ╠═15177ca0-9e14-4be5-b845-14979d09c301
# ╟─f30f3b03-99ea-49fa-ba04-6b160d71b497
# ╠═ae944a0a-0ac2-4b6e-a45f-27083675ac86
# ╠═f1f28dca-f739-4c4d-b103-186ac19be4fe
# ╟─91b43131-a5a0-44b6-9c69-73e104b9e5df
# ╠═a34d8280-6673-48da-b9c6-bf7c7f262915
# ╠═582becdd-10e7-4e89-8cb3-28235ff4cb1b
# ╟─e2aa1ab3-ef06-427f-b958-5fc28d8a22b3
# ╠═209e64c1-c47e-43d8-a63c-21fcc19dfb30
# ╠═9b2da7b6-a5e6-441c-8488-89a8efdc0064
# ╟─8f1e0d9a-016b-408e-a408-7c612a4dd2c4
# ╠═b0b551d5-c484-4499-b307-7c809d444a8a
# ╠═373f8fdd-2363-4f03-be52-974d1fd781a7
# ╠═454392a0-8ab7-4ee4-a0fb-fa7e8a20f8de
# ╠═ea52b823-b72a-40ae-be9b-d6151ba70c77
# ╠═938ecac8-ccd9-4f86-b14c-a190f4fe2514
# ╠═c81dee2c-37a1-4cd3-b9da-4edd7b3187be
# ╠═cbd745cc-5f1e-4eb5-934f-dd4d5fe95f9a
# ╠═939f931c-fb72-4cc8-81cf-87ed9941f2ae
# ╠═62388672-c0cc-4eff-a948-a9ed6f9dd0c1
# ╠═503f44a3-e798-4ccf-bed5-256448f2c3f8
# ╠═dfe3f7d1-57f2-4e35-8c39-530606d9ac57
# ╠═fc1305a1-09e1-4e9e-83bc-755e5650a26b
# ╠═26624f7f-878b-48a8-a9e7-05d650b5628c
# ╠═56b9af8d-774b-4aa7-9351-8986257d2ced
# ╠═585e662e-a1d2-436c-9332-d0845157cbee
# ╠═4aafd4e9-bc06-4d1f-b334-aeadc994e959
# ╟─18a9a913-8f20-406e-8383-cc4e243f84f7
# ╠═d3d1ffc2-dd29-429c-9910-272daa025854
# ╠═bd1b2a8c-3ada-42b2-be21-09e3207c9cb0
# ╠═3b1ce0cd-2dbb-405d-887c-8cbfeab047ae
# ╠═0f4e9a29-8945-47ce-bc73-56119469ee1a
# ╠═8b2cc7b0-6f25-4c63-8735-0c0822a9c9b0
# ╠═f3fddd28-6ae7-42e4-a253-4fcce5c8e670
# ╠═cebfb82b-a3f2-4220-97e4-2ae25e6ee434
# ╠═2dc6a92b-bb1b-44d3-ac58-841a50309ba8
# ╟─61df8c87-91d3-4fd2-ba5e-7cad2aa1d63b
# ╠═45dbb744-21b4-4fcb-a74c-01ac5fb30ef3
# ╠═02554976-765c-4306-ac09-7bdda637f921
# ╠═4ef65e6e-6cb4-46fe-9077-7eee660c423f
# ╠═7725e477-cd0b-4374-a08e-e099aab8953b
# ╠═9659a65a-6b91-4e0d-b3d6-9bf6c170e8c3
# ╠═0bdaf75e-5f94-439f-a60f-c79ee7be2c84
# ╠═48c8a0a8-a386-4188-8482-9ce9c7e83468
# ╠═c7248266-be03-47b9-9afe-680228b3bca8
# ╠═82624ece-7b63-4ffd-a92f-fd45fb0fb7c7
# ╠═38cde096-558e-446a-bdb1-e7ac8fc65cc0
# ╠═e76f1a9d-ba63-46c4-9818-49d2f1ab5054
# ╠═8e9b9cf7-2f31-42c6-9868-100093ef13f4
# ╠═ef235d66-30e9-4b88-9bc5-440f88aa948c
# ╠═533027ae-df7a-4c93-822c-df378b87ef99
