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

# ╔═╡ a16793d7-0725-4bc5-bb5c-b7c44fd59375
begin
	using Revise
	using CarbonNetworks
end

# ╔═╡ 76b72e95-e8fd-4c49-94fc-72ed94f5ba10
theme(:default, label=nothing, 
		tickfont=(:Times, 8), guidefont=(:Times, 8), legendfont=(:Times, 8))

# ╔═╡ 15dc8f16-49f9-40c5-bfa4-2cbc1d892e52
md"""
## Generate data
"""

# ╔═╡ 50d152c6-937f-41f8-868e-83fbc46a5a94
begin
	net, d_peak = load_synthetic_network("case30.m")
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
