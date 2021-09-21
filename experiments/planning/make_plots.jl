### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ ac4ef27c-7562-455e-8ed5-397d57634afd
using Pkg; Pkg.activate("")

# ╔═╡ ae6537c7-b714-402a-b113-606e7556b63a
using Plots

# ╔═╡ 5ef99391-cc00-4767-a91b-eb539cd5ad02
theme(:default,
	label=nothing,
	tickfont=(:Times, 8),
	guidefont=(:Times, 8),
	legendfont=(:Times, 8),
	titlefont=(:Times, 8, :bold),
	dpi=100,
)

# ╔═╡ ae663449-d931-4c3d-a689-6298e2cddc2a
begin
	in_to_pix(x) = 100*x
	in_to_pix(x, y) = (100*x, 100*y)
end

# ╔═╡ 649ce882-ba32-4f68-9c36-9d9dc454064e
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ 0853bb22-1636-11ec-25b3-3350b9d91e8c
begin
	util = ingredients("util.jl")
	
	load_results = util.load_results
	solve_problem = util.formulate_and_solve_problem
end

# ╔═╡ 1c020762-2035-46b9-8456-18391a973c70
md"""
## Figure 1: Convergence plot
"""

# ╔═╡ 2a444a9f-6530-481d-adb6-81c85b6a05e8
P, config, θ, history = util.load_results("base");

# ╔═╡ 25a090d1-acdf-4501-aa5e-807fd8c19f34
begin
	seed_results = []
	for i in 1:length(util.SEED_CONFIG)
		_, _, _, history = load_results("seed$(i)")
		θ = history.θ[argmin(history.loss)-1]
		push!(seed_results, (θ=θ, history=history))
	end
end

# ╔═╡ 2e9c26bd-771c-4d4f-9343-4d37a439c8c4
ξ = P.ξ

# ╔═╡ 183ed6d5-884b-446d-acc2-e75f2d86abe9
md"""
## Figure 2: Pareto curve + carbon tax
"""

# ╔═╡ 3d1e4c4e-538b-40a2-bbc9-adafa9490d94
begin
	tax_results = []
	for i in 1:length(util.CARBON_TAX_CONFIG)
		P, _, _, history = load_results("tax$(i)")
		θ = history.θ[end]
		push!(tax_results, (θ=θ, history=history, P=P))
	end
end

# ╔═╡ fa3bf131-14a1-4afb-a1ef-93d32eb59404
begin
	Ms_tax = []
	Es_tax = []
	for r in tax_results
		J, x, _, _, _ = solve_problem(r.θ, r.P)
		
		push!(Ms_tax, J + r.P.ξ'*(r.θ-r.P.θ_min))
		push!(Es_tax, sum([r.P.c'xt for xt in x]))
	end
end

# ╔═╡ 7f529f04-a945-4b88-a646-79f88e534683
J_base, x_base, _, _, _ = solve_problem(θ, P);

# ╔═╡ be3cb9fc-8b95-4287-b089-252a378067b1
M_base = J_base + P.ξ'*(θ-P.θ_min)

# ╔═╡ d0ba32cf-58a8-42fe-b3bd-a8fbe8ba084c
E_base = sum([P.c'xt for xt in x_base])

# ╔═╡ 58e66183-8859-429c-8dd1-377f17c57571
seed = 5

# ╔═╡ 961c7c6a-f73c-49d1-891b-18ffc3ea1e0e
begin
	dim_lam, dim_seed = size(util.PARETO_CONFIG)
	
	pareto_results = []
	
	for i in 1:dim_lam
		_, _, _, history = load_results("pareto$(i)_$(seed)")
		θ = history.θ[argmin(history.loss)-1]
		push!(pareto_results, (θ=θ, history=history))
	end
	
	pareto_results[1] = tax_results[1]
end

# ╔═╡ ced988c3-e644-4bfa-85a4-cf6ab075b0f6
begin
	Ms = []
	Es = []
	for r in pareto_results[[1:7; 11:21]]
		J, x, Dx, dnet, opf = solve_problem(r.θ, P)
		
		push!(Ms, J + P.ξ'*(r.θ-P.θ_min))
		push!(Es, sum([P.c'xt for xt in x]))
	end
end

# ╔═╡ b3021867-a100-470d-8cd2-04e005ae9cdf
net_gen = (
	sum(sum(P.d_dyn))  # Demand per day
	* P.Z  # Time horizon
)

# ╔═╡ 5019aabc-525f-4c50-a205-8826c28b0300
k = 6

# ╔═╡ f07689f4-f861-45f6-a40d-a3e1a7368fa3
md"""
## Figure 2: Network changes
With, and without, tax
"""

# ╔═╡ 20092653-2b80-4bc5-a1d9-d26930a4e9c1
function make_rect(z)
	verts = [(0, 0), (0, z), (1, z), (1, 0), (0, 0)]
	return Shape(verts)
end

# ╔═╡ 274b4c81-46d1-4c95-9f4d-c7be0d2b1389
make_font(s) = text(s, :Times, 8)

# ╔═╡ 6297c39b-3ea3-435c-83a9-e8d2b634c2ad
n, m = size(P.net.A)

# ╔═╡ cee5b0b4-0ac1-40aa-8fbb-ca7942b23a78
positions = [
	(0, 0),
	(1, -5),
	(7, -5.5),
	(6.5, -2),
	(2, -1),
	(2.5, 0.5),
	(7.5, -0.5),
	(8, 1),
	(7, 2),
	(5, 2),
	(3.5, 3),
	(1, 3.5),
	(2.5, 4.5),
	(5, 4),
]

# ╔═╡ dc074b9a-3021-4140-bdaa-1d7f9d78ab1b
x_pos, y_pos = collect.(zip(positions...))

# ╔═╡ 6c624b96-bb83-4318-bb5c-9e8f1c297517
# begin
# 	k = 3
# 	yticks3 = 0 : 0.5 : 1.5
	
# 	θA = θ .* (θ .> 0.1)
# 	θB = tax_results[k].θ .* (tax_results[k].θ .> 0.1)
	
# 	plt3 = plot(
# 		bar(θA[1:m] .- P.net.pmax, 
# 			yticks=(yticks3, yticks3), xticks=false, 
# 			label="Penalty", title="Transmission", 
# 			ylim=(0, 1.55), legend=(0.75, 0.8),
# 			yguide="Added Capacity [MW]", color=:Green,
# 		),
# 		bar(θA[m+1:end], 
# 			xticks=false, ylim=(0, 1.55), title="Storage", yticks=false, 
# 			color=:Green,
# 		),
# 		bar(θB[1:m] .- P.net.pmax, 
# 			c=:Red, label="Tax", ylim=(0, 1.55), yticks=(yticks3, yticks3),
# 			legend=(0.75, 0.8), xticks=false, xlabel="Line",
# 			yguide="Added Capacity [MW]"
# 		), 
# 		bar(θB[m+1:end], 
# 			c=:Red, ylim=(0, 1.55), yticks=false, xticks=false, xlabel="Bus"
# 		),
# 		layout=(2, 2), grid=false, frame=:box
# 	)
# 	plot!(size=in_to_pix(0.6*textwidth, figheight), link=:y)
# 	annotate!(plt3[1], -2, 1.6, text("C", :Times, :left, :bottom, 14))
# end

# ╔═╡ 2fda99db-6bf2-4a90-82b6-2006a5dad386
# function make_graph(θ, legend=false)
# 	C = θ[m+1:end]
# 	C_rel = ((C .- minimum(C)) / maximum(C .- minimum(C)))
	
	
# 	plt = plot()
	
# 	# Plot edges
# 	edge_lw = sqrt.(P.net.pmax)
# 	edge_lw .*= 4 / maximum(edge_lw)
# 	edge_lw .+= 2
# 	for col in 1:size(P.net.A, 2)
# 		fr = findfirst(x -> x == -1, P.net.A[:, col])
# 		to = findfirst(x -> x == 1, P.net.A[:, col])
# 		lw = sqrt(P.net.pmax[col])
		
# 		plot!(x_pos[[fr, to]], y_pos[[fr, to]], c=:Gray, lw=edge_lw[col], alpha=0.5)
# 	end
	
# 	# Plot nodes with labels
# 	is_gen = (P.d_dyn[1] .== 0) .+ 1
# 	node_c = [:Black, :Yellow][is_gen]
# 	scatter!(x_pos, y_pos, shape=:rect, ms=4, group=is_gen, c=node_c, markerborder=10, label=["load" "generator"])
# 	#?annotate!(x_pos .- 0.5, y_pos .+ 0.5, make_font.(["$i" for i in 1:n]))
	
# 	# Plot changes to storage
# 	cap_grad = cgrad([:honeydew2, :green])
# 	scatter!(x_pos .- 0.7, y_pos .- 0.25, 
# 		shape=make_rect.(2*C), ms=4, c=cap_grad[C_rel], line_z=cap_grad[C_rel], seriescolor=cap_grad)
	
# 	# Adjust plot
# 	plot!(xlim=(minimum(x_pos)-1, maximum(x_pos)+1))
# 	plot!(ylim=(minimum(y_pos)-1, maximum(y_pos)+1))
# 	plot!(xticks=nothing, yticks=nothing)
# 	plot!(frame=:none, legend=legend)
# 	plot!(size=in_to_pix(0.3*textwidth, 3), margin=0Plots.pt, right_margin=0Plots.pt)
	
# 	return plt
# end

# ╔═╡ 9cef002c-12e3-408c-95ff-13c40fae6583
textheight = 556.47656 / 72

# ╔═╡ 9b7bf97c-4895-4a1b-b97a-359dbd1b4678
textwidth = 430.00462 / 72

# ╔═╡ 21790d04-851a-4f9f-acd0-dfde51f4c7fe
figheight = 2.5

# ╔═╡ 215d236b-0b90-4c5a-8151-98efc29dff73
begin
	xticks = 0 : 200 : 400
	yticks = 26 : 2 : 32
	
	plt1 = plot()  #plot(history.loss / 1e3, lw=3, alpha=0.75)
	for i in [1, 3, 5]
		plot!(seed_results[i].history.loss / 1e3, lw=3, alpha=0.75, color=1)
	end
	
	plot!(xlabel="Iteration", ylabel="Loss")
	plot!(xticks=(xticks, xticks), yticks=(yticks, yticks), frame=:box, grid=false)
	plot!(size=in_to_pix(0.4textwidth, figheight/2))
	plot!(xlim=(0, 500), ylim=(24, 33))
	annotate!(0, 24, text("A", :Times, :left, :bottom, 14))
end

# ╔═╡ 577a49a5-6c18-4176-b624-9d0bdf094fa4
begin
	yticks2 = 35 : 5 : 45
	
	plt2 = plot()
	
	plot!(Ms/net_gen, P.Z*Es, lw=3, label="Penalty", color=1, alpha=0.7)
	plot!(Ms_tax/net_gen, P.Z*Es_tax, lw=3, label="Tax", color=1, ls=:dot)
	
	scatter!([M_base/net_gen], P.Z*[E_base], color=:Green, ms=5)
	
	# $125 and $500 tax
	scatter!([Ms_tax[k]/net_gen], P.Z*[Es_tax[k]], color=:Red, ms=5)
	scatter!([Ms_tax[end]/net_gen], P.Z*[Es_tax[end]], color=:Yellow, ms=5)
	
	plot!(xlabel="Levelized Cost [\$/MWh]", ylabel="Emissions [mTCO2]")
	plot!(frame=:box, grid=false)
	plot!(xlim=(20, 150), ylim=(30, 50))
	plot!(size=in_to_pix(0.4*textwidth, figheight/2), bottom_margin=8Plots.pt)
	plot!(legend=(0.8, 1))
	
	
	annotate!(20, 30, text("B", :Times, :left, :bottom, 14))
end

# ╔═╡ dfe68663-7c4e-48c6-94c0-1d2e27ecb515
begin
	yticks3 = 0 : 0.5 : 1.5
	bw = 0.4
	
	θA = θ .* (θ - P.θ_min .> 0.1)
	θB = tax_results[k].θ .* (tax_results[k].θ - P.θ_min .> 0.1)
	
	
	plt3_top = plot(
		xlabel="Line", xlim=(0, m+1), xticks=[1, 2, 4, 5, 6, 7, 8, 11, 12, 13, 18, 19],
		yticks=(yticks3, yticks3), ylim=(0, 1.8), ylabel="Transmission [MW]",
		legend=(0.75, 0.8),	
	)
	bar!(1:m, θA[1:m] .- P.net.pmax, 
		bar_width=bw, color=:Green, label="Penalty")
	bar!((1:m) .+ bw, θB[1:m] .- P.net.pmax, 
		bar_width=bw, c=:Red, label="Tax")
	annotate!(0.25, 1.5, text("C", :Times, :left, :top, 14))
	
	
	plt3_bot = plot(
		xlabel="Bus", xlim=(0, n+1), xticks=[1, 2, 4, 7, 8, 9, 10, 11, 14],
		ylim=(0, 1.55), ylabel="Storage [MWh]",
		
	)
	bar!(1:n, θA[m+1:end], bar_width=bw, color=:Green)
	bar!((1:n) .+ bw, θB[m+1:end], bar_width=bw, c=:Red)
	annotate!(0.25, 1.5, text("D", :Times, :left, :top, 14))
	
	
	
	plt3 = plot(
		plt3_bot, plt3_top,
		layout=(2, 1), grid=false, frame=:box
	)
	
	plot!(plt3[2], bottom_margin=8Plots.pt)
	plot!(size=in_to_pix(0.6*textwidth, figheight), link=:y)
	
end

# ╔═╡ b45669ed-f4c8-42fd-9f78-4d98806884ea
begin
	l = Plots.@layout [
		[a{0.5h}; b] c{0.6w}
	]
	
	main_plt = plot(
		plt1, plt2, plt3, 
		#deepcopy(plt3), 
		layout=l,
		size=in_to_pix(textwidth, figheight),
	)
	
	savefig("../../img/planning_figure.pdf")
	main_plt
end

# ╔═╡ Cell order:
# ╠═ac4ef27c-7562-455e-8ed5-397d57634afd
# ╠═ae6537c7-b714-402a-b113-606e7556b63a
# ╠═5ef99391-cc00-4767-a91b-eb539cd5ad02
# ╠═ae663449-d931-4c3d-a689-6298e2cddc2a
# ╟─649ce882-ba32-4f68-9c36-9d9dc454064e
# ╠═0853bb22-1636-11ec-25b3-3350b9d91e8c
# ╟─1c020762-2035-46b9-8456-18391a973c70
# ╠═2a444a9f-6530-481d-adb6-81c85b6a05e8
# ╠═25a090d1-acdf-4501-aa5e-807fd8c19f34
# ╠═2e9c26bd-771c-4d4f-9343-4d37a439c8c4
# ╠═215d236b-0b90-4c5a-8151-98efc29dff73
# ╠═183ed6d5-884b-446d-acc2-e75f2d86abe9
# ╠═3d1e4c4e-538b-40a2-bbc9-adafa9490d94
# ╠═fa3bf131-14a1-4afb-a1ef-93d32eb59404
# ╠═7f529f04-a945-4b88-a646-79f88e534683
# ╠═be3cb9fc-8b95-4287-b089-252a378067b1
# ╠═d0ba32cf-58a8-42fe-b3bd-a8fbe8ba084c
# ╠═58e66183-8859-429c-8dd1-377f17c57571
# ╠═961c7c6a-f73c-49d1-891b-18ffc3ea1e0e
# ╠═ced988c3-e644-4bfa-85a4-cf6ab075b0f6
# ╠═b3021867-a100-470d-8cd2-04e005ae9cdf
# ╠═5019aabc-525f-4c50-a205-8826c28b0300
# ╠═577a49a5-6c18-4176-b624-9d0bdf094fa4
# ╟─f07689f4-f861-45f6-a40d-a3e1a7368fa3
# ╟─20092653-2b80-4bc5-a1d9-d26930a4e9c1
# ╟─274b4c81-46d1-4c95-9f4d-c7be0d2b1389
# ╠═6297c39b-3ea3-435c-83a9-e8d2b634c2ad
# ╟─cee5b0b4-0ac1-40aa-8fbb-ca7942b23a78
# ╠═dc074b9a-3021-4140-bdaa-1d7f9d78ab1b
# ╠═dfe68663-7c4e-48c6-94c0-1d2e27ecb515
# ╟─6c624b96-bb83-4318-bb5c-9e8f1c297517
# ╟─2fda99db-6bf2-4a90-82b6-2006a5dad386
# ╠═9cef002c-12e3-408c-95ff-13c40fae6583
# ╠═9b7bf97c-4895-4a1b-b97a-359dbd1b4678
# ╠═21790d04-851a-4f9f-acd0-dfde51f4c7fe
# ╠═b45669ed-f4c8-42fd-9f78-4d98806884ea
