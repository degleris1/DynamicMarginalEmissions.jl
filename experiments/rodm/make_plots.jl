### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ f83c0c3a-200d-11ec-2135-1fa17a0df82a
begin
	using Pkg; Pkg.activate("")
	using Convex
	using Statistics
end

# ╔═╡ 5d212716-9661-4d67-9382-4b0c775d6d85
using Plots

# ╔═╡ 58555d65-d842-423c-8380-f6d38149f7d5
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

# ╔═╡ dcb739ab-7837-4807-8df1-b471492f17a2
util = ingredients("util.jl")

# ╔═╡ 3065f8e6-eebc-41ab-8516-e3f019bcfba2
dpi = 100

# ╔═╡ 3b69d7a4-2ca3-484d-b86f-4801414fce53
theme(:default,
	label=nothing,
	tickfont=(:Times, 8),
	guidefont=(:Times, 8),
	legendfont=(:Times, 8),
	titlefont=(:Times, 8, :bold),
	dpi=dpi,
)

# ╔═╡ 6dcd2b42-ecfc-4d7c-afc9-13fda3807827
md"""
## Load data

Comparing RODM, RODM + MDT, and RODM + ramping
"""

# ╔═╡ d0a30c4e-b675-439b-8bc0-adec6126360e
dyn_config, df, data, dynamic_data, results = util.load_results_dynamic("ramp");

# ╔═╡ 98c7c6d5-8c08-42ef-9652-798f9c26137d
static_config, _, _, static_results = util.load_results_static("base");

# ╔═╡ a359e195-6121-4438-aae8-0a3b2a1f1aa5
md"""
#### Compute true (historical) MEFs
"""

# ╔═╡ 9c78666f-f03c-4b1b-9a8c-8b784eeb0b65
bar(results[25].g[36] ./ data[25*7*24].θ.gmax)

# ╔═╡ 9c886bfe-7931-4cb7-8ee5-e3ef18b3de8f
demand = df.carbon.demand

# ╔═╡ 8c22b4e7-6ecf-4188-a34b-2b8e4d5f9b58
co2 = df.carbon.co2_tot

# ╔═╡ 93a455cc-c620-463e-90da-654bb96dca8b
mef_hist = (diff(co2) ./ diff(demand))[util.TIMES]

# ╔═╡ 33d7aace-51f7-411c-bef3-79e6c672152d
valid_times = findall(x -> 0 < x < Inf, mef_hist)

# ╔═╡ 0c58311b-3001-4307-b8d3-e9365a0829c1
get_err = mef_est -> (
	mean(abs, (mef_hist - mef_est)[valid_times]) 
	/ mean(mef_hist[valid_times])
)

# ╔═╡ 3330f77a-73a3-4a8e-90b3-27824d356514
md"""
#### Load MEF estimates
"""

# ╔═╡ 002f2a27-8489-4521-8ed9-660ec6cd5760
mefs_ed = vcat(getfield.(static_results, :dq)...)[:, 1]

# ╔═╡ 94ccf5d2-fa1f-4b08-aee5-51340bad7e56
mefs_mdt = df.carbon.co2_marg[util.TIMES]

# ╔═╡ 54db5151-a50c-43a0-aa2e-f873a1de30b1
mefs_ramp = vcat([sum.(r.dq[1]) for r in results]...)

# ╔═╡ bdac9b46-6773-4546-8f03-4a4e1e3c4b88
md"""
#### Compute errors
"""

# ╔═╡ dc592388-31ab-4428-9bac-9671f12b0c5c
err_ed = get_err(mefs_ed)

# ╔═╡ e9623b71-b6de-455d-8a26-c41be6d58f8c
err_mdt = get_err(mefs_mdt)

# ╔═╡ 9a65cc16-734c-4101-b65c-cda75aa9b75e
err_ramp = get_err(mefs_ramp)

# ╔═╡ 000fd4ee-c2a7-4924-ac34-b76f0ca80fe7
md"""
## Panel A: Prediction error of each model
"""

# ╔═╡ e467b463-0228-4eae-9833-7b0afcb14bf3
begin
	plta = plot(size = dpi .* (2, 1.5))
	plot!(frame=:box)
	
	bar!(["ROD", "MDT", "RMP"], 100*[err_ed, err_mdt, err_ramp]) 
	plot!(ylabel="Relative Absolute Error", yticks=([0, 25, 50], ["0%", "25%", "50%"]))
end

# ╔═╡ Cell order:
# ╠═f83c0c3a-200d-11ec-2135-1fa17a0df82a
# ╟─58555d65-d842-423c-8380-f6d38149f7d5
# ╠═dcb739ab-7837-4807-8df1-b471492f17a2
# ╠═5d212716-9661-4d67-9382-4b0c775d6d85
# ╠═3065f8e6-eebc-41ab-8516-e3f019bcfba2
# ╠═3b69d7a4-2ca3-484d-b86f-4801414fce53
# ╟─6dcd2b42-ecfc-4d7c-afc9-13fda3807827
# ╠═d0a30c4e-b675-439b-8bc0-adec6126360e
# ╠═98c7c6d5-8c08-42ef-9652-798f9c26137d
# ╟─a359e195-6121-4438-aae8-0a3b2a1f1aa5
# ╠═9c78666f-f03c-4b1b-9a8c-8b784eeb0b65
# ╠═9c886bfe-7931-4cb7-8ee5-e3ef18b3de8f
# ╠═8c22b4e7-6ecf-4188-a34b-2b8e4d5f9b58
# ╠═93a455cc-c620-463e-90da-654bb96dca8b
# ╠═33d7aace-51f7-411c-bef3-79e6c672152d
# ╠═0c58311b-3001-4307-b8d3-e9365a0829c1
# ╟─3330f77a-73a3-4a8e-90b3-27824d356514
# ╠═002f2a27-8489-4521-8ed9-660ec6cd5760
# ╠═94ccf5d2-fa1f-4b08-aee5-51340bad7e56
# ╠═54db5151-a50c-43a0-aa2e-f873a1de30b1
# ╟─bdac9b46-6773-4546-8f03-4a4e1e3c4b88
# ╠═dc592388-31ab-4428-9bac-9671f12b0c5c
# ╠═e9623b71-b6de-455d-8a26-c41be6d58f8c
# ╠═9a65cc16-734c-4101-b65c-cda75aa9b75e
# ╟─000fd4ee-c2a7-4924-ac34-b76f0ca80fe7
# ╠═e467b463-0228-4eae-9833-7b0afcb14bf3
