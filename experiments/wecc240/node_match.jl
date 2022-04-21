### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 87226227-885e-4e0c-a0c0-59f9027a8e08
using Dates

# ╔═╡ fc5263ff-e382-4f42-bad1-df958cbf553b
using LinearAlgebra

# ╔═╡ 0446fb81-0a2c-42ec-8e88-9fe77caef65d
using DataFrames

# ╔═╡ aba69cce-c0f7-11ec-0796-ffd0516acfcd
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
end;

# ╔═╡ 188c1787-1374-4604-9c8a-0a6b9bae6ed7
begin
	util = ingredients("util.jl")
end

# ╔═╡ e990d01b-59c9-461b-bc3f-8270f65b0500
md"""## Loading 2004 nodes"""

# ╔═╡ 1b468574-c563-466a-b7f3-83e558d0950b
begin
	# 2004 dataset
	df = util.load_wecc_240_dataset()
	# Network structure
	node_names_04, node_ids_04, nicknames_04 = util.get_node_info(df.branch)
end;

# ╔═╡ 320d0bbb-298d-449c-b725-736f1ce0d52f
length(node_ids_04), length(unique(node_ids_04))

# ╔═╡ 688edfa7-a878-4129-99d0-00f567da2cda
sort(node_ids_04)

# ╔═╡ 88cdbfb3-e891-40ce-921b-59886b44554d
length(node_names_04), length(unique(node_names_04))

# ╔═╡ 420b28e7-6be5-4643-a042-4b353c554d45
unique(node_names_04)

# ╔═╡ 2ed35de2-1dd4-47ff-a316-6c9d6ce0c74f
md"""## Loading 2018 nodes"""

# ╔═╡ b90d7ba1-8090-481f-94d8-af34bbd06d34
begin
	# 2018 dataset
	DATE_FORMAT = "yy-mm-dd-HH"
	test_date = "2018-01-01-01"
	date = DateTime(test_date, DATE_FORMAT)
end

# ╔═╡ b218b2d7-65e8-4931-b2c8-d456ba3b59ec
begin
	params = util.get_nrel_data(date)
	(; name, lat, lon) = params.node
end

# ╔═╡ 1522792a-ba57-4ece-aacd-6280d79a6b30
params.node.name

# ╔═╡ 0a10ce03-7cd0-4186-95e0-9eec661d7f52
length(params.node.name), length(unique(params.node.name))

# ╔═╡ 0683c6f1-6afb-4ec4-a3c3-bea511dbff2b
params.node.id_map

# ╔═╡ 594124be-ffe6-461d-9479-d04ecad34b5b
collect(keys(params.node.id_map)) # keys of the dict

# ╔═╡ fafedee0-2080-4df6-a413-9217e3ae3a58
collect(values(params.node.id_map))

# ╔═╡ 0644b7f9-8280-44fb-a1fe-d046e5b6e77c
params.node.nickname

# ╔═╡ ec94cfb9-51d6-4289-b142-0dfa41c2403b
length(params.node.nickname), length(unique(params.node.nickname))

# ╔═╡ f880f159-595a-4d58-a92c-b96ae5d69288
md"""
## Building a matching dataframe
"""

# ╔═╡ c87e4809-da9d-4671-b092-5a749a49eb72
md"""
It seems that the node_id from 004 matches the node name from 2018. 

let us check if the nicknames match
"""

# ╔═╡ da43a648-6308-4c95-b4b6-823801e4ae71
I(5)

# ╔═╡ 8d6ce266-f068-4be3-acdd-5e288cd7edb4
begin
map_nodes, nickname_04 = get_map_nodes(node_ids_04, params.node.name, nicknames_04)
end

# ╔═╡ 3554d619-22a8-47c8-bb92-e1743cc022da
sum(sum(map_nodes, dims=2).==0)

# ╔═╡ 67336d03-be80-48ce-9410-a1d54daaad3e
size(map_nodes)

# ╔═╡ 090e4c52-b7eb-41d2-ac39-9347a6bcf2fa
df_name = DataFrame(Number=params.node.name, Name04=nickname_04, Name18=params.node.nickname);

# ╔═╡ 9f18bc61-6e03-4370-9e92-b7e73bc2e789
df_name.Equal = df_name.Name04 .== df_name.Name18

# ╔═╡ c63bb825-557e-4592-bb15-65de6dee3369
df_name

# ╔═╡ db7f2a23-1fba-44a0-b88e-f97894f4bd49
md""" We see that 240/243 nodes are matching the name exactly:"""

# ╔═╡ ec05ddf6-2105-4791-be13-16cc541f6e0f
sum(skipmissing(df_name.Equal))

# ╔═╡ 277274a0-0e7b-4571-bee5-eecd75e05d36
md"""And that only 3 are missing:"""

# ╔═╡ 15b6fa3d-c5e3-4eb3-a6de-3f830b8deda9
idx_missing = findall(ismissing, df_name.Equal)

# ╔═╡ bd22b4ff-2fb9-41ee-b62f-34de2f9bda21
df_name[idx_missing, :]

# ╔═╡ 3b117431-c016-4c19-b20a-b6f5e1df31dd
md"""
We have the mapping between the nodes from 04 and 18. There is only three nodes missing. The only thing to be done is to map the node_ids_04 to the params.nodes.name

Looking at the code, it seems that the generator data is matched via node_ids, already (which are the matching elements here). So there is not much to do I think: 
- load the network parameters (topology, constraints, ... ) from NREL
- handle the `missing` possibility in `get_generator_data` for 2004
"""

# ╔═╡ 47656f59-56c8-4313-b5ff-46807e4e1e3c
sum(map_nodes[idx_missing, :], dims=2) # those are the nodes that don't have matching with the map!

# ╔═╡ 4d96331b-b89b-4f84-bc9c-f09c97c10e9d
md"""
## How to integrate that in the code? 
"""

# ╔═╡ 437ab621-b63e-4af2-b469-aa086e131e54
begin
	year=2004
	month=1
	day=1
	hour=1
	demand_map = util.get_demand_map(hour, day, month, year, df.demand)
    d = util.make_demand_vector(demand_map, node_names_04, df.participation)
end

# ╔═╡ 893863dd-1e75-4227-85ed-df587ccd44db
demand_map

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[compat]
DataFrames = "~1.3.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "b153278a25dd42c65abbf4e62344f9d22e59191b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.43.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "6c19003824cbebd804a51211fd3bbd81bf1ecad5"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.3"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "28ef6c7ce353f0b35d0df0d5930e0d072c1f5b9b"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.1"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═aba69cce-c0f7-11ec-0796-ffd0516acfcd
# ╠═188c1787-1374-4604-9c8a-0a6b9bae6ed7
# ╟─e990d01b-59c9-461b-bc3f-8270f65b0500
# ╠═1b468574-c563-466a-b7f3-83e558d0950b
# ╠═320d0bbb-298d-449c-b725-736f1ce0d52f
# ╠═688edfa7-a878-4129-99d0-00f567da2cda
# ╠═88cdbfb3-e891-40ce-921b-59886b44554d
# ╠═420b28e7-6be5-4643-a042-4b353c554d45
# ╟─2ed35de2-1dd4-47ff-a316-6c9d6ce0c74f
# ╠═87226227-885e-4e0c-a0c0-59f9027a8e08
# ╠═b90d7ba1-8090-481f-94d8-af34bbd06d34
# ╠═b218b2d7-65e8-4931-b2c8-d456ba3b59ec
# ╠═1522792a-ba57-4ece-aacd-6280d79a6b30
# ╠═0a10ce03-7cd0-4186-95e0-9eec661d7f52
# ╠═0683c6f1-6afb-4ec4-a3c3-bea511dbff2b
# ╠═594124be-ffe6-461d-9479-d04ecad34b5b
# ╠═fafedee0-2080-4df6-a413-9217e3ae3a58
# ╠═0644b7f9-8280-44fb-a1fe-d046e5b6e77c
# ╠═ec94cfb9-51d6-4289-b142-0dfa41c2403b
# ╟─f880f159-595a-4d58-a92c-b96ae5d69288
# ╟─c87e4809-da9d-4671-b092-5a749a49eb72
# ╠═fc5263ff-e382-4f42-bad1-df958cbf553b
# ╠═da43a648-6308-4c95-b4b6-823801e4ae71
# ╠═8d6ce266-f068-4be3-acdd-5e288cd7edb4
# ╠═3554d619-22a8-47c8-bb92-e1743cc022da
# ╠═67336d03-be80-48ce-9410-a1d54daaad3e
# ╠═0446fb81-0a2c-42ec-8e88-9fe77caef65d
# ╠═090e4c52-b7eb-41d2-ac39-9347a6bcf2fa
# ╠═9f18bc61-6e03-4370-9e92-b7e73bc2e789
# ╠═c63bb825-557e-4592-bb15-65de6dee3369
# ╟─db7f2a23-1fba-44a0-b88e-f97894f4bd49
# ╠═ec05ddf6-2105-4791-be13-16cc541f6e0f
# ╟─277274a0-0e7b-4571-bee5-eecd75e05d36
# ╠═15b6fa3d-c5e3-4eb3-a6de-3f830b8deda9
# ╠═bd22b4ff-2fb9-41ee-b62f-34de2f9bda21
# ╟─3b117431-c016-4c19-b20a-b6f5e1df31dd
# ╠═47656f59-56c8-4313-b5ff-46807e4e1e3c
# ╠═4d96331b-b89b-4f84-bc9c-f09c97c10e9d
# ╠═437ab621-b63e-4af2-b469-aa086e131e54
# ╠═893863dd-1e75-4227-85ed-df587ccd44db
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
