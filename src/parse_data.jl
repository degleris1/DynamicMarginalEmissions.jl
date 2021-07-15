get_data_path(datapath, f) = joinpath(datapath, f)

SKIP_RESOURCES = ["WND", "SUN", "BIO", "GEO", "OTH", "WAT"]
BASE_RESOURCES = ["COL", "NG", "NUC", "OIL"]


function parse_network_data(datapath; num_generators=1)
    # Open dataframes
    df_node, df_branch, df_resource = open_datasets(datapath)
    
    DEMAND = "DEMAND"
    num_resource = nrow(df_resource)
    num_location = nrow(df_node)
    
    # Preprocess
    df_node = coalesce.(df_node, 0.0)  # Fill missing values with 0.0

    # Initialize outputs
    nodes = []
    branches = []
    gmax = Float64[]
    pmax = Float64[]
    f = Float64[]
    
    # Parse nodes
    for iso in eachrow(df_node)
        k = length(nodes)+1
        push!(nodes, (k, iso.id, DEMAND))
        push!(gmax, 0.0)
        push!(f, 0.0)

        for r in eachrow(df_resource)
            for _ in 1:num_generators
                (r.id in SKIP_RESOURCES) && continue  # Skip some resources
        
                i = length(nodes)+1
                capacity = iso[r.id * "_max"] / num_generators
                push!(nodes, (i, iso.id, r.id))
                push!(gmax, capacity)
                push!(f, r.emission_factor)

                push!(branches, (i, k))
                push!(pmax, Inf)
            end
        end
    end
    
    # Parse cross-ISO branches
    for row in eachrow(df_branch)
        if row.sink_node == "ALL"
            continue
        end
        i = findfirst(x -> x[2:3] == (row.source_node, DEMAND), nodes)
        j = findfirst(x -> x[2:3] == (row.sink_node, DEMAND), nodes)
        push!(branches, (i, j))
        push!(pmax, row.max_exchanges)
    end
    
    # Construct graph
    G = SimpleWeightedGraph(length(nodes))
    for ((i, j), w) in zip(branches, pmax)
        add_edge!(G, i, j, w)
    end
    A = incidence_matrix(G, oriented=true)
    pmax = [e.weight for e in edges(G)]
    
    return A, gmax, pmax, f, G, nodes
end

function open_datasets(datapath)
    df_node = DataFrame(CSV.File(get_data_path(datapath, "node_data.csv")))
    df_branch = DataFrame(CSV.File(get_data_path(datapath, "branch_data.csv"), 
        drop=(i, name)->(i==1)))
    df_resource = DataFrame(CSV.File(get_data_path(datapath, "resource_data.csv"), 
        header=[:id, :name, :emission_factor], datarow=2))
    return df_node, df_branch, df_resource
end

function create_generation_map(nodes, resources=BASE_RESOURCES)
    n = length(nodes)
    bas = unique([node[2] for node in nodes])
    
    agg_nodes = []
    for (i, (ba, res)) in enumerate(product(bas, resources))
        push!(agg_nodes, (i, ba, res))
    end
    
    agg_map = spzeros(length(agg_nodes), length(nodes))
    for (i, ba, res) in nodes
        !(res in resources) && continue
        k = findfirst(x -> x[2:3] == (ba, res), agg_nodes)
        agg_map[k, i] = 1
    end

    return agg_nodes, agg_map
end

function load_case(name, agg_nodes, B, nodes; resources=BASE_RESOURCES)
    ng, n = size(B)
    
    case = DataFrame(CSV.File(name))
    case = rename!(case, "Column1" => "ba")
    case = coalesce.(case, 0.0)

    bas = case[:, "ba"]
    demands = case[:, "demand"]
    for res in SKIP_RESOURCES
        demands -= case[:, res]
    end

    # Load demands
    d = zeros(n)
    for (dem, ba) in zip(demands, bas)
        ind_d = findfirst(x -> x[2:3] == (ba, "DEMAND"), nodes)
        d[ind_d] = dem
    end  
    
    # Load generation profile
    g = zeros(ng)
    for res in resources
        for ba in bas
            ind_g = findfirst(x -> x[2:3] == (ba, res), agg_nodes)
            g[ind_g] = case[case.ba .== ba, res][1]
        end
    end
    
    return d, g, case
end

function load_synthetic_network(case_name)
    # Load network data
    network_data = parse_file(joinpath("../data", case_name));
	net = make_basic_network(network_data)
	make_per_unit!(net)

    # Parse network data
    base_mva = net["baseMVA"]
	gen, load, branch = net["gen"], net["load"], net["branch"]
	
	# Dimensions
	l = length(net["gen"])
	n = length(net["bus"])
	m = length(net["branch"])
	
	# Topology
	A = SparseMatrixCSC(calc_basic_incidence_matrix(net)')
	B = _make_B(gen, n, l)
	
	@assert size(A) == (n, m)
	@assert size(B) == (n, l)
	
	# Capacities
	pmax = [branch[i]["rate_a"] for i in string.(1:m)]
	gmax = [gen[i]["mbase"] * gen[i]["pmax"] for i in string.(1:l)] / base_mva
	
	# Generator costs
	fq = [gen[i]["cost"][1] for i in string.(1:l)]
	fl = [gen[i]["cost"][2] for i in string.(1:l)]
	
    θ = PowerNetwork(fq, fl, pmax, gmax, A, B, TAU)

    # Demand
	d = _make_d(load, n)

    return θ, d, net
end


"""
	_make_d(load, n)
"""
function _make_d(load, n)
	n_load = length(load)
	nodes = [load[i]["load_bus"] for i in string.(1:n_load)]
	#dp = [load[i]["pd"] for i in string.(1:n_load)]
	
	d = zeros(n)
	for (ind_d, ind_n) in enumerate(nodes)
		d[ind_n] = load[string(ind_d)]["pd"]
	end
	
	return d
end

"""
	_make_B(gen, n, l)
"""
function _make_B(gen, n, l)
	nodes = [gen[i]["gen_bus"] for i in string.(1:l)]
	
	B = spzeros(n, l)
	for (ind_g, ind_n) in enumerate(nodes)
		B[ind_n, ind_g] = 1.0
	end
	
	return B
end