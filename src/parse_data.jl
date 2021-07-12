get_data_path(datapath, f) = joinpath(datapath, f)

SKIP_RESOURCES = ["WND", "SUN", "BIO", "GEO", "OTH", "WAT"]
BASE_RESOURCES = ["COL", "NG", "NUC", "OIL"]
DEMAND = "DEMAND"


function parse_network_data(datapath; num_generators=1)
    # Open dataframes
    df_node, df_branch, df_resource = open_datasets(datapath)
    
    num_resource = length(BASE_RESOURCES)
    n = nrow(df_node)
    l = n * num_resource * num_generators
    # l' = n * num_resource

    # Preprocess
    df_node = coalesce.(df_node, 0.0)  # Fill missing values with 0.0

    # Initialize outputs
    nodes = []
    branches = []
    gmax = Float64[]
    gen_labels = []
    pmax = Float64[]
    f = Float64[] #Carbon cost

    # Create node generation map
    generator_per_node = num_generators * num_resource 
    B = hcat([I(n) for _ in 1:generator_per_node]...)
    @assert size(B) == (n, l)
    
    # Parse nodes
    nodes = df_node[:, :id]
    for res in BASE_RESOURCES
        for i in 1:num_generators
            for (iso_ind, iso) in enumerate(eachrow(df_node))
                capacity = iso[res * "_max"] / num_generators
                push!(gmax, capacity)
                ind = findfirst(r -> r.id == res, eachrow(df_resource))
                push!(f, df_resource[ind, :emission_factor])
                push!(gen_labels, (iso_ind, res))
            end
        end
    end

    # Parse cross-ISO branches
    for row in eachrow(df_branch)
        if row.sink_node == "ALL"
            continue
        end
        i = findfirst(x -> x == row.source_node, nodes)
        j = findfirst(x -> x == row.sink_node, nodes)
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
    
    return A, B, gmax, pmax, f, G, nodes, gen_labels
end

function open_datasets(datapath)
    df_node = DataFrame(CSV.File(get_data_path(datapath, "node_data.csv")))
    sort!(df_node, [:id])
    df_branch = DataFrame(CSV.File(get_data_path(datapath, "branch_data.csv"), 
        drop=(i, name)->(i==1)))
    df_resource = DataFrame(CSV.File(get_data_path(datapath, "resource_data.csv"), 
        header=[:id, :name, :emission_factor], datarow=2))
    return df_node, df_branch, df_resource
end

function create_generation_map(gen_labels)
    l = length(gen_labels)  #n * num_resource * num_generators

    agg_generators = unique(gen_labels)
    lp = length(agg_generators)

    M = zeros(lp, l)
    for i in 1:l
       idx = findfirst(x -> x == gen_labels[i], agg_generators)
       M[idx, i] = 1.0
    end

    return agg_generators, M
end

function load_case(name, B, agg_generators; resources=BASE_RESOURCES)
    n, l = size(B)
    
    case = DataFrame(CSV.File(name))
    case = rename!(case, "Column1" => "ba")
    case = coalesce.(case, 0.0)

    bas = case[:, "ba"]
    demands = case[:, "demand"]
    for res in SKIP_RESOURCES
        demands -= case[:, res]
    end
 
    d = demands
    
    # Load generation profile
    g = zeros(n*length(resources))
    for (node, res) in product(1:n, resources)
        idx = findfirst(x -> x == (node, res), agg_generators)
        g[idx] = case[node, res]
    end
    
    return d, g, case
end