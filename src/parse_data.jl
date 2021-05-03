using DataFrames, CSV
using LightGraphs

_DATA_PATH = ENV["CARBON_NETWORKS_DATA"]

_get_data_path(f) = joinpath(_DATA_PATH, f)

function parse_network_data()
    # Open dataframes
    df_node, df_branch, df_resource = _open_datasets()
    
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
            i = length(nodes)+1
            push!(nodes, (i, iso.id, r.id))
            push!(gmax, iso[r.id * "_max"])
            push!(f, r.emission_factor)

            push!(branches, (i, k))
            push!(pmax, Inf)
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
    G = SimpleGraph()
    add_vertices!(G, length(nodes))
    for e in branches
        add_edge!(G, e...)
    end
    A = incidence_matrix(G, oriented=true)
    
    return A, gmax, pmax, f, G, nodes
end

function _open_datasets()
    df_node = DataFrame(CSV.File(_get_data_path("node_data.csv")))
    df_branch = DataFrame(CSV.File(_get_data_path("branch_data.csv"), 
        drop=(i, name)->(i==1)))
    df_resource = DataFrame(CSV.File(_get_data_path("resource_data.csv"), 
        header=[:id, :name, :emission_factor], datarow=2))
    return df_node, df_branch, df_resource
end