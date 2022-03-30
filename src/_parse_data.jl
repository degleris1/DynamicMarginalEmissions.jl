
"""
    load_synthetic_network(case_name)

Load a synthetic network file. Returns a `Network` object `θ`, a demand
vector `d`, and the original data dictionary `data`.
"""
function load_synthetic_network(case_name; line_outage=nothing)
    # Load network data
    network_data = parse_file(joinpath(@__DIR__, "../data", case_name));
    net = make_basic_network(network_data)
    make_per_unit!(net)

    # Parse network data
    base_mva = net["baseMVA"]
    gen, load, branch = net["gen"], net["load"], net["branch"]

    # Dimensions
    l = length(net["gen"])
    n = length(net["bus"])
    m = length(net["branch"])

    # Susceptance
    β = [1 / branch[i]["br_x"] for i in string.(1:m)]
    if line_outage !== nothing
        β[line_outage] = 0
    end

    # Topology
    A = SparseMatrixCSC(calc_basic_incidence_matrix(net)')
    B = _make_B(gen, n, l)
    F = make_pfdf_matrix(A, β)

    @assert size(A) == (n, m)
    @assert size(B) == (n, l)

    # Capacities
    pmax = [branch[i]["rate_a"] for i in string.(1:m)]
    gmax = [gen[i]["mbase"] * gen[i]["pmax"] for i in string.(1:l)] / base_mva

    # Generator costs
    fq = [gen[i]["cost"][1] for i in string.(1:l)]
    fl = [gen[i]["cost"][2] for i in string.(1:l)]

    θ = PowerNetwork(fq, fl, pmax, gmax, A, B, F, TAU)

    # Demand
    d = _make_d(load, n)

    return θ, d, net
end

"""
    load_demand_data(case_name::String; source="caiso", normalize_rows=false)

Load CAISO demand data from date `case_name`. Returns a `(T, n)` matrix,
where `T` is the number of timesteps (usually 24) and `n` is the number
of nodes. 
"""
function load_demand_data(case_name::String; source="caiso", normalize_rows=false)
    # Load dataframe
    filename = join([source, "demand", case_name, ".csv"], "_", "") 
    path = joinpath(@__DIR__, "../data/", filename)
    df = DataFrame(CSV.File(path))

    # Parse data
    groups = unique(df.GROUP)
    hours = unique(df.OPR_HR)

    n, T = length(groups), length(hours)

    data = zeros(T, n)
    for gr in groups
        df_gr = filter(row -> row.GROUP == gr, df)
        data[:, gr] = sort(df_gr, :OPR_HR).MW
    end

    if normalize_rows
        data ./= maximum(data, dims=1)
    end

    return data
end

"""
    load_renewable_data(case_name::String; source="caiso", normalize_rows=false)

Load CAISO renewable data from date `case_name`. Returns a `(T, n)` matrix,
where `T` is the number of timesteps (usually 24) and `n` is the number
of nodes.
"""
function load_renewable_data(case_name::String; source="caiso", normalize_rows=false)
    # Load dataframe
    filename = join([source, "renewables", case_name, ".csv"], "_", "") 
    path = joinpath(@__DIR__, "../data/", filename)
    df = DataFrame(CSV.File(path))

    # Parse data
    groups = unique(df.GROUP)
    hours = unique(df.OPR_HR)

    n, T = length(groups), length(hours)

    data = zeros(T, n)
    labels = []
    for gr in groups
        df_gr = filter(row -> row.GROUP == gr, df)
        data[:, gr] = sort(df_gr, :OPR_HR).MW

        push!(labels, lowercase(df_gr.RENEWABLE_TYPE[1]))
    end

    if normalize_rows
        data ./= maximum(data, dims=1)
    end

    return data, labels
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











using Distributions, Random

function make_dynamic(net::PowerNetwork, T, P, C, dyn_gmax, η_c, η_d)
    fqs = [net.fq for t in 1:T]
    fls = [net.fl for t in 1:T]
    pmaxs = [net.pmax for t in 1:T]
    return DynamicPowerNetwork(fqs, fls, pmaxs, dyn_gmax, net.A, net.B, net.F, P, C, T; η_c=η_c, η_d=η_d)
end

make_dynamic(net::PowerNetwork, T, P, C, η) = make_dynamic(net, T, P, C, [net.gmax for _ in 1:T], η, η);

function generate_random_data(n, l, ns, T)
    Random.seed!(2)

    # Make graph
    if n > 3
    G = watts_strogatz(n, 3, 0.3)
    else
        G = Graph(3)
        add_edge!(G, 1, 2)
        add_edge!(G, 1, 3)
        add_edge!(G, 2, 3)
    end


    # Convert to incidence matrix
    A = incidence_matrix(G, oriented=true)
    m = size(A, 2)

    # Create generator-node mapping
    node_map = vcat(collect(1:n), rand(1:n, l-n))
    B = spzeros(n, l)
    for (gen, node) in enumerate(node_map)
        B[node, gen] = 1
    end

    # Generate carbon costs
    cq = rand(Exponential(2), l)
    cl = rand(Exponential(2), l)

    # Generate demands
    d = rand(Bernoulli(0.8), n) .* rand(Gamma(3.0, 3.0), n)

    # Generate generation and flow capacities
    gmax = rand(Gamma(4.0, 3.0), l) + (B'd)  # This is just to make the problem easier
    pmax = rand(Gamma(1.0, 0.1), m);
    cq_dyn = [cq for _ in 1:T]
    cl_dyn = [cl for _ in 1:T]
    pmax_dyn = [pmax for _ in 1:T]
    gmax_dyn = [gmax for _ in 1:T]
    d_dyn = [d for _ in 1:T]

    
    C = rand(Exponential(10), ns)
    P = 0.25 * C

    S = zeros(n, ns)
    nodes_storage = sample(1:n,  ns)
    for i in 1:ns
        S[nodes_storage[i], i] = 1
    end
    
    return A, B, cq_dyn, cl_dyn, d_dyn, gmax_dyn, pmax_dyn, P, C, S
end

function generate_network(n, l, T, ns)

    A, B, cq_dyn, cl_dyn, d_dyn, gmax_dyn, pmax_dyn, P, C, S = generate_random_data(n, l, ns, T)
    β = rand(Exponential(10), n)
    F = make_pfdf_matrix(A, β)

    net_dyn = DynamicPowerNetwork(
        cq_dyn, cl_dyn, pmax_dyn, gmax_dyn, A, B, F, S, P, C, T
    )
    return net_dyn, d_dyn
end