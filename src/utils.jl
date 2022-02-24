using Distributions, Random

function make_dynamic(net::PowerNetwork, T, P, C, dyn_gmax, η_c, η_d)
	fqs = [net.fq for t in 1:T]
	fls = [net.fl for t in 1:T]
	pmaxs = [net.pmax for t in 1:T]
	return DynamicPowerNetwork(fqs, fls, pmaxs, dyn_gmax, net.A, net.B, net.F, P, C, T; η_c=η_c, η_d=η_d)
end

make_dynamic(net::PowerNetwork, T, P, C, η) = make_dynamic(net, T, P, C, [net.gmax for _ in 1:T], η, η);

function generate_random_data(n, l, T)
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

    
    C = rand(Exponential(10), n)
    P = 0.25 * C
    ;
    
    return A, B, cq_dyn, cl_dyn, d_dyn, gmax_dyn, pmax_dyn, P, C
end

function generate_network(n, l, T)

    A, B, cq_dyn, cl_dyn, d_dyn, gmax_dyn, pmax_dyn, P, C = generate_random_data(n, l, T)
    β = rand(Exponential(10), n)
    F = make_pfdf_matrix(A, β)

    net_dyn = DynamicPowerNetwork(
        cq_dyn, cl_dyn, pmax_dyn, gmax_dyn, A, B, F, P, C, T
    )
    return net_dyn, d_dyn
end