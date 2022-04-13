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

function generate_network(n, l, T, ns; η=1.)

    A, B, cq_dyn, cl_dyn, d_dyn, gmax_dyn, pmax_dyn, P, C, S = generate_random_data(n, l, ns, T)
    β = rand(Exponential(10), n)
    F = make_pfdf_matrix(A, β)

    net_dyn = DynamicPowerNetwork(
        cq_dyn, cl_dyn, pmax_dyn, gmax_dyn, A, B, F, S, P, C, T;
        η_c=η, η_d=η
    )
    return net_dyn, d_dyn
end