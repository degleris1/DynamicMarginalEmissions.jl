using CarbonNetworks
using Convex, ECOS
using SparseArrays
using Zygote
using LinearAlgebra

OPT = () -> ECOS.Optimizer(verbose=false)

n_threads = 10

η_c = .99;
η_d = .82;

T = 20;

n_vec = [4, 10, 50, 100, 200, 500];
norms = zeros(size(n_vec));
time_Zygote = zeros(size(n_vec));
time_manual = zeros(size(n_vec));

Threads.@threads for i = 1:length(n_vec)
    n = n_vec[i]
    l = n

    A, B, cq_dyn, cl_dyn, d_dyn, gmax_dyn, pmax_dyn, P, C = generate_random_data(n, l, T);

    dnet = DynamicPowerNetwork(
        cq_dyn, cl_dyn, pmax_dyn, gmax_dyn, A, B, P, C, T; η_c=η_c, η_d=η_d
    )
    dmin = DynamicPowerManagementProblem(dnet, d_dyn);

    solve!(dmin, OPT, verbose=true);

    x = flatten_variables_dyn(dmin);

    tz = @elapsed _, ∂K_xT = Zygote.forward_jacobian(x -> kkt_dyn(x, dnet, d_dyn), x);
    tm = @elapsed Jac_manual = compute_jacobian_kkt_dyn(x, dnet, d_dyn);
    
    time_Zygote[i] = tz;
    time_manual[i] = tm;

    _jac = sparse(∂K_xT');

    norms[i] = norm(_jac - Jac_manual);
end

@show norms
@show time_Zygote
@show time_manual
@show time_Zygote./time_manual