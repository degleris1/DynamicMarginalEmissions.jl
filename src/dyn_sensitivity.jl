# Computing sensitivities for dynamic power management problems

"""
    sensitivity_demand_dyn(P::PowerManagementProblem, net::DynamicPowerNetwork, d, ∇C, t)

Compute `∇_{d_t} C( x_opt(d_t) )`, where `x_opt(d_t)` is the optimal solution (primal
and dual variables) of the power management problem `P` with parameters 
`(fq, fl, d, pmax, gmax, A, B)` and demand `d_t` at time `t`, and where
`∇C` is the gradient `∇_x C(x)`.
"""
function sensitivity_demand_dyn(P::PowerManagementProblem, net::DynamicPowerNetwork, d, ∇C)

    _, _, l, ns, T = get_problem_dims(net)
    x = flatten_variables_dyn(P)

    # Get partial Jacobians of KKT operator
    ∂K_xT = sparse(adjoint(compute_jacobian_kkt_dyn(x, net, d)))

    v = ∂K_xT \ ∇C


    if norm(∂K_xT*v - ∇C) / length(v) > 1e-6
        @warn "KKT Jacobian is ill-conditioned. Solutions may be innaccurate."
        @show norm(∂K_xT*v - ∇C), length(v)
    end

    ∇C_θ = []
    for t in 1:T
        _, ∂K_θT = Zygote.forward_jacobian(
            dt -> kkt_dyn(x, net, [tp == t ? dt : d[tp] for tp in 1:T]), 
            d[t]
        )
        _, ncols = size(∂K_θT)
        to_delete = _to_delete(ncols, ns, l)
        ∂K_θT_ = prune_matrix(∂K_θT, to_delete, axis="cols")

        # Now compute ∇C(g*(θ)) = -∂K_θ' * inv(∂K_x') * v
        push!(∇C_θ, -∂K_θT_ * v)
    end

    return ∇C_θ
end


"""
    compute_mefs(P::PowerManagementProblem, net::DynamicPowerNetwork, d, c)

Compute the marginal emission factors given carbon costs `c`
and demands `d` across the entire time horizon.
"""
function compute_mefs(P::PowerManagementProblem, net::DynamicPowerNetwork, d, c)
    ∇C_dyn = _make_∇C(net, c)

    return sensitivity_demand_dyn(P, net, d, ∇C_dyn)
end

"""
    _make_∇C(net::DynamicPowerNetwork, d, c)

Constructs carbon cost gradient to be propagated for mef computation for 
the DynamicPowerNetwork `net`. `d` is the demand and `c` are the carbon costs. 
"""
function _make_∇C(net::DynamicPowerNetwork, c, cq=0, g=0)
    # Extract dimensions
    n, m, l, ns, T = get_problem_dims(net)
    static_dim = kkt_dims(n, m, l)
 
    # Construct ∇_x C(x)
    ∇C_dyn = zeros(kkt_dims_dyn(n, m, l, ns, T) - 3ns, T)
    idx = 0
    for t in 1:T
        if g == 0
            ∇C_dyn[idx+1 : idx+l, t] .= c
        else
            ∇C_dyn[idx+1 : idx+l, t] .= c .+ (cq .* g[t])
        end
         idx += static_dim
    end

    # Return sensitivity
    return ∇C_dyn
end

"""
    compute_jacobian_kkt_dyn(x, net, d_dyn)

Constructs the Jacobian for a dynamic network `net`, with input variables `x` and
demand `d_dyn`.

Notes
-----
The Jacobian in this system is computed as... 
# TODO: add
"""
function compute_jacobian_kkt_dyn(x, net, d_dyn)

    n, m, l, ns, T = get_problem_dims(net)

    dim_t = kkt_dims(n, m, l)
    dim_s = storage_kkt_dims(ns, l)

    # decompose `x` in arrays of T variables
    g, p, s, ch, dis, λpl, λpu, λgl, λgu, ν, νE, λsl, λsu, λchl, λchu, λdisl, λdisu, λrampl, λrampu, νs = 
        unflatten_variables_dyn(x, n, m, l, ns, T)

    # Compute individual Jacobians for the static system
    Kτ1 = [
        compute_jacobian_kkt(
            net.fq[t], net.fl[t], d_dyn[t], net.pmax[t], net.gmax[t], net.A, net.B, net.F,
            [g[t]; p[t]; λpl[t]; λpu[t]; λgl[t]; λgu[t]; ν[t]; νE[t]]; τ=TAU)
        for t in 1:T
    ]
    Kτ2 = [
        compute_jacobian_kkt_charge_discharge_ramp(dim_t, ns, m, l, net.F, net.S)
        for _ in 1:T
    ]
    Kτ3 = [
        compute_jacobian_kkt_future_ramp(dim_t, ns, l)
        for _ in 1:(T-1)
    ]

    # Stack matrices
    Kτ = [
        (t == T) ?
        hcat(
            spzeros(dim_t, (t-1)*dim_t),
            Kτ1[t],
            spzeros(dim_t, (T-t)*dim_t),
            spzeros(dim_t, (t-1)*dim_s),
            Kτ2[t],
            spzeros(dim_t, (T-t)*dim_s)
        ) :
        hcat(
            spzeros(dim_t, (t-1)*dim_t),
            Kτ1[t],
            spzeros(dim_t, (T-t)*dim_t),
            spzeros(dim_t, (t-1)*dim_s),
            Kτ2[t],
            Kτ3[t],
            spzeros(dim_t, (T-t-1)*dim_s)
        )
        for t in 1:T
    ]

    # Compute individual Jacobians for the dynamic system
    g_prev = t -> (t == 1) ? zeros(l) : g[t-1]
    ρ_ = t -> (t==1) ? net.gmax[1] : net.ρ

    # we do not need the below I think because of unflatten....()
    # s_crt = t -> (t==T) ? FINAL_COND.*C : s[t]
    Ks = [
        compute_jacobian_kkt_dyn_t(
            s[t], ch[t], dis[t], 
            λsl[t], λsu[t], λchl[t], λchu[t], λdisl[t], λdisu[t], λrampl[t], λrampu[t], 
            g[t], g_prev(t), ρ_(t),
            net, t
        )
        for t in 1:T
    ]
    
    ∂K = [vcat(Kτ...) ; vcat(Ks...)] # standard form of the Jacobian

    # eliminate the rows and columns that correspond to either s[T] or the dual variables
    # for the inequality constraints
    _, ncols = size(∂K)
    to_delete = _to_delete(ncols, ns, l)
    ∂K = prune_matrix(∂K, to_delete)

    return ∂K
end

"""
Returns the rows/columns to eliminate in the jacobian
"""
function _to_delete(size, ns, l)
    idx_starts = [
        size - storage_kkt_dims(ns, l), 
        size - storage_kkt_dims(ns, l) + 3ns, 
        size - storage_kkt_dims(ns, l) + 4ns 
        ]
    to_delete = [] 
    for idx in idx_starts
        for i in 1:ns
            push!(to_delete, idx+i)
        end
    end
    return to_delete
end

"""
Eliminate rows and columns in a matrix
"""
function prune_matrix(∂K, to_delete; axis="both")
    _, ncols = size(∂K)
    to_keep = [k for k in 1:ncols if !(k in to_delete)]
    if axis=="both" || axis=="rows"
        ∂K = ∂K[to_keep, :]
    end
    if axis=="both" || axis=="cols"
        ∂K = ∂K[:, to_keep]
    end
    return ∂K
end


"""
    compute_jacobian_kkt_charge_discharge(dims, n)

Compute the part of the Jacobian of the static problem associated with dynamic variables (primal and dual), 
with `dims` being the dimension of the static system, `ns` the number of storage nodes, and `l` 
the number of generators.
"""
function compute_jacobian_kkt_charge_discharge_ramp(dims, ns, m, l, F, S)
    dKdch = [spzeros(dims-m-1, ns); F*S; -ones(1, ns)]
    dKddis = [spzeros(dims-m-1, ns); -F*S; ones(1, ns)]
    dKdλl = [-I(l); spzeros(dims-l, l)]
    dKdλu = [I(l); spzeros(dims-l, l)]

    # each entry in the return vector corresponds to a primal and then dual variable
    return [spzeros(dims, ns) dKdch dKddis spzeros(dims, 6ns) dKdλl dKdλu spzeros(dims, ns)]
end


"""
    compute_jacobian_kkt_future_ramp(dims, ns, l)

Compute the part of the Jacobian (dStatic / dRamp), where `dims` is the 
dimension of the static system, `n` is the number of nodes, and `l` is
the number of generators.

Note: this function computes it only for the future timestep. 
The ramping constraints associated with the current timestep are accounted for 
in `compute_jacobian_kkt_charge_discharge_ramp`.

"""
function compute_jacobian_kkt_future_ramp(dims, ns, l)
    dKdλl = [I(l); spzeros(dims-l, l)]
    dKdλu = [-I(l); spzeros(dims-l, l)]

    return [spzeros(dims, 9ns) dKdλl dKdλu spzeros(dims, ns)]
end

"""
    compute_jacobian_kkt_dyn_t(s, ch, dis, λsl, λsu, λchl, λchu, λdisl, λdisu, net, t)

Compute the Jacobian for a given timestep of the storage part of the problem. 
Input variables are the primal and dual variables, `net` is the dynamic network and 
`t` is the timestep at which the Jacobian is to be computed.
"""
function compute_jacobian_kkt_dyn_t(
    s, ch, dis, 
    λsl, λsu, λchl, λchu, λdisl, λdisu, λrampl, λrampu, 
    gt, gt_prev, ρ,
    net, t
)

    # extract some variables
    F = net.F
    η_c = net.η_c
    η_d = net.η_d
    C = net.C
    P = net.P
    S = net.S
    n, m, l, ns, T = get_problem_dims(net)

    # Relevant sizes for the problem
    kdims = kkt_dims(n, m, l)
    sdims = storage_kkt_dims(ns, l)

    ##########################
    #
    # 1. Constructing the dynamic part of the KKTs
    #
    ##########################

    # D_x ∇_x L 
    # dim = 3n, 3n
    K11 = spzeros(3ns, 3ns) # dim = ns for s, ch, and dis 
    # D_x F^T
    # dim = 3ns, (6ns+2l)
    K12 = [
            -I(ns) I(ns) spzeros(ns, 4ns+2l);
            spzeros(ns, 2ns) -I(ns) I(ns) spzeros(ns, 2ns+2l);
            spzeros(ns, 4ns) -I(ns) I(ns) spzeros(ns, 2l)
        ]
    # D_x H^T
    # dim = 3ns, ns
    K13 = [
        -I(ns);
        η_c*I(ns);
        -(1/η_d)*I(ns)
    ]
    # Diag(λ) * D_xF
    # dim = 8ns, 8ns
    D_xF = K12'
    D = Diagonal([λsl; λsu; λchl; λchu; λdisl; λdisu; λrampl; λrampu])
    K21 = D * D_xF
    
    # diag(F)
    # dim = 8ns, 8ns
    K22 = Diagonal([-s; s-C; -ch; ch-P; -dis; dis-P; (gt_prev - ρ - gt); (gt - gt_prev - ρ)])
    
    # diagonal part of the Jacobian
    K_storage_t = [
        K11 K12 K13;
        K21 K22 spzeros(6ns+2l, ns);
        K13' spzeros(ns, 7ns+2l) 
    ]

    ##########################
    #
    # 2. Constructing the non-diagonal parts of the Jacobian
    #    i.e. includes cross-terms with the static variables
    #       - with the ν dual variable
    #       - with λgl, λgu variables
    #
    ##########################

    Dλramp = (
        (t == 1) ? 
        [
            spzeros(l, (t-1)*kdims) -Diagonal(λrampl) spzeros(l, kdims-l+(T-t)*kdims);
            spzeros(l, (t-1)*kdims) Diagonal(λrampu) spzeros(l, kdims-l+(T-t)*kdims)
        ] : 
        [
            spzeros(l, (t-2)*kdims) Diagonal(λrampl) spzeros(l, kdims-l) -Diagonal(λrampl) spzeros(l, kdims-l+(T-t)*kdims);
            spzeros(l, (t-2)*kdims) -Diagonal(λrampu) spzeros(l, kdims-l) Diagonal(λrampu) spzeros(l, kdims-l+(T-t)*kdims)
        ]
    )

    # The KKT of the static problem associated with the dynamic variables, 
    # namely: s, ch, dis, ramp
    K_static = [
        spzeros(ns, T*kdims); # s
        spzeros(ns, t*kdims - m - 1) -(F*S)' ones(ns) spzeros(ns, (T-t)*kdims); # ch
        spzeros(ns, t*kdims - m - 1) (F*S)' -ones(ns) spzeros(ns, (T-t)*kdims); # dis
        spzeros(6ns, T*kdims); # dual variables : in dynamic problem
        Dλramp;
        spzeros(ns, T*kdims); # coupling between s and ch/dis : in dynamic problem
    ]
    
    if t > 1
        K_storage_left = [
            spzeros(sdims-ns, (t-1)*sdims);
            spzeros(ns, (t-2)*sdims) I(ns) spzeros(ns, sdims-ns);
        ]
    else
        K_storage_left = spzeros(sdims, 0);
    end
    if t < T
        K_storage_right = [
            spzeros(ns, sdims-ns) I(ns) spzeros(ns, (T-(t+1))*sdims);
            spzeros(sdims-ns, (T-t)*sdims);
        ]
    else
        K_storage_right = spzeros(sdims, 0);
    end

    K_storage = [K_storage_left K_storage_t K_storage_right];

    return [K_static K_storage]
   
end





