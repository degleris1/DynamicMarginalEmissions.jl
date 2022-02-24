# Computing sensitivities for dynamic power management problems

"""
    sensitivity_demand_dyn(P::PowerManagementProblem, net::DynamicPowerNetwork, d, ∇C, t)

Compute `∇_{d_t} C( x_opt(d_t) )`, where `x_opt(d_t)` is the optimal solution (primal
and dual variables) of the power management problem `P` with parameters 
`(fq, fl, d, pmax, gmax, A, B)` and demand `d_t` at time `t`, and where
`∇C` is the gradient `∇_x C(x)`.
"""
function sensitivity_demand_dyn(P::PowerManagementProblem, net::DynamicPowerNetwork, d, ∇C, t)
    T = net.T
    x = flatten_variables_dyn(P)

    # Get partial Jacobians of KKT operator
    _, ∂K_xT = Zygote.forward_jacobian(x -> kkt_dyn(x, net, d), x)

    # Now compute ∇C(g*(θ)) = -∂K_θ' * inv(∂K_x') * ∇C 
    v = ∂K_xT \ ∇C
    ∇C_θ = -∂K_θT * v

    return ∇C_θ
end

"""
    sensitivity_demand_dyn(P::PowerManagementProblem, net::DynamicPowerNetwork, d, ∇C)
"""

function sensitivity_demand_dyn(P::PowerManagementProblem, net::DynamicPowerNetwork, d, ∇C)
    T = net.T
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

        # Now compute ∇C(g*(θ)) = -∂K_θ' * inv(∂K_x') * v
        push!(∇C_θ, -∂K_θT * v)
    end

    return ∇C_θ
end


"""
    compute_mefs(P::PowerManagementProblem, net::DynamicPowerNetwork, d, c, t)

Compute the marginal emission factors at time `t` given carbon costs `c`
and demands `d`.
"""
function compute_mefs(P::PowerManagementProblem, net::DynamicPowerNetwork, d, c, t)
    ∇C_dyn = _make_∇C(net, c)
    
    return sensitivity_demand_dyn(P, net, d, ∇C_dyn, t)
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
    n, m, l, T = get_problem_dims(net)
    static_dim = kkt_dims(n, m, l)
 
    # Construct ∇_x C(x)
    ∇C_dyn = zeros(kkt_dims_dyn(n, m, l, T), T)
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
"""
function compute_jacobian_kkt_dyn(x, net, d_dyn)
    n, m, l, T = get_problem_dims(net)

    dim_t = kkt_dims(n, m, l)
    dim_s = storage_kkt_dims(n, l)

    # decompose `x` in arrays of T variables
    g, p, s, ch, dis, λpl, λpu, λgl, λgu, ν, νE, λsl, λsu, λchl, λchu, λdisl, λdisu, λrampl, λrampu, νs = 
        unflatten_variables_dyn(x, n, m, l, T)

    # Compute individual Jacobians for the static system

    Kτ1 = [
        compute_jacobian_kkt(
            net.fq[t], net.fl[t], d_dyn[t], net.pmax[t], net.gmax[t], net.A, net.B, net.F,
            [g[t]; p[t]; λpl[t]; λpu[t]; λgl[t]; λgu[t]; ν[t]; νE[t]]; τ=TAU)
        for t in 1:T
    ]
    Kτ2 = [
        compute_jacobian_kkt_charge_discharge_ramp(dim_t, n, m, l, net.F)
        for t in 1:T
    ]
    Kτ3 = [
        compute_jacobian_kkt_future_ramp(dim_t, n, l)
        for t in 1:(T-1)
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
    Ks = [
        compute_jacobian_kkt_dyn_t(
            s[t], ch[t], dis[t], 
            λsl[t], λsu[t], λchl[t], λchu[t], λdisl[t], λdisu[t], λrampl[t], λrampu[t], 
            g[t], g_prev(t), net.ρ,
            net, t
        )
        for t in 1:T
    ]

    return [vcat(Kτ...) ; vcat(Ks...)]
end


"""
    compute_jacobian_kkt_charge_discharge(dims, n)

Compute the part of the Jacobian associated with charge and discharge, 
with `dims` being the dimension of the static system, `n` the number of nodes, and `l` 
the number of generators.
"""
function compute_jacobian_kkt_charge_discharge_ramp(dims, n, m, l, F)
    dKdch = [spzeros(dims-m-1, n); F; -ones(1, n)]
    dKddis = [spzeros(dims-m-1, n); -F; ones(1, n)]
    dKdλl = [-I(l); spzeros(dims-l, l)]
    dKdλu = [I(l); spzeros(dims-l, l)]

    return [spzeros(dims, n) dKdch dKddis spzeros(dims, 6n) dKdλl dKdλu spzeros(dims, n)]
end


"""
    compute_jacobian_kkt_future_ramp(dims, n, l)

Compute the part of the Jacobian (dStatic / dRamp), where `dims` is the 
dimension of the static system, `n` is the number of nodes, and `l` is
the number of generators.
"""
function compute_jacobian_kkt_future_ramp(dims, n, l)
    dKdλl = [I(l); spzeros(dims-l, l)]
    dKdλu = [-I(l); spzeros(dims-l, l)]

    return [spzeros(dims, 9n) dKdλl dKdλu spzeros(dims, n)]
end

###################################################################################
#
# Below contains functions for automated testing of the sensitivity
#
# TODO: update, and make sure they work
#
####################################################################################

"""
    compute_obj_sensitivity(P::PowerManagementProblem, net::DynamicPowerNetwork, d, t)

Compute the sensitivity of monetary costs (i.e. marginal costs) at time `t` for a given demand
profile `d`.
"""
function compute_obj_sensitivity(P::PowerManagementProblem, net::DynamicPowerNetwork, d, t)
    # Extract dimensions

    n, m, l, T = get_problem_dims(net)
    static_dim = kkt_dims(n, m, l)

    # Construct ∇_x C(x)
    ∇C_dyn = zeros(kkt_dims_dyn(n, m, l, T))
    idx = (t-1)*static_dim
    ∇C_dyn[idx+1 : idx+l] .=  net.fq[t] .* vec(P.g[t].value) + net.fl[t]

    # Return sensitivity
    return sensitivity_demand_dyn(P, net, d, ∇C_dyn, t)
end

"""
    compute_var_sensitivity(P::PowerManagementProblem, net::DynamicPowerNetwork, d, varName, unit, t)

Compute the sensitivity of a given primal variable `varName` (either 'g', 'p' or 's') wrt demand `d` at time `t`. 
The `unit` parameter defines which unit is considered (e.g. generator number, edge number, or storage node number).
"""
function compute_var_sensitivity(P::PowerManagementProblem, net::DynamicPowerNetwork, d, varName, unit, t)
    # Extract dimensions
    n, m, l, T = get_problem_dims(net)
    static_dim = kkt_dims(n, m, l) 

    # J = sensitivity_demand_dyn(P, net, d, I(kkt_dims_dyn(n, m, l, T)), t)

    id_vec = zeros(kkt_dims_dyn(n, m, l, T))
    ref_val = 0
    if varName == 'g'
        idx = (t-1)*static_dim + unit
        ref_val = P.g[t].value[unit]
    elseif varName == 'p'
        idx = (t-1)*static_dim + l + unit
        ref_val = P.p[t].value[unit]
    elseif varName == 's'
        idx = T*static_dim + (t-1) * storage_kkt_dims(n) + unit
        ref_val = P.s[t].value[unit]
    end
    id_vec[idx] = 1

    return sensitivity_demand_dyn(P, net, d, id_vec, t), ref_val
end


"""
    sensitivity_demand_check(dnet::DynamicPowerNetwork, d_dyn, node, t; npoints=10, rel_inc=1e-1)

Compares the objective value with that estimated with implicit differentiation. Estimates are made wrt demand
 at a given `node` at time `t`. 

Parameters:
- `npoints` gives the number of increments in both directions (i.e. positive and negative) of the reference value of demand
- `rel_inc` specifies how large each increment is
"""
function sensitivity_demand_check(dnet::DynamicPowerNetwork, d_dyn, node, t; npoints=10, rel_inc=1e-1)
    
    dmin = DynamicPowerManagementProblem(dnet, d_dyn);
    # solve the problem
    solve!(dmin, OPT, verbose=true);
    # obtain the optimal objective value
    obj_ref = dmin.problem.optval
    # compute sensitivity of Objective to demand
    ∂O∂d = compute_obj_sensitivity(dmin, dnet, d_dyn, t)
    ∂O∂d = ∂O∂d[node] #extract for the node of interest
    
    obj_vals = zeros(2npoints+1)
    for i in -npoints:npoints
        d_crt = [copy(d) for d in d_dyn]
        d_crt[t][node] = d_dyn[t][node] * (1+ i * rel_inc)
        dmin = DynamicPowerManagementProblem(dnet, d_crt)
        solve!(dmin, OPT, verbose=true)
        obj_vals[i+1+npoints] = dmin.problem.optval
    end
    
    # sensitivity-based estimation
    estimated_obj_vals = [obj_ref + ∂O∂d*rel_inc*d_dyn[t][node]*i for i in -npoints:npoints]
    x = [1+i*rel_inc for i in -npoints:npoints]
    
    return obj_vals, estimated_obj_vals, x
end

"""
    sensitivity_var_check(dnet::DynamicPowerNetwork, d_dyn, node, varName, unit, t; npoints=10, rel_inc=1e-1)

Compares the value of a primal variable with that estimated with implicit differentiation. Estimates are made wrt demand
at a given `node` at time `t`. 'varName' specifies which primal variable whose sensitivity is estimated, 'unit' specifies
the unit number (e.g. generator number).

Parameters:
- `npoints` gives the number of increments in both directions (i.e. positive and negative) of the reference value of demand
- `rel_inc` specifies how large each increment is
"""
function sensitivity_var_check(dnet::DynamicPowerNetwork, d_dyn, node, varName, unit, t; npoints=10, rel_inc=1e-1)

    dmin = DynamicPowerManagementProblem(dnet, d_dyn);
    # solve the problem
    solve!(dmin, OPT, verbose=true);
    # compute sensitivity of Variable to demand
    ∂V∂d, ref_val = compute_var_sensitivity(dmin, dnet, d_dyn, varName, unit, t)
    ∂V∂d = ∂V∂d[node] #extract for the node of interest
    
    opt_vals = zeros(2npoints+1)
    for i in -npoints:npoints
        d_crt = [copy(d) for d in d_dyn]
        d_crt[t][node] = d_dyn[t][node] * (1+ i * rel_inc)
        dmin = DynamicPowerManagementProblem(dnet, d_crt)
        solve!(dmin, OPT, verbose=true)
        if varName == 'g'
            opt_vals[i+1+npoints] = dmin.g[t].value[unit]
        elseif varName == 'p'
            opt_vals[i+1+npoints] = dmin.p[t].value[unit]
        elseif varName == 's'
            opt_vals[i+1+npoints] = dmin.s[t].value[unit]
        end
    end
    
    # sensitivity-based estimation
    estimated_vals = [ref_val + ∂V∂d*rel_inc*d_dyn[t][node]*i for i in -npoints:npoints]
    x = [1+i*rel_inc for i in -npoints:npoints]
    
    return opt_vals, estimated_vals, x
end



"""
    compute_jacobian_kkt_dyn_t(s, ch, dis, λsl, λsu, λchl, λchu, λdisl, λdisu, net, t)

Compute the Jacobian for a given timestep of the storage part of the problem. Input variables are the primal and dual variables, 
`net` is the dynamic network and `t` is the timestep at which the Jacobian is to be computed.
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
    n, m, l, T = get_problem_dims(net)

    # Relevant sizes for the problem
    kdims = kkt_dims(n, m, l)
    sdims = storage_kkt_dims(n, l)

    # D_x ∇_x L 
    # dim = 3n, 3n
    K11 = spzeros(3n, 3n) # dim = n for s, ch, and dis 
    # D_x F^T
    # dim = 3n, (6n+2l)
    K12 = [
            -I(n) I(n) spzeros(n, 4n+2l);
            spzeros(n, 2n) -I(n) I(n) spzeros(n, 2n+2l);
            spzeros(n, 4n) -I(n) I(n) spzeros(n, 2l)
        ]
    # D_x H^T
    # dim = 3n, n
    K13 = [
        -I(n);
        η_c*I(n);
        -(1/η_d)*I(n)
    ]
    # Diag(λ) * D_xF
    # dim = 6n, 6n
    D_xF = K12'
    D = Diagonal([λsl; λsu; λchl; λchu; λdisl; λdisu; λrampl; λrampu])
    K21 = D * D_xF
    
    # diag(F)
    # dim = 6n, 6n
    K22 = Diagonal([-s; s-C; -ch; ch-P; -dis; dis-P; (gt_prev - ρ - gt); (gt - gt_prev - ρ)])
    
    # diagonal part of the Jacobian
    K_storage_t = [
        K11 K12 K13;
        K21 K22 spzeros(6n+2l, n);
        K13' spzeros(n, 7n+2l) 
    ]

    # static part
    # includes "cross terms" 
    # - with the ν dual variable
    # - with λgl, λgu variables
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
    K_static = [
        spzeros(n, T*kdims);
        spzeros(n, t*kdims - m - 1) -F' ones(n) spzeros(n, (T-t)*kdims);
        spzeros(n, t*kdims - m - 1) F' -ones(n) spzeros(n, (T-t)*kdims);
        spzeros(6n, T*kdims);
        Dλramp;
        spzeros(n, T*kdims);
    ]
    
    if t > 1
        K_storage_left = [
            spzeros(sdims-n, (t-1)*sdims);
            spzeros(n, (t-2)*sdims) I(n) spzeros(n, sdims-n);
        ]
    else
        K_storage_left = spzeros(sdims, 0);
    end
    if t < T
        K_storage_right = [
            spzeros(n, sdims-n) I(n) spzeros(n, (T-(t+1))*sdims);
            spzeros(sdims-n, (T-t)*sdims);
        ]
    else
        K_storage_right = spzeros(sdims, 0);
    end

    K_storage = [K_storage_left K_storage_t K_storage_right];

    return [K_static K_storage]
   
end
