###################################################################################
#
# Below contains older functions for automated testing of the sensitivity
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