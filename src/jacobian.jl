
# =====
# IMPLICIT DIFFERENTIATION
# =====

function pullback_pmin(pmp, pmp_result, demand_devices, dz)
    (; devices) = pmp

    Jz = jacobian_kkt_z(pmp, pmp_result)

    K(θ, k) = kkt_local_params(pmp, pmp_result, k, θ, devices[k].pmax, devices[k].cost)

    # Note: Zygote actually computes the Jacobian transposed
    # How misleading!
    _jacs = [forward_jacobian(θ -> K(θ, k), devices[k].pmin)[2] for k in demand_devices]
    JθT = reduce(vcat, _jacs)

    dK = Jz' \ dz

    return -JθT * (Jz' \ dz)
end


# =====
# SYSTEM JACOBIAN
# =====

function jacobian_kkt_z(pmp::DynamicPowerManagementProblem, pmp_result) 
    dim = get_dims(pmp)
    J = spzeros(dim, dim)
    
    return jacobian_kkt_z!(J, pmp, pmp_result)
end

function jacobian_kkt_z!(J, pmp::DynamicPowerManagementProblem, pmp_result)
    (; devices, F, T) = pmp
    (; p, devc, devv) = pmp_result
    
    n = length(devices)
    m = size(F, 1)

    dim = get_dims(pmp)
    J = spzeros(dim, dim)
    
    # Global part
    _jacobian_kkt_z_network!(J, pmp, pmp_result)

    # Local parts
    J_devs = sparse.([_get_jac_kkt_z_device(pmp, p, devc, devv, i) for i in 1:n])
    J[(T+T*m+1):end, (T+T*m+1):end] = blockdiag(J_devs...)

    return J
end

function _get_jac_kkt_z_device(pmp, p, devc, devv, i)
    inds = get_device_inds(pmp, i)
    p = p[i, :]
    μ = [make_vec(c.dual) for c in devc[i]]

    return jacobian_kkt_z_device(pmp.devices[i], p, μ, devv[i])
end

function _jacobian_kkt_z_network!(J, pmp, pmp_result)
    (; devices, F, fmax, T) = pmp
    (; p, netc, devc, problem) = pmp_result

    n = length(devices)
    m = size(F, 1)
    
    ν, λ = reshape(netc[1].dual, :), netc[2].dual
    λt = [λ[:, t] for t in 1:T]
    
    # Energy balance constraint K_ν (first T rows)
    for i in 1:n
        inds = get_device_primary_var_inds(pmp, i)
        J[1:T, inds] .= I(T)
    end

    # Flow constraint K_λ
    # Part 1: dK_λ/dλt
    for t in 1:T
        block = (T+(t-1)*m+1):(T+t*m)
        J[block, block] = Diagonal(F*p[:, t] - fmax[:, t])
    end
    
    # Part 2: dK_λ/dpi
    for i in 1:n
        inds = get_device_primary_var_inds(pmp, i)
        fi = F[:, i]
        
        for t in 1:T
            block = (T+(t-1)*m+1):(T+t*m)
            @. J[block, inds[t]] = λt[t] * fi  # LITTLE SLOW, TODO: use blockdiag
        end
    end

    # Last but not least: dK_pi/dν and dK_pi/dλ
    for i in 1:n
        inds = get_device_primary_var_inds(pmp, i)
        fi = sparse(F[:, i])

        # dK_pi/dν
        J[inds, 1:T] = Diagonal(ν)

        # dK_pi/dλ
        J[inds, (T+1):(T+T*m)] = blockdiag([sparse(fi') for _ in inds]...)
    end
    
    return J
end




# =====
# PARAMETER JACOBIANS
# =====

function kkt_local_params(pmp::DynamicPowerManagementProblem, pmp_result, i, params...)
    d = typeof(pmp.devices[i])(params...)
    p = pmp_result.p[i, :]
    μ = [make_vec(c.dual) for c in pmp_result.devc[i]]

    dims = get_dims(pmp)
    inds = get_device_inds(pmp, i)

    n_before = minimum(inds) - 1
    n_after = dims - maximum(inds)
    
    return [
        zeros(n_before)
        get_device_dL_dx(d, p, μ);
        get_device_comp_slack(d, p, μ)
        zeros(n_after)
    ]
end