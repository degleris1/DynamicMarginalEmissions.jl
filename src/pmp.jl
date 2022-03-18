mutable struct DynamicPowerManagementProblem
    devices::Vector{<:AbstractDevice}
    F
    fmax
    T::Int
end




# =====
# SOLVER
# =====

function solve!(pmp::DynamicPowerManagementProblem, solver=ECOS.Optimizer)
    (; devices, F, fmax, T) = pmp

    @assert all([get_time_horizon(d) == T for d in devices])
    
    n = length(devices)
    p = Variable(n, T)

    # Instantiate problem
    costs = [get_cost(devices[i], p[i, :]') for i in 1:n]
    problem = minimize(sum(costs))

    # Add global network constraints
    network_constraints = [0 == ones(n)' * p, F*p <= fmax]
    add_constraints!(problem, network_constraints)

    # Add local variables
    device_vars = make_aux_vars.(devices)

    # Add local constraints
    device_constraints = [make_constraints(devices[i], p[i, :]', device_vars[i]) for i in 1:n]
    for c in device_constraints
        add_constraints!(problem, c)
    end

    # Solve
    solve!(problem, solver)

    p = evaluate(p)
    device_vars = map(evaluate, device_vars)

    return (p=p, netc=network_constraints, devc=device_constraints, devv=device_vars, problem=problem)
end




# =====
# KKT CONDITIONS
# =====

function kkt(pmp, pmp_result)
    (; devices, F, fmax, T) = pmp
    (; p, netc, devc, devv, problem) = pmp_result

    n = length(devices)
    p = make_vec(p)
    ν, λ = reshape(netc[1].dual, :), netc[2].dual

    # Energy balance constraint
    resid = make_vec(p' * ones(n))

    # Transmission constraints
    K_λ = λ .* (F*p - fmax)
    for t in 1:T
        append!(resid, K_λ[:, t])
    end

    # Device objectives and constraints
    for i in 1:n
        μi = [make_vec(c.dual) for c in devc[i]]
        pi = p[i, :]
        fi = F[:, i]

        K_pi = get_device_dL_dx(devices[i], pi, μi, devv[i])
        K_pi += ν + λ'*fi
        append!(resid, K_pi)

        K_ui = get_device_dL_daux(devices[i], pi, μi, devv[i])
        append!(resid, K_ui)

        K_μi = get_device_comp_slack(devices[i], pi, μi, devv[i])
        append!(resid, K_μi)
    end

    return resid
end




# =====
# DIMENSIONS - indexing, flattening, etc
# =====

get_time_horizon(pmp::DynamicPowerManagementProblem) = pmp.T

get_num_lines(pmp::DynamicPowerManagementProblem) = size(pmp.F, 1)

get_num_devices(pmp::DynamicPowerManagementProblem) = length(pmp.devices)

function get_dims(pmp::DynamicPowerManagementProblem)
    T, m = get_time_horizon(pmp), get_num_lines(pmp)
    return T + m*T + sum([get_dims(d) for d in pmp.devices])
end

function get_device_inds(pmp, i)
    T = pmp.T
    m = size(pmp.F, 1)
    
    offset = T + m*T + sum(get_dims.(pmp.devices[1:i-1]))
    dim = get_dims(pmp.devices[i])
    
    return (offset+1):(offset+dim)
end

function get_device_primary_var_inds(pmp, i)
    T, m = get_time_horizon(pmp), get_num_lines(pmp)
    offset = T + m*T + sum(Int[get_dims(d) for d in pmp.devices[1:i-1]])
    return (offset+1):(offset+T)
end

function get_device_aux_var_inds(pmp, i)
    T, m = get_time_horizon(pmp), get_num_lines(pmp)
    offset = T + m*T + sum([get_dims(d) for d in pmp.devices[1:i-1]])

    num_aux = get_num_aux(pmp.devices[i])
    return (offset+T+1):(offset+T+num_aux)
end

function get_device_constraint_inds(pmp, i)
    T = pmp.T
    m = size(pmp.F, 1)
    
    offset = T + m*T + sum([get_dims(d) for d in pmp.devices[1:i-1]])
    dim = get_dims(pmp.devices[i])
    return (offset+T+1):(offset+dim)
end

function flatten_pmp_result(pmp, pmp_result)
    (; devices, F, fmax, T) = pmp
    (; p, netc, devc, devv, problem) = pmp_result

    n = length(devices)
    ν, λ = reshape(netc[1].dual, :), netc[2].dual

    # Energy balance constraint
    z = make_vec(ν)

    # Transmission constraint
    for t in 1:T
        append!(z, λ[:, t])
    end

    # Device variables and constraints
    for i in 1:n
        append!(z, p[i, :])
        for v in something(devv[i], [])
            append!(z, v)
        end
        for c in devc[i]
            append!(z, c.dual)
        end
    end

    return z
end

function unflatten_pmp_result(pmp, z)
    n, m, T = get_num_devices(pmp), get_num_lines(pmp), get_time_horizon(pmp)

    net_duals = (z[1:T], reshape(z[(T+1):(T+T*m)], m, T))

    p = [z[get_device_primary_var_inds(pmp, i)] for i in 1:n]
    p = adjoint(reduce(hcat, p))

    dev_vars = [z[get_device_aux_var_inds(pmp, i)] for i in 1:n]
    dev_duals = [z[get_device_constraint_inds(pmp, i)] for i in 1:n]

    return (p=p, net_duals=net_duals, dev_duals=dev_duals, dev_vars=dev_vars)
end




# =====
# HELPERS
# =====
make_vec(x::Array) = x
make_vec(x) = [x]

reshape(x::Real, ::Colon) = [x]

evaluate(x::Nothing) = nothing
evaluate(x::Tuple) = map(evaluate, x)
evaluate(x::NamedTuple) = map(evaluate, x)