# Code for computing sensitivities with respect to various parameters

"""
    get_problem_dims(net::PowerNetwork)

Return `(n, m, l)`, where `n` is the number of nodes in the network,
`m` is the number of edges, and `l` is the number of generators.
"""
function get_problem_dims(net::PowerNetwork)
    n, m = size(net.A)
    n, l = size(net.B)

    return n, m, l
end

"""
    get_problem_dims(net::DynamicPowerNetwork)

Return `(n, m, l, T)`, where `n` is the number of nodes in the network,
`m` is the number of edges, and `l` is the number of generators,
and `T` is the number of timesteps in the dynamic problem.
"""
function get_problem_dims(net::DynamicPowerNetwork)
    n, m = size(net.A)
    _, l = size(net.B)
    _, ns = size(net.S)

    return n, m, l, ns, net.T
end


# ===
# SENSITIVITIES
# ===

"""
    sensitivity_demand(P::PowerManagementProblem, ∇C, fq, fl, d, pmax, gmax, A, B, F)

Compute `∇_d C( x_opt(d) )`, where `x_opt(d)` is the optimal solution (primal
and dual variables) of the power management problem `P` with parameters 
`(fq, fl, d, pmax, gmax, A, B, F)` and `∇C` is the gradient `∇_x C(x)`.
"""
function sensitivity_demand(P::PowerManagementProblem, ∇C, fq, fl, d, pmax, gmax, A, B, F)
    x = flatten_variables(P)

    # Get partial Jacobians of KKT operator
    ∂K_xT = adjoint(compute_jacobian_kkt(fq, fl, d, pmax, gmax, A, B, F, x))
    _, ∂K_θT = Zygote.forward_jacobian(d -> kkt(x, fq, fl, d, pmax, gmax, A, B, F), d)

    # Now compute ∇C(g*(θ)) = -∂K_θ' * inv(∂K_x') * v
    v = ∂K_xT \ ∇C
    ∇C_θ = -∂K_θT * v

    return ∇C_θ
end

sensitivity_demand(P::PowerManagementProblem, ∇C, net::PowerNetwork, d) = 
    sensitivity_demand(
        P::PowerManagementProblem, ∇C, net.fq, net.fl, d, net.pmax, 
        net.gmax, net.A, net.B, net.F
    )

"""
    compute_mefs(P::PowerManagementProblem, net::PowerNetwork, d, c)

Compute the marginal emission factors given carbon costs `c` and demand `d`.
"""
function compute_mefs(P::PowerManagementProblem, net::PowerNetwork, d, c)
    # Extract dimensions
    n, m, l = get_problem_dims(net)

    # Construct ∇_x C(x)
    ∇C = zeros(kkt_dims(n, m, l), size(c, 2))
	∇C[1:l, :] .= c

    # Return sensitivity
	return sensitivity_demand(P, ∇C, net, d)
end

"""
    sensitivity_price(P::PowerManagementProblem, ∇C, fq, fl, d, pmax, gmax, A)

Returns `∇C_θq` and `∇C_θl`, the sensitivity of `∇C` with respect to
the quadratic generation cost coefficients `fq` and the linear cost
coefficients `fl`, respectively. 
"""
function sensitivity_price(P::PowerManagementProblem, ∇C, fq, fl, d, pmax, gmax, A, B, F)
    x = flatten_variables(P)

    # Get partial Jacobians of KKT operator
    # _, ∂K_xT = Zygote.forward_jacobian(x -> kkt(x, f, d, params), x)
    ∂K_x = compute_jacobian_kkt(fq, fl, d, pmax, gmax, A, B, F, x)
    ∂K_xT = ∂K_x'

    # _, ∂K_θT = Zygote.forward_jacobian(f -> kkt(x, f, d, params), f)
    ∂K_θq = compute_jacobian_quad_price(A, B, x)
    ∂K_θqT = ∂K_θq'

    ∂K_θl = compute_jacobian_lin_price(A, B)
    ∂K_θlT = ∂K_θl'

    # Now compute ∇C(g*(θ)) = -∂K_θ' * (∂K_x' * v)
    v = ∂K_xT \ ∇C
    ∇C_θq = -∂K_θqT * v
    ∇C_θl = -∂K_θlT * v

    return ∇C_θq, ∇C_θl
end

sensitivity_price(P::PowerManagementProblem, ∇C, net::PowerNetwork, d) = 
    sensitivity_price(
        P::PowerManagementProblem, ∇C, net.fq, net.fl, 
        d, net.pmax, net.gmax, net.A, net.B, net.F
    )


# ===
# JACOBIANS
# ===

"""
    compute_jacobian_quad_price(A, x)

Compute the Jacobian of the KKT operator with respect to the quadratic
generation cost coefficients.
"""
function compute_jacobian_quad_price(A, B, x)
    n, m = size(A)
    n, l = size(B)
    g, _, _, _, _, _, _, _ = unflatten_variables(x, n, m, l)

    return [diagm(g); spzeros(3m + 2l + n, l)]
end

"""
    compute_jacobian_lin_price(A)

Compute the Jacobian of the KKT operator with respect to the linear
generation cost coefficients.
"""
function compute_jacobian_lin_price(A, B)
    n, m = size(A)
    n, l = size(B)

    return [I(l); spzeros(3m + 2l + n, l)]
end

"""
    compute_jacobian_kkt(fq, fl, d, pmax, gmax, A, x; τ=TAU)

Compute the Jacobian of the KKT operator.

Notes 
-----
The Jacobian below is computed by decomposing it in several blocks, following Barrat 2018. 
"""
function compute_jacobian_kkt(fq, fl, d, pmax, gmax, A, B, F, x; τ=TAU)
    n, m = size(A)
    n, l = size(B)

    g, p, λpl, λpu, λgl, λgu, _, _ = unflatten_variables(x, n, m, l)

    # Jacobian of stationarity conditions wrt primal variables
    # Here primal variables include (p, g) only
    K11 = [
        Diagonal(fq)     spzeros(l, m);
        spzeros(m, l)    τ * I(m)
    ]

    # Jacobian of stationarity conditions wrt dual variables of inequality constraints
    K12 = [
        spzeros(l, 2 * m)   -I(l)       I(l);
        -I(m)             I(m)        spzeros(m, 2l)
    ]

    # Jacobian of complementary slackness wrt primal variables
    K21 = [
        spzeros(m, l) -Diagonal(λpl);
        spzeros(m, l) Diagonal(λpu);
        -Diagonal(λgl) spzeros(l, m);
        Diagonal(λgu) spzeros(l, m)
    ]
    
    # Jacobian of complementary slackness wrt dual variables of inequality constraints
    K22 = Diagonal([-p - pmax; p - pmax; -g; g - gmax])
    
    # Jacobian of stationarity wrt dual variables of equality constraints
    K31 = [
        -F*B I(m);
        ones(n)'B spzeros(1, m);
    ]

    return [
        K11 K12 K31';
        K21 K22 spzeros(2 * (m + l), m + 1);
        K31 spzeros(m + 1, 2 * (m + l) + m + 1)
    ]
end