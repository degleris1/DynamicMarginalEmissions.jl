# Code for computing sensitivities with respect to various parameters




# ===
# SENSITIVITIES
# ===

"""
    sensitivity_demand(P::PowerManagementProblem, ∇C, fq, fl, d, pmax, gmax, A, B)

Compute `∇_d C( x_opt(d) )`, where `x_opt(d)` is the optimal solution (primal
and dual variables) of the power management problem `P` with parameters 
`(fq, fl, d, pmax, gmax, A, B)` and `∇C` is the gradient `∇_x C(x)`.
"""
function sensitivity_demand(P::PowerManagementProblem, ∇C, fq, fl, d, pmax, gmax, A, B)
    x = flatten_variables(P)

    # Get partial Jacobians of KKT operator
    _, ∂K_xT = Zygote.forward_jacobian(x -> kkt(x, fq, fl, d, pmax, gmax, A, B), x)
    _, ∂K_θT = Zygote.forward_jacobian(d -> kkt(x, fq, fl, d, pmax, gmax, A, B), d)

    # Now compute ∇C(g*(θ)) = -∂K_θ' * (∂K_x' * v)
    v = ∂K_xT \ ∇C
    ∇C_θ = -∂K_θT * v

    return ∇C_θ
end

sensitivity_demand(P::PowerManagementProblem, ∇C, net::PowerNetwork, d) = 
    sensitivity_demand(P::PowerManagementProblem, ∇C, net.fq, net.fl, d, net.pmax, net.gmax, net.A, net.B)

"""
    sensitivity_price(P::PowerManagementProblem, ∇C, fq, fl, d, pmax, gmax, A)

Returns `∇C_θq` and `∇C_θl`, the sensitivity of `∇C` with respect to
the quadratic generation cost coefficients `fq` and the linear cost
coefficients `fl`, respectively. 
"""
function sensitivity_price(P::PowerManagementProblem, ∇C, fq, fl, d, pmax, gmax, A, B)
    x = flatten_variables(P)

    # Get partial Jacobians of KKT operator
    # _, ∂K_xT = Zygote.forward_jacobian(x -> kkt(x, f, d, params), x)
    ∂K_x = compute_jacobian_kkt(fq, fl, d, pmax, gmax, A, B, x)
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
    sensitivity_price(P::PowerManagementProblem, ∇C, net.fq, net.fl, d, net.pmax, net.gmax, net.A, net.B)



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
    g, _, _, _, _, _, _ = unflatten_variables(x, n, m, l)

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
"""
function compute_jacobian_kkt(fq, fl, d, pmax, gmax, A, B, x; τ=TAU)
    n, m = size(A)
    n, l = size(B)

    g, p, λpl, λpu, λgl, λgu, _ = unflatten_variables(x, n, m, l)

    K11 = [
        Diagonal(fq)     spzeros(l, m);
        spzeros(m, l)    τ * I(m)
    ]
   
    K12 = [
        spzeros(l, 2*m)   -I(l)       I(l);
        -I(m)             I(m)        spzeros(m, 2l)
    ]

    K13 = [-B'; A']

    K21 = [
        spzeros(m, l) -Diagonal(λpl);
        spzeros(m, l) Diagonal(λpu);
        Diagonal(λgu) spzeros(l, m);
        -Diagonal(λgl) spzeros(l, m)
    ]

    K22 = Diagonal([-p-pmax; p-pmax; -g; g-gmax])
    
    return [
        K11 K12 K13;
        K21 K22 spzeros(2*(m+l), n);
        K13' spzeros(n, 2*(m+l)+n)
    ]
end