# Code for computing sensitivities with respect to various parameters




# ===
# SENSITIVITIES
# ===

function sensitivity_demand(P::PowerManagementProblem, ∇C, fq, fl, d, pmax, gmax, A)
    x = flatten_variables(P)

    # Get partial Jacobians of KKT operator
    _, ∂K_xT = Zygote.forward_jacobian(x -> kkt(x, fq, fl, d, pmax, gmax, A), x)
    _, ∂K_θT = Zygote.forward_jacobian(d -> kkt(x, fq, fl, d, pmax, gmax, A), d)

    # Now compute ∇C(g*(θ)) = -∂K_θ' * (∂K_x' * v)
    v = ∂K_xT \ ∇C
    ∇C_θ = -∂K_θT * v

    return ∇C_θ
end

"""
    sensitivity_price(P::PowerManagementProblem, ∇C, fq, fl, d, pmax, gmax, A)

Returns `∇C_θq` and `∇C_θl`, the sensitivity of `∇C` with respect to
the quadratic generation cost coefficients `fq` and the linear cost
coefficients `fl`, respectively. 
"""
function sensitivity_price(P::PowerManagementProblem, ∇C, fq, fl, d, pmax, gmax, A)
    x = flatten_variables(P)

    # Get partial Jacobians of KKT operator
    # _, ∂K_xT = Zygote.forward_jacobian(x -> kkt(x, f, d, params), x)
    ∂K_x = compute_jacobian_kkt(fq, fl, d, pmax, gmax, A, x)
    ∂K_xT = ∂K_x'

    # _, ∂K_θT = Zygote.forward_jacobian(f -> kkt(x, f, d, params), f)
    ∂K_θq = compute_jacobian_quad_price(A, x)
    ∂K_θqT = ∂K_θq'

    ∂K_θl = compute_jacobian_lin_price(A)
    ∂K_θlT = ∂K_θl'

    # Now compute ∇C(g*(θ)) = -∂K_θ' * (∂K_x' * v)
    v = ∂K_xT \ ∇C
    ∇C_θq = -∂K_θqT * v
    ∇C_θl = -∂K_θlT * v

    return ∇C_θq, ∇C_θl
end




# ===
# JACOBIANS
# ===

"""
    compute_jacobian_quad_price(A, x)

Compute the Jacobian of the KKT operator with respect to the quadratic
generation cost coefficients.
"""
function compute_jacobian_quad_price(A, x)
    n, m = size(A)
    g, _, _, _, _, _, _ = unflatten_variables(x, n, m)

    return [diagm(g); spzeros(3m + 3n, n)]
end

"""
    compute_jacobian_lin_price(A)

Compute the Jacobian of the KKT operator with respect to the linear
generation cost coefficients.
"""
function compute_jacobian_lin_price(A)
    n, m = size(A)

    return [I(n); spzeros(3m + 3n, n)]
end

"""
    compute_jacobian_kkt(fq, fl, d, pmax, gmax, A, x; τ=TAU)

Compute the Jacobian of the KKT operator.
"""
function compute_jacobian_kkt(fq, fl, d, pmax, gmax, A, x; τ=TAU)
    n, m = size(A)

    g, p, λpl, λpu, λgl, λgu, _ = unflatten_variables(x, n, m)

    K11 = [
        Diagonal(fq)     spzeros(n, m);
        spzeros(m, n)    τ * I(m)
    ]
   
    K12 = [
        spzeros(n, 2*m)   -I(n)       I(n);
        -I(m)             I(m)        spzeros(m, 2n)
    ]

    K13 = [-I(n); A']

    K21 = [
        spzeros(m, n) -Diagonal(λpl);
        spzeros(m, n) Diagonal(λpu);
        Diagonal(λgu) spzeros(n, m);
        -Diagonal(λgl) spzeros(n, m)
    ]

    K22 = Diagonal([-p-pmax; p-pmax; -g; g-gmax])
    
    return [
        K11 K12 K13;
        K21 K22 spzeros(2*(m+n), n);
        K13' spzeros(n, 2*(m+n)+n)
    ]
end