# Computing sensitivities for dynamic power management problems

"""
    sensitivity_demand_dyn(P::PowerManagementProblem, net::DynamicPowerNetwork, d, ∇C, t)

Compute `∇_{d_t} C( x_opt(d_t) )`, where `x_opt(d_t)` is the optimal solution (primal
and dual variables) of the power management problem `P` with parameters 
`(fq, fl, d, pmax, gmax, A, B)` and demand `d_t` at time `t`, and where
`∇C` is the gradient `∇_x C(x)`.
"""
function sensitivity_demand_dyn(P::PowerManagementProblem, net::DynamicPowerNetwork, d, ∇C, t)
    T = length(d)
    x = flatten_variables_dyn(P)

    # Get partial Jacobians of KKT operator
    _, ∂K_xT = Zygote.forward_jacobian(x -> kkt_dyn(x, net, d), x)
    _, ∂K_θT = Zygote.forward_jacobian(
        dt -> kkt_dyn(x, net, [tp == t ? dt : d[tp] for tp in 1:T]), 
        d[t]
    )

    # Now compute ∇C(g*(θ)) = -∂K_θ' * inv(∂K_x') * v
    v = ∂K_xT \ ∇C
    ∇C_θ = -∂K_θT * v

    return ∇C_θ
end