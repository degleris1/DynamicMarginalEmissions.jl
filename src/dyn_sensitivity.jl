# Computing sensitivities in the dynamic PMP setting

function sensitivity_demand_dyn(P::PowerManagementProblem, net::DynamicPowerNetwork, d, ∇C)

    x = flatten_variables_dyn(P)

    # Get partial Jacobians of KKT operator
    _, ∂K_xT = Zygote.forward_jacobian(x -> kkt_dyn(x, net, d), x)
    _, ∂K_θT = Zygote.forward_jacobian(d -> kkt_dyn(x, net, d), d)

    # Now compute ∇C(g*(θ)) = -∂K_θ' * (∂K_x' * v)
    v = ∂K_xT \ ∇C
    ∇C_θ = -∂K_θT * v

    return ∇C_θ
end