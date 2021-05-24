
mutable struct PowerManagementProblem
    problem::Problem
    p
    g
end

"""
    PowerManagementProblem(f, d, pmax, gmax, A; τ=1e-4)

Set up a power management problem with generation costs `f`, demand `d`,
maximum power flow `pmax`, maximum generation `gmax`, and network
incidence matrix `A`.

The parameter `τ` is a regularization weight used to make the problem
strongly convex by adding τ ∑ᵢ pᵢ² to the objective.
"""
function PowerManagementProblem(f, d, pmax, gmax, A; τ=1e-5)
    n, m = size(A)
    g = Variable(n)
    p = Variable(m)

    problem = minimize(f'g + τ*sumsquares(p))
    add_constraints!(problem, [
        -p <= pmax,
        p <= pmax,
        -g <= 0, 
        g <= gmax,
        A*p - g + d == 0,
    ])

    return PowerManagementProblem(problem, p, g)
end

Convex.solve!(prob::PowerManagementProblem, opt; verbose=false) = 
    solve!(prob.problem, opt; verbose=verbose)

get_lmps(P::PowerManagementProblem) = -P.problem.constraints[5].dual[:, 1]

function sensitivity_demand(P::PowerManagementProblem, ∇C, params)
    x = flatten_variables(P)
    f, d = params[1:2]

    # Get partial Jacobians of KKT operator
    _, ∂K_xT = Zygote.forward_jacobian(x -> kkt(x, f, d, params), x)
    _, ∂K_θT = Zygote.forward_jacobian(d -> kkt(x, f, d, params), d)

    # Now compute ∇C(g*(θ)) = -∂K_θ' * (∂K_x' * v)
    v = ∂K_xT \ ∇C
    ∇C_θ = -∂K_θT * v

    return ∇C_θ
end

function sensitivity_price(P::PowerManagementProblem, ∇C, params)
    x = flatten_variables(P)
    f, d = params[1:2]

    # Get partial Jacobians of KKT operator
    _, ∂K_xT = Zygote.forward_jacobian(x -> kkt(x, f, d, params), x)
    _, ∂K_θT = Zygote.forward_jacobian(f -> kkt(x, f, d, params), f)

    # Now compute ∇C(g*(θ)) = -∂K_θ' * (∂K_x' * v)
    v = ∂K_xT \ ∇C
    ∇C_θ = -∂K_θT * v

    return ∇C_θ
end


# ===
# DIFFERENTIATION UTILS
# ===

function kkt(x, f, d, params; τ=1e-5)
    (_, _, pmax, gmax, A) = params
    n, m = size(A)

    g, p, λpl, λpu, λgl, λgu, ν = unflatten_variables(x, n, m)

    return [
        f + ν - λgl + λgu
        -A'ν + λpu - λpl + 2τ*p;
        λpu .* (p - pmax);
        -λpl .* (-p - pmax);
        -λgl .* g;
        λgu .* (g - gmax);
        A*p - g + d;
    ]
end

function kkt_dims(n, m)
    return 4n + 3m
end

function flatten_variables(P::PowerManagementProblem)
    x = [evaluate(P.g); evaluate(P.p)]
    λ = vcat([c.dual for c in P.problem.constraints]...)
    return [x; λ][:, 1]
end

function unflatten_variables(x, n, m)
    i = 0
    
    g = x[i+1:i+n]
    i += n
    
    p = x[i+1:i+m]
    i += m
    
    λpl = x[i+1:i+m]
    i += m
    λpu = x[i+1:i+m]
    i += m
    
    λgl = x[i+1:i+n]
    i += n
    λgu = x[i+1:i+n]
    i += n
    
    ν = x[i+1:i+n]
    
    return g, p, λpl, λpu, λgl, λgu, ν
end

function loss_and_grad(f̂, B, case_list, pmax, gmax, A)
    n, m = size(A)
    L = 0.0
    df = zeros(n)
    
    T = length(case_list)
    _∇L = zeros(kkt_dims(n, m))

    for case in case_list
        params = (f̂, case.d, pmax, gmax, A)
        opf = PowerManagementProblem(params...)
        solve!(opf, () -> ECOS.Optimizer(verbose=false), verbose=false)
        ĝ = evaluate(opf.g)
        
        L += (1/2) * norm(B*ĝ - case.g)^2 / T
        
        _∇L[1:n] .+= B' * (B*ĝ - case.g)
        df += sensitivity_price(opf, _∇L, params) / T
    end
    
    return L, df
end

function stochastic_loss_and_grad(f̂, B, case_list, pmax, gmax, A, sample)
    return loss_and_grad(f̂, B, case_list[sample], pmax, gmax, A)
end