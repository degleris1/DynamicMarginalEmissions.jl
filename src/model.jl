# Code for building and running a power management problem


# ===
# POWER MANAGEMENT PROBLEM
# ===

TAU = 1e-5

mutable struct PowerManagementProblem
    problem::Problem
    p
    g
    params
end

"""
    PowerManagementProblem(f, d, pmax, gmax, A; τ=1e-5)

Set up a power management problem with generation costs `f`, demand `d`,
maximum power flow `pmax`, maximum generation `gmax`, and network
incidence matrix `A`.

The parameter `τ` is a regularization weight used to make the problem
strongly convex by adding τ ∑ᵢ pᵢ² to the objective.
"""
function PowerManagementProblem(fq, fl, d, pmax, gmax, A; τ=TAU)
    n, m = size(A)
    g = Variable(n)
    p = Variable(m)

    problem = minimize(
        (1/2)*quadform(g, Diagonal(fq))
        + fl'g
        + (τ/2)*sumsquares(p)
    )
    add_constraints!(problem, [
        -p <= pmax,
        p <= pmax,
        -g <= 0, 
        g <= gmax,
        0 == A*p - g + d,
    ])

    params = (fq=fq, fl=fl, d=d, pmax=pmax, gmax=gmax, A=A, τ=τ)

    return PowerManagementProblem(problem, p, g, params)
end

"""
    Convex.solve!(P::PowerManagementProblem, opt; verbose=false)
"""
Convex.solve!(P::PowerManagementProblem, opt; verbose=false) = 
    solve!(P.problem, opt; verbose=verbose)

"""
    get_lmps(P::PowerManagementProblem)

Return the locational marginal prices of `P` (assumes the problem has 
already been solved). The LMPs are the dual variables of the power
conservation constraint.
"""
get_lmps(P::PowerManagementProblem) = 
    -P.problem.constraints[5].dual[:, 1]




# ===
# KKT OPERATOR
# ===

"""
    kkt_dims(n, m)

Compute the dimensions of the input / output of the KKT operator for
a network with `n` nodes and `m` edges.
"""
kkt_dims(n, m) = 4n + 3m

"""
    kkt(x, fq, fl, d, pmax, gmax, A; τ=TAU)

Compute the KKT operator applied to `x`, with parameters given by `fq`,
`fl`, `d`, `pmax`, `gmax`, `A`, and `τ`.
"""
function kkt(x, fq, fl, d, pmax, gmax, A; τ=TAU)
    n, m = size(A)

    g, p, λpl, λpu, λgl, λgu, ν = unflatten_variables(x, n, m)

    # Lagragian is
    # L = J + λpl'(-p - pmax) + ... + λgu'(g - gmax) + v'(Ap - g - d)
    return [
        Diagonal(fq)*g +  fl - ν - λgl + λgu; 
        A'ν + λpu - λpl + τ*p;
        λpl .* (-p - pmax);
        λpu .* (p - pmax);
        -λgl .* g;
        λgu .* (g - gmax);
        A*p - g + d;
    ]
end

function flatten_variables(P::PowerManagementProblem)
    x = [evaluate(P.g); evaluate(P.p)]
    λ = vcat([c.dual for c in P.problem.constraints]...)
    return [x; λ][:, 1]
end

function unflatten_variables(x, n, m)
    i = 0
    
    g =  x[i+1:i+n]
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
