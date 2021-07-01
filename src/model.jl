# Code for building and running a power management problem


# ===
# POWER MANAGEMENT PROBLEM
# ===

TAU = 0.0

mutable struct PowerNetwork
    fq
    fl
    pmax
    gmax
    A
    B
    τ
end

PowerNetwork(fq, fl, pmax, gmax, A, B; τ=TAU) =
    PowerNetwork(fq, fl, pmax, gmax, A, B, τ)

mutable struct PowerManagementProblem
    problem::Problem
    p
    g
    s
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
function PowerManagementProblem(fq, fl, d, pmax, gmax, A, B; τ=TAU, ds=0)
    n, m = size(A)
    n, l = size(B)
    g = Variable(l)
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
        0 == A*p - B*g + d + ds,
    ])

    params = (fq=fq, fl=fl, d=d, pmax=pmax, gmax=gmax, A=A, B=B, τ=τ)

    return PowerManagementProblem(problem, p, g, zeros(n), params)
end

PowerManagementProblem(net::PowerNetwork, d) =
    PowerManagementProblem(net.fq, net.fl, d, net.pmax, net.gmax, net.A, net.B; τ=net.τ)

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
    P.problem.constraints[end].dual[:, 1]




# ===
# KKT OPERATOR
# ===

"""
    kkt_dims(n, m)

Compute the dimensions of the input / output of the KKT operator for
a network with `n` nodes and `m` edges.

TODO: add `l` in docs
"""
kkt_dims(n, m, l) = 3m + 3l + n

"""
    kkt(x, fq, fl, d, pmax, gmax, A; τ=TAU)

Compute the KKT operator applied to `x`, with parameters given by `fq`,
`fl`, `d`, `pmax`, `gmax`, `A`, and `τ`.
"""
function kkt(x, fq, fl, d, pmax, gmax, A, B; τ=TAU, ds=0)
    n, m = size(A)
    n, l = size(B)

    g, p, λpl, λpu, λgl, λgu, ν = unflatten_variables(x, n, m, l)

    # Lagragian is
    # L = J + λpl'(-p - pmax) + ... + λgu'(g - gmax) + v'(Ap - g - d)
    return [
        Diagonal(fq)*g +  fl - B'ν - λgl + λgu; 
        A'ν + λpu - λpl + τ*p;
        λpl .* (-p - pmax);
        λpu .* (p - pmax);
        -λgl .* g;
        λgu .* (g - gmax);
        A*p - B*g + d .+ ds;
    ]
end

kkt(x, net::PowerNetwork, d) =
    kkt(x, net.fq, net.fl, d, net.pmax, net.gmax, net.A, net.B; τ=net.τ, ds=zeros(size(net.A)[1]))

function flatten_variables(P::PowerManagementProblem)
    x = [evaluate(P.g); evaluate(P.p)]
    λ = vcat([c.dual for c in P.problem.constraints]...)
    return [x; λ][:, 1]
end

function unflatten_variables(x, n, m, l)
    i = 0
    
    g =  x[i+1:i+l]
    i += l
    p = x[i+1:i+m]
    i += m
    
    λpl = x[i+1:i+m]
    i += m
    λpu = x[i+1:i+m]
    i += m
    
    λgl = x[i+1:i+l]
    i += l
    λgu = x[i+1:i+l]
    i += l
    
    ν = x[i+1:i+n]
    return g, p, λpl, λpu, λgl, λgu, ν
end
