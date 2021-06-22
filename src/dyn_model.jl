# Code for building and running a power management problem


# ===
# POWER MANAGEMENT PROBLEM
# ===

TAU = 1e-5

mutable struct PowerManagementProblem
    problem::Problem
    p
    g
    s
    params
end

"""
    PowerManagementProblem(f, d, pmax, gmax, A, ds; τ=1e-5)

Set up a power management problem with generation costs `f`, demand `d`,
maximum power flow `pmax`, maximum generation `gmax`, and network
incidence matrix `A`. `ds` is the change in battery state of charge. It should
be set to 0 if no battery is in the network. 

The parameter `τ` is a regularization weight used to make the problem
strongly convex by adding τ ∑ᵢ pᵢ² to the objective. Parameters fq and fl 
are the Quadratic and Linear terms in the generator costs. P is the maximum battery charge and 
discharge power (i.e. the maximum amount of energy that can be exchanged with a nodal battery
during a given timestep).

TODO: check that we're not messing up the ordering!
"""
function PowerManagementProblem(fq, fl, d, pmax, gmax, A, ds; τ=TAU)
    n, m = size(A)
    g = Variable(n)
    p = Variable(m)

    # add variable-specific constraints
    add_constraints(p, [
        -p <= pmax, 
        p <= pmax, 
    ])
    add_constraints!(g, [
        -g <= 0, 
        g <= gmax
    ]
    )

    # build problem
    problem = minimize(
        (1/2)*quadform(g, Diagonal(fq))
        + fl'g
        + (τ/2)*sumsquares(p)
    )
    # add network constraint
    add_constraint!(problem, 0 == A*p - g + d + ds)

    params = (fq=fq, fl=fl, d=d, pmax=pmax, gmax=gmax, A=A, τ=τ)

    # individual PMPs do not have s assigned
    return PowerManagementProblem(problem, p, g, nothing, params)
end

# ===
# Dynamic PowerManagementProblem
# ===
"""
Wrapper and PowerManagementProblem to build a sequence of these. 

The arguments to DynPowerManagementProblem are the same as for PowerManagementProblem
except they are arrays, with dimension equal to the time horizon (i.e. number of timesteps) T. 

P is the maximum charge/discharge power. 
C is the maximum SOC of the batteries.

TODO: parameterize the management of initial and terminal constraints
"""
function DynPowerManagementProblem(
    fq, fl, d, pmax, gmax, A, P, C; τ=TAU
)
    T = length(fq)
    n, _ = size(A)

    # Define a storage variable for n nodes over T timesteps
    # TODO: for now the edge constraints are neglected
    s = Variable(n, T)
    add_constraints!(s,[
        0 <= s <= C
        # s[:, 1] == 0, 
        # s[:, T+1] == 0
    ])
    add_constraints!(s, [
        s[t+1] - s[t] <= P for t in 1:T
    ])
    add_constraints!(s, [
        s[t] - s[t+1] <= P for t in 1:T
    ])

    subproblems = [
        # first treating the initial constraint explicitly to avoid having
        # specify another variable for s_0
        PowerManagementProblem(
            fq[t], fl[t], d[t], pmax[t], gmax[t], A, (s[1] - C)
            ), 
        PowerManagementProblem(
            fq[t], fl[t], d[t], pmax[t], gmax[t], A, (s[t] - s[t-1])
            ) for t in 2:T
        ]

    objective = sum([sub.problem.objective for sub in subproblems])
    g_T = [sub.g for sub in subproblems]
    p_T = [sub.p for sub in subproblems]

    dynProblem = minimize(
        objective
    )

    params = (fq=fq, fl=fl, d=d, pmax=pmax, gmax=gmax, A=A, P=P, C=C, τ=τ)

    return PowerManagementProblem(dynProblem, p_T, g_T, s, params)
end

# ===
# Utils
# ===

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
    P.problem.constraints[5].dual[:, 1]




# ===
# KKT OPERATOR for the PMP
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

# ===
# KKT OPERATOR for the DynPMP
# ===

"""
    kkt_dims(n, m)

Compute the dimensions of the input / output of the KKT operator for
a network with `n` nodes and `m` edges.
"""
kkt_dims_dyn(n, m) = 4n + 3m

"""
    kkt(x, fq, fl, d, pmax, gmax, A; τ=TAU)

Compute the KKT operator applied to `x`, with parameters given by `fq`,
`fl`, `d`, `pmax`, `gmax`, `A`, and `τ`.
"""
function kkt_dyn(x, fq, fl, d, pmax, gmax, A; τ=TAU)
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

function flatten_variables_dyn(P::PowerManagementProblem)
    x = [evaluate(P.g); evaluate(P.p)]
    λ = vcat([c.dual for c in P.problem.constraints]...)
    return [x; λ][:, 1]
end

function unflatten_variables_dyn(x, n, m)
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
