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
    F
    τ
end

"""
    PowerNetwork(fq, fl, pmax, gmax, A, B, F; τ=TAU)

Create a power network with quadratic and linear prices `fq` and `fl`, line capacities
`pmax`, generation capacities `gmax`, incidence matrix `A`, node-generator matrix `B`,
and PFDF matrix `F`.

Optionally, quadratically penalize power flows with weight `(1/2) τ^2`.
"""
PowerNetwork(fq, fl, pmax, gmax, A, B, F; τ=TAU) =
    PowerNetwork(fq, fl, pmax, gmax, A, B, F, τ)

mutable struct PowerManagementProblem
    problem::Problem
    p
    g
    s
    ch
    dis
    params
end


"""
    PowerManagementProblem(fq, fl, d, pmax, gmax, A, B, F; τ=TAU, ch=0, dis=0, S=0, η_c=1.0, η_d=1.0)

Set up a static power management problem with the following parameters: 
- `fq`: quadratic generator costs
- `fl`: linear generator costs
- `d`: nodal demand
- `pmax`: maximum power flow
- `gmax`: generation capacities
- `A`: network incidence matrix
- `B`: generator-to-node matrix
- `F`: PFDF matrix
- `τ`: regularization weight used to make the problem strongly convex by adding τ ∑ᵢ pᵢ² to the objective
- `ch`: additional charge over one timestep in the batteries of the network
- `dis`: additional discharge over onetimestep in the batteries of the network
- `S`: battery-to-node matrix
- `η_c` and `η_d`: battery charge and discharge efficiencies
"""
function PowerManagementProblem(fq, fl, d, pmax, gmax, A, B, F; τ=TAU, ch=0, dis=0, S=0, η_c=1.0, η_d=1.0)
    n, m = size(A)
    n, l = size(B)
    g = Variable(l)
    p = Variable(m)

    problem = minimize(
        (1/2)*sumsquares(dot(*)(sqrt.(fq), g))
        + fl'g
        + (τ/2)*sumsquares(p)
    )
    add_constraints!(problem, [
        -p <= pmax,  # λpl
        p <= pmax,  # λpu
        -g <= 0,  # λgl
        g <= gmax,  # λgu
        0 == p - F*(B*g - d - S*ch + S*dis),  # ν
        0 == ones(n)' * (B*g - d - S*ch + S*dis),  # νE
    ])

    params = (fq=fq, fl=fl, d=d, pmax=pmax, gmax=gmax, A=A, B=B, F=F, S=S, τ=τ, η_c=η_c, η_d=η_d)

    return PowerManagementProblem(problem, p, g, zeros(n), zeros(n), zeros(n), params)
end

PowerManagementProblem(net::PowerNetwork, d) =
    PowerManagementProblem(net.fq, net.fl, d, net.pmax, net.gmax, net.A, net.B, net.F; τ=net.τ)

"""
    Convex.solve!(P::PowerManagementProblem, opt; verbose=false)
"""
Convex.solve!(P::PowerManagementProblem, opt; verbose=false) = 
    solve!(P.problem, opt; silent_solver=!verbose, verbose=verbose)

"""
    get_lmps(P::PowerManagementProblem)

Return the locational marginal prices of `P` (assumes the problem has 
already been solved). The LMPs are the dual variables of the power
conservation constraint.
"""
get_lmps(P::PowerManagementProblem) = -P.problem.constraints[end-1].dual[:, 1] .- P.problem.constraints[end].dual[1]

# ===
# KKT OPERATOR
# ===

"""
    kkt_dims(n, m, l)

Compute the dimensions of the input / output of the KKT operator for
a network with `n` nodes and `m` edges and `l` generator per node

"""
kkt_dims(n, m, l) = 4m + 3l + 1

"""
    kkt(x, fq, fl, d, pmax, gmax, A; τ=TAU)

Compute the KKT operator applied to `x`, with: 
- `fq`: quadratic generator costs
- `fl`: linear generator costs
- `d`: nodal demand
- `pmax`: maximum power flow
- `gmax`: generation capacities
- `A`: network incidence matrix
- `τ`: regularization weight used to make the problem strongly convex by adding τ ∑ᵢ pᵢ² to the objective
"""
function kkt(x, fq, fl, d, pmax, gmax, A, B, F; τ=TAU, ch=0, dis=0,  S=0)
    n, m = size(A)
    n, l = size(B)

    g, p, λpl, λpu, λgl, λgu, ν, νE = unflatten_variables(x, n, m, l)

    # Lagragian is
    # L = J + λpl'(-p - pmax) + ... + λgu'(g - gmax) + v'(Ap - g - d) + νF'(p - F*(B*g - d))
    return [
        Diagonal(fq)*g + fl - λgl + λgu - B'*F'*ν + νE[1]*B'*ones(n); # stationarity: ∇_g L
        ν + λpu - λpl + τ*p; # stationarity: ∇_p L
        λpl .* (-p - pmax); # complementary slackness
        λpu .* (p - pmax); # complementary slackness
        -λgl .* g; # complementary slackness
        λgu .* (g - gmax); # complementary slackness
        p - F*(B*g - d .- S*ch .+ S*dis); # equality constraints
        ones(n)' * (B*g - d .- S*ch .+ S*dis); # equality constraints
    ]
end

kkt(x, net::PowerNetwork, d) =
    kkt(x, net.fq, net.fl, d, net.pmax, net.gmax, net.A, net.B, net.F; τ=net.τ)


"""
    flatten_variables(P::PowerManagementProblem)

Concatenates primal and dual variables from `P` into a single vector.
"""
function flatten_variables(P::PowerManagementProblem)
    x = [evaluate(P.g); evaluate(P.p)]
    λ = vcat([c.dual for c in P.problem.constraints]...)
    return [x; λ][:, 1]
end


"""
    unflatten_variables(x, n, m, l)

Extracts primal and dual variables from `x`.
"""
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
    
    ν = x[i+1:i+m]
    i += m

    νE = x[i+1:i+1]
    i += 1

    return g, p, λpl, λpu, λgl, λgu, ν, νE
end
