# Code for building and running a dynamic power management problem


# ===
# POWER MANAGEMENT PROBLEM
# ===

INIT_COND = 0

"""
Here elements are assumed to be arrays of length equal to the time
horizon of the problem (T)

Parameters:
-----------
fq: quadratic costs
fl: linear costs
pmax: maximum power through edges
gmax: maximum generation 
A: incidence matrix
B: generator matrix
P: maximum charge/discharge power for storage
C: maximum SOC
τ: quadratic power penalty
"""
mutable struct DynamicPowerNetwork
    fq
    fl
    pmax
    gmax
    A
    B
    P
    C
    τ
end
DynamicPowerNetwork(fq, fl, pmax, gmax, A, B, P, C; τ=TAU) = 
    DynamicPowerNetwork(fq, fl, pmax, gmax, A, B, P, C, τ)

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
function DynamicPowerManagementProblem(
    fq, fl, d, pmax, gmax, A, B, P, C; τ=TAU
)
    T = length(fq)
    n, _ = size(A)

    # Define a storage variable for n nodes over T timesteps
    # TODO: for now the edge constraints are neglected
    s = [Variable(n) for _ in 1:T]

    subproblems = vcat(
        # first treating the initial constraint explicitly to avoid having to
        # specify another variable for s_0
        [PowerManagementProblem(
            fq[1], fl[1], d[1], pmax[1], gmax[1], A, B; ds=s[1] - INIT_COND
        )]
        ,
        # then iterating over all the timesteps 
        [PowerManagementProblem(
            fq[t], fl[t], d[t], pmax[t], gmax[t], A, B; ds=s[t] - s[t-1]
        ) for t in 2:T]
    )

    objective = sum([sub.problem.objective for sub in subproblems])
    g_T = [sub.g for sub in subproblems]
    p_T = [sub.p for sub in subproblems]

    dynProblem = minimize(objective)
    # adding the constraints of individual subproblems
    for sub in subproblems
        add_constraints!(dynProblem, sub.problem.constraints)
    end

    # storage constraints
    # initial conditions
    add_constraints!(dynProblem, [
        s[1] - INIT_COND <= P,
        -(s[1] - INIT_COND) <= P, 
        s[1] <= C, 
        0 <= s[1]
    ])
    # running condition
    for t in 2:T
        add_constraints!(dynProblem,[
            0 <= s[t], 
            s[t] <= C, 
            s[t] - s[t-1] <= P,
            s[t-1] - s[t] <= P
        ])
    end

    params = (fq=fq, fl=fl, d=d, pmax=pmax, gmax=gmax, A=A, B=B, P=P, C=C, τ=τ)

    return PowerManagementProblem(dynProblem, p_T, g_T, s, params)
end

DynamicPowerManagementProblem(net::DynamicPowerNetwork, d) =
    DynamicPowerManagementProblem(
        net.fq, net.fl, d, net.pmax, net.gmax, net.A, net.B, net.P, net.C; τ=net.τ
        )

# ===
# KKT OPERATOR for the DynPMP
# ===

"""
    kkt_dims_dyn(n, m,T)

Compute the dimensions of the input / output of the KKT operator for
a network with `n` nodes and `m` edges for `T` timesteps
"""
kkt_dims_dyn(n, m, T) = T*(9n+3m)

"""
    kkt(x, fq, fl, d, pmax, gmax, A; τ=TAU)

Compute the KKT operator applied to `x`, with parameters given by `fq`,
`fl`, `d`, `pmax`, `gmax`, `A`, and `τ`.


WIP/TODO:
---------
The steps should be:
- compute the KKTs for individual problems
- include the KKTs linked to S_t
- include the KKTs linked to the coupling between different timesteps
"""
function kkt_dyn(x, fq, fl, d, pmax, gmax, A, P, C; τ=TAU)

    n, m = size(A)
    T = length(fq)

    # Assume that it is giving back what we want for now
    # the format we assume is that each variable (both primal and dual)
    # is returned as an array of T elements
    # λdsi is the dual for charge/discharge of the battery, with i being 1 or 2
    g, p, s, λpl, λpu, λgl, λgu, λds1, λds2, λsl, λsu, ν = unflatten_variables_dyn(x, n, m, T)

    # compute the KKTs for individual problems
    KKT_tot = []
    for t in 1:T
        x = [g[t]; p[t]; λpl[t]; λpu[t]; λgl[t]; λgu[t]; ν[t]]
        KKT = kkt(x, fq[t], fl[t], d[t], pmax[t], gmax[t], A; τ=TAU)
        if t == 1
            ds = s[t] - INIT_COND
        else
            ds = s[t] - s[t-1]
        end
        KKT_s = kkt_storage(λds1[t], λds2[t], λsl[t], λsu[t], ds, s[t]) 
        KKT_tot = vcat(KKT_tot, KKT, KKT_s) # append the subproblem kkt matrix to the total KKT matrix
    end

    return KKT_tot
end

"""
Compute the terms kkt matrix that are only related to the storage
"""
function kkt_storage(λds1, λds2, λsl, λsu, ds, s)
    
    return [
        λds1 .* (ds - P);
        λds2 .* (-ds - P);
        λsl .* (-s);
        λsu .* (s - C)
    ]
end

"""
add doc
"""
function flatten_variables_dyn(P::PowerManagementProblem)
    x = [evaluate(P.g); evaluate(P.p)]
    λ = vcat([c.dual for c in P.problem.constraints]...)
    return [x; λ][:, 1]
end

"""
add doc
"""
function unflatten_variables_dyn(x, n, m, T)
    #TODO
end
