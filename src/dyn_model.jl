# Code for building and running a dynamic power management problem


# TODO: correct to include the proper size of `B` (and therefore of `g`)

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
        0 <= s[1],
        s[1] <= C, 
        s[1] - INIT_COND <= P,
        -(s[1] - INIT_COND) <= P
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

Details: 
- size of variables: T*(l + m + n)
- static subproblem constraints: (n + 2l + 2m)*T
- storage constraints: T*(4n)

"""
kkt_dims_dyn(n, m, l, T) = T*(6n+3m+3l)

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
function kkt_dyn(x, fq, fl, d, pmax, gmax, A, B, P, C; τ=TAU)

    n, m = size(A)
    _, l = size(B)
    T = length(fq)

    # Assume that it is giving back what we want for now
    # the format we assume is that each variable (both primal and dual)
    # is returned as an array of T elements
    # λdsi is the dual for charge/discharge of the battery, with i being 1 or 2
    g, p, s, λpl, λpu, λgl, λgu, ν, λsl, λsu, λdsl, λdsu = unflatten_variables_dyn(x, n, m, l, T)

    # compute the KKTs for individual problems
    KKT_tot = []
    for t in 1:T
        # extract primal/dual variables for the constraints of instance at time t
        # g, p, s, λpl, λpu, λgl, λgu, ν, λdsl, λdsu, λsl, λsu = extract_vars_t(P, t)
        x = [g[t]; p[t]; λpl[t]; λpu[t]; λgl[t]; λgu[t]; ν[t]]
        # compute the KKTs for the static subproblem
        KKT = kkt(x, fq[t], fl[t], d[t], pmax[t], gmax[t], A, B; τ=TAU)
        # add the KKts for the storage
        if t == 1
            ds = s[t] .- INIT_COND
            ν_prev = zeros(n)
            λdsl_prev = zeros(n)
            λdsu_prev = zeros(n)
        else
            ds = s[t] - s[t-1]
            ν_prev = ν[t-1]
            λdsl_prev = λdsl[t-1]
            λdsu_prev = λdsu[t-1]
        end
        KKT_s = kkt_storage(
            ν_prev, ν[t], λsl[t], λsu[t], λdsl_prev, λdsl[t], λdsu_prev, λdsu[t], ds, s[t], P, C
            ) 
        KKT_tot = vcat(KKT_tot, KKT, KKT_s) # append the subproblem kkt matrix to the total KKT matrix
    end

    @assert length(KKT_tot) == kkt_dims_dyn(n, m, l, T)

    return KKT_tot
end

kkt_dyn(x, net::DynamicPowerNetwork, d) =
    kkt_dyn(x, net.fq, net.fl, d, net.pmax, net.gmax, net.A, net.B, net.P, net.C; τ=net.τ)

"""
Compute the terms kkt matrix that are only related to the storage
"""
function kkt_storage(ν_prev, ν_t, λsl, λsu, λdsl_prev, λdsl_t, λdsu_prev, λdsu_t, ds, s, P, C)
    
    return [
        λsu - λsl + (λdsu_t-λdsl_t) - (λdsu_prev-λdsl_prev) + ν_t - ν_prev;
        λsl .* (-s);
        λsu .* (s - C);
        λdsl_t .* (-ds - P);
        λdsu_t .* (ds - P);
    ]
end

function extract_vars_t(P::PowerManagementProblem, t)

    n, m = size(P.params(A))
    _, l = size(P.params(B))

    T = length(P.g)

    g = P.g[t].value[:]
    p = P.p[t].value[:]
    s = P.s[t].value[:]
    # Constraints go in the following order: 
    # pl, pu, gl, gu, equality
    # ...
    # dsl, dsu, sl, su 
    start_index = (t-1)*5+1
    λpl = P.problem.constraints[start_index].dual  
    λpu = P.problem.constraints[start_index+1].dual 
    λgl = P.problem.constraints[start_index+2].dual
    λgu = P.problem.constraints[start_index+3].dual 
    ν = P.problem.constraints[start_index+4].dual  

    # checking size consistency of primal variables
    @assert length(λpl) == m
    @assert length(λpu) == m
    @assert length(λgl) == n
    @assert length(λgu) == n
    @assert length(ν) == n

    storage_index = 5*T + 1
    λdsl = P.problem.constraints[storage_index].dual
    λdsu = P.problem.constraints[storage_index+1].dual
    λsl = P.problem.constraints[storage_index+2].dual 
    λsu = P.problem.constraints[storage_index+3].dual 

    # checking size consistency of dual variables
    @assert length(λdsl) == n
    @assert length(λdsu) == n
    @assert length(λsl) == n
    @assert length(λsu) == n

    return g, p, s, λpl, λpu, λgl, λgu, ν, λdsl, λdsu, λsl, λsu
end

"""
    flatten_variables_dyn(P)

flattens the variables from a PMP

The structure of the resulting array should be:
- 
- 
- 
TODO: fill in above
"""
function flatten_variables_dyn(P::PowerManagementProblem)

    T = length(P.g)
    # extracting primal variables
    g = vcat([P.g[t].value[:] for t in 1:T]...)
    p = vcat([P.p[t].value[:] for t in 1:T]...)
    s = vcat([P.s[t].value[:] for t in 1:T]...)
    x = [g; p; s]

    # extracting dual variables
    # a. from the subproblems: 5 constraints - first 5T elements of 
    #    λ are for the static constraints
    # b. from storage: 4 constraints - last 4T elements of λ are for
    #    the storage coupling between timesteps
    λ = vcat([c.dual for c in P.problem.constraints]...)

    vars = [x; λ][:, 1]
    n, m = size(P.params.A)
    _, l = size(P.params.B)

    # checking that the size is consistent with expectations
    @assert length(vars) == kkt_dims_dyn(n, m, l, T)

    return vars
end

"""
add doc
"""
function unflatten_variables_dyn(x, n, m, l, T)

    # retrieving primal variables
    g = x[1:l*T]
    p = x[l*T+1: l*T+m*T]
    s = x[l*T+m*T+1: l*T+m*T+n*T]

    g = [g[(t-1)*l+1:t*l] for t in 1:T]
    p = [p[(t-1)*m+1:t*m] for t in 1:T]
    s = [s[(t-1)*n+1:t*n] for t in 1:T]

    # retrieving dual variables
    dual_start_index = n*T+m*T+l*T
    size_dual_static = (n+2m+2l)#(7n+2m)
    size_dual_storage = 4n

    λpl, λpu, λgl, λgu, ν, λdsl, λdsu, λsl, λsu = [], [], [], [], [], [], [], [] ,[]

    for t in 1:T

        crt_index = dual_start_index + (t-1)*size_dual_static
        # retrieving dual variables of the static problems
        λpl = vcat(λpl, [x[crt_index+1:crt_index+m]])
        λpu = vcat(λpu, [x[crt_index+m+1:crt_index+2m]])
        λgl = vcat(λgl, [x[crt_index+2m+1:crt_index+2m+n]])
        λgu = vcat(λgu, [x[crt_index+2m+n+1:crt_index+2m+2n]])
        ν = vcat(ν, [x[crt_index+2m+2n+1: crt_index + 2m+3n]])

        crt_index = dual_start_index + T*size_dual_static + (t-1) * size_dual_storage
        λsl = vcat(λsl, [x[crt_index+1: crt_index+n]])
        λsu = vcat(λsu, [x[crt_index+n+1: crt_index+2n]])
        λdsl = vcat(λdsl, [x[crt_index+2n+1: crt_index+3n]])
        λdsu = vcat(λdsu, [x[crt_index+3n+1: crt_index+4n]])

    end

    # sanity check for the sizes
    for t in 1:T
        @assert length(g[t]) == l
        @assert length(p[t]) == m
        @assert length(s[t]) == n
        @assert length(λpl[t]) == m
        @assert length(λpu[t]) == m
        @assert length(λgl[t]) == n
        @assert length(λgu[t]) == n
        @assert length(ν[t]) == n
        @assert length(λsl[t]) == n
        @assert length(λsu[t]) == n
        @assert length(λdsl[t]) == n
        @assert length(λdsu[t]) == n
    end

    return g, p, s, λpl, λpu, λgl, λgu, ν, λdsl, λdsu, λsl, λsu
end 

