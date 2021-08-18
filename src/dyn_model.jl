# Code for building and running a dynamic power management problem

# ===
# POWER MANAGEMENT PROBLEM
# ===

# TODO: 
# - add edge constraints

INIT_COND = 0

"""
A dynamic power network

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
    DynamicPowerManagementProblem(fq, fl, d, pmax, gmax, A, B, P, C; τ=TAU)

Define a PowerManagementProblem over `T` timesteps, where timesteps can be coupled
via the presence of nodal storage. Arguments are arrays of dimension `T`. `fq` and `fl`
are the quadratic and linear components of the generator costs, respectively. `d` is the 
nodal demand. `pmax`, `gmax` are the maximum transmission and generation power. 
`A` is the incidence matrix, `B` is a generator-to-node mapping. `P` and `C` are the maximum
charge/discharge power for batteries and the maximum SOC, respectively. 

Note:
-----
- Nodal storage variables are constrained to start at a value given by INIT_COND. 
- There is currently no final condition on the storage. 
"""
function DynamicPowerManagementProblem(
    fq, fl, d, pmax, gmax, A, B, P, C; τ=TAU
)
    T = length(fq)
    n, _ = size(A)

    # Define a storage variable for n nodes over T timesteps
    s = [Variable(n) for _ in 1:T]
    # Define a charge and discharge variable for n nodes over T timesteps
    ch = [Variable(n) for _ in 1:T]
    dis = [Variable(n) for _ in 1:T]

    #TODO: need to handle the initial condition
    subproblems = vcat(
        # first treating the initial constraint explicitly to avoid having to
        # specify another variable for s_0
        [PowerManagementProblem(
            fq[1], fl[1], d[1], pmax[1], gmax[1], A, B; ch=ch[1], dis=dis[1]
        )],
        # then iterating over all the timesteps 
        [PowerManagementProblem(
            fq[t], fl[t], d[t], pmax[t], gmax[t], A, B; ch=ch[t], dis=dis[t]
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
        0 <= s[1], # λsl
        s[1] <= C, # λsu
        ch[1] >= 0, #λchl
        ch[1] <= P, # λchu
        dis[1] >= 0, #λdisl
        dis[1] <= P, # λdisu
        s[1] == INIT_COND + ch[1] - dis[1], #νs (ν for storage)
    ])
    # running condition
    for t in 2:T
        add_constraints!(dynProblem,[
        0 <= s[t], # λsl
        s[t] <= C, # λsu
        ch[t] >= 0, #λchl
        ch[t] <= P, # λchu
        dis[t] >= 0, #λdisl
        dis[t] <= P, # λdisu
        0 == - s[t] + s[t-1] + ch[t] - dis[t], #νs (ν for storage)
        ])
    end

    params = (
        fq = fq, fl = fl, d = d, pmax = pmax, gmax = gmax, 
        A = A, B = B, P = P, C = C, τ = τ
    )

    return PowerManagementProblem(dynProblem, p_T, g_T, s, ch, dis, params)
end

DynamicPowerManagementProblem(net::DynamicPowerNetwork, d) =
    DynamicPowerManagementProblem(
        net.fq, net.fl, d, net.pmax, net.gmax, net.A, net.B, net.P, net.C; τ=net.τ
    )

# ===
# KKT OPERATOR for the DynPMP
# ===

"""
    kkt_dims_dyn(n, m, l, T)

Compute the dimensions of the input / output of the KKT operator for
a network with `n` nodes, `m` edges, `l` generators for `T` timesteps.

Details: 
--------
- size of variables: T*(l + m + n)
- static subproblem constraints: (n + 2l + 2m)*T
- storage constraints: T*(4n)
"""
kkt_dims_dyn(n, m, l, T) = T * (kkt_dims(n, m, l) + storage_kkt_dims(n))

"""
    storage_kkt_dims(n)

Compute the dimensions of the KKT operator associated with storage constraints.
"""
storage_kkt_dims(n) = 3n + 7n

"""
    kkt_dyn(x, fq, fl, d, pmax, gmax, A, B, P, C; τ=TAU)

Compute the KKT operator applied to `x`, with parameters given by `fq`,
`fl`, `d`, `pmax`, `gmax`, `A`, `B`, `P`, `C`, and `τ`.
"""
function kkt_dyn(x, fq, fl, d, pmax, gmax, A, B, P, C; τ=TAU)

    n, m = size(A)
    _, l = size(B)
    T = length(fq)

    # decompose `x` in arrays of T variables
    g, p, s, ch, dis, λpl, λpu, λgl, λgu, ν, λsl, λsu, λchl, λchu, λdisl, λdisu, νs = 
        unflatten_variables_dyn(x, n, m, l, T)

    # compute the KKTs for individual problems
    KKT_static_tot = Float64[]
    KKT_storage_tot = Float64[]
    for t in 1:T

        # extract primal/dual variables for the constraints of instance at time t
        x = [g[t]; p[t]; λpl[t]; λpu[t]; λgl[t]; λgu[t]; ν[t]]

        # handle edge cases
        # ?? TODO: figure out if there is an edge case for ch and dis? 
        # are they actually variables in n x T or n x (T-1)
        t==1 ? s_prev = INIT_COND : s_prev = s[t-1]
        t < T ? νs_next = νs[t + 1] : νs_next = zeros(n)

        # compute the KKTs for the static subproblem
        KKT = kkt(x, fq[t], fl[t], d[t], pmax[t], gmax[t], A, B; τ=TAU, ch=ch[t], dis=dis[t])

        # add the KKts for the storage
        KKT_s = kkt_storage(
            s[t], s_prev, ch[t], dis[t], λsu[t], λsl[t], λchu[t], λchl[t],  
            λdisu[t], λdisl[t], ν[t], νs[t], νs_next, P, C
            )

        # check the sizes of both matrices computed above
        @assert length(KKT) == kkt_dims(n, m, l)
        @assert length(KKT_s) == storage_kkt_dims(n)

        # append both matrices to the general KKT operator
        KKT_static_tot = vcat(KKT_static_tot, KKT) 
        KKT_storage_tot = vcat(KKT_storage_tot, KKT_s) 
    end

    KKT_tot = vcat(KKT_static_tot, KKT_storage_tot)
    @assert length(KKT_tot) == kkt_dims_dyn(n, m, l, T)

    return KKT_tot
end

kkt_dyn(x, net::DynamicPowerNetwork, d) = kkt_dyn(
    x, net.fq, net.fl, d, net.pmax, net.gmax, net.A, net.B, net.P, net.C; τ=net.τ
)

"""
Compute the terms in the kkt matrix that are only related to the storage

Notes:
------
- Dimensions should be 10n, as there are 3 variables of size n + 7 constraints of size n

"""
function kkt_storage(s, s_prev, ch, dis, λsu, λsl, λchu, λchl, λdisu, λdisl, ν, νs_t, νs_next, P, C)
    return [
        (λsu - λsl) + (νs_next - νs_t); #∇_s L
        (λchu - λchl) + ν/η_c + νs_t ; #∇_ch L
        (λdisu - λdisl) - η_d * ν - νs_t; #∇_dis L
        λsl .* (-s);
        λsu .* (s - C);
        λchl .* (-ch);
        λchu .* (ch-P);
        λdisl .* (-dis);
        λdisu .* (dis-P);
        s .- s_prev - ch + dis;
    ]
end

"""
    extract_vars_t(P::PowerManagementProblem, t)

Extracts both the primal and dual variables for a dynamic problem at timestep `t`.
"""
function extract_vars_t(P::PowerManagementProblem, t)

    n, m = size(P.params.A)
    _, l = size(P.params.B)
    T = length(P.g)
    
    n_constraints_static = 5
    n_constraints_storage = 7

    g = P.g[t].value[:]
    p = P.p[t].value[:]
    s = P.s[t].value[:]
    ch = P.ch[t].value[:]
    dis = P.dis[t].value[:]

    start_index = (t - 1) * n_constraints_static
    λpl = P.problem.constraints[start_index + 1].dual  
    λpu = P.problem.constraints[start_index + 2].dual 
    λgl = P.problem.constraints[start_index + 3].dual
    λgu = P.problem.constraints[start_index + 4].dual 
    ν = P.problem.constraints[start_index + 5].dual  

    # checking size consistency of primal variables
    @assert length(λpl) == m
    @assert length(λpu) == m
    @assert length(λgl) == l
    @assert length(λgu) == l
    @assert length(ν) == n

    storage_index = n_constraints_static * T + (t - 1) * n_constraints_storage
    λsl = P.problem.constraints[storage_index + 1].dual 
    λsu = P.problem.constraints[storage_index + 2].dual 
    λchl = P.problem.constraints[storage_index + 3].dual
    λchu = P.problem.constraints[storage_index + 4].dual
    λdisl = P.problem.constraints[storage_index + 5].dual
    λdisu = P.problem.constraints[storage_index + 6].dual
    νs = P.problem.constraints[storage_index + 7].dual

    # checking size consistency of dual variables
    @assert length(λsl) == n
    @assert length(λsu) == n
    @assert length(λchl) == n
    @assert length(λchu) == n
    @assert length(λdisl) == n
    @assert length(λdisu) == n
    @assert length(νs) == n

    return g, p, s, ch, dis, λpl, λpu, λgl, λgu, ν, λsl, λsu, λchl, λchu, λdisl, λdisu, νs
end

"""
    flatten_variables_dyn(P::PowerManagementProblem)

Flattens the variables from a `PowerManagementProblem`, i.e.
all variables over all timesteps are condensed into one array.

Details:
--------
The variables are laid out in the following order: 
- First variables from the static problem, i.e.
[g, p, λpl, λpu, λgl, λgu, ν] at timestep t. Variables 
from susbsequent timesteps are concatenated. 
- Then variables from the dynamic problem, i.e.
[s, λsl, λsu, λdsu, λdsl]. Variables from susbsequent
timesteps are concatenated.
"""
function flatten_variables_dyn(P::PowerManagementProblem)

    T = length(P.g)
    n, m = size(P.params.A)
    _, l = size(P.params.B)

    x_static = Float64[]
    x_storage = Float64[]

    # extract the variables at each timestep
    for t in 1:T
        g, p, s, ch, dis, λpl, λpu, λgl, λgu, ν, λsl, λsu, λchl, λchu, λdisl, λdisu, νs = extract_vars_t(P, t)
        x_static = vcat(x_static, g, p, λpl, λpu, λgl, λgu, ν)
        x_storage = vcat(x_storage, s, ch, dis, λsl, λsu, λchl, λchu, λdisl, λdisu, νs)
    end
    vars = vcat(x_static, x_storage)

    @assert length(x_static) == T * (3l + 3m + n) 
    @assert length(x_storage) == 10 * n * T

    # checking that the size is consistent with expectations
    @assert length(vars) == kkt_dims_dyn(n, m, l, T)

    return vars
end

"""
    unflatten_variables_dyn(x, n, m, l, T)

Extract the primal and dual variables from `x`, an array containing them all. 
`n`, `m`, `l`, `T` are the dimensions of the problem. 
"""
function unflatten_variables_dyn(x, n, m, l, T)

    static_length = 3l + 3m + n
    storage_length = 10n

    # variables
    g, p, s, ch, dis = [], [], [], [], []
    # static constriaints
    λpl, λpu, λgl, λgu, ν = [], [], [], [], []
    # dynamic constraints
    λsl, λsu, λchl, λchu, λdisl, λdisu, νs = [], [], [], [], [], [], []

    # loop through timesteps, and retrieve variables in the order that 
    # x is constructed in flatten_variables_dyn()
    # Step 1: going through the static variables
    crt_idx, i = 0, 0
    for t in 1:T
        crt_idx = static_length * (t - 1)
        i = 0
        g = vcat(g, [x[crt_idx + i + 1:crt_idx + i + l]])
        i += l
        p = vcat(p, [x[crt_idx + i + 1:crt_idx + i + m]])
        i += m
        λpl = vcat(λpl, [x[crt_idx + i + 1:crt_idx + i + m]])
        i += m
        λpu = vcat(λpu, [x[crt_idx + i + 1:crt_idx + i + m]])
        i += m
        λgl = vcat(λgl, [x[crt_idx + i + 1:crt_idx + i + l]])
        i += l
        λgu = vcat(λgu, [x[crt_idx + i + 1:crt_idx + i + l]])
        i += l
        ν = vcat(ν, [x[crt_idx + i + 1:crt_idx + i + n]])
        i += n
    end
    # check that we have gone through all the static subproblems
    @assert crt_idx + i == T * static_length
    # Step 2: going throught the storage variables
    for t in 1:T
        crt_idx = T * static_length + (t - 1) * storage_length
        i = 0
        s = vcat(s, [x[crt_idx + i + 1:crt_idx + i + n]])
        i += n
        ch = vcat(ch, [x[crt_idx + i + 1:crt_idx + i + n]])
        i += n
        dis = vcat(dis, [x[crt_idx + i + 1:crt_idx + i + n]])
        i += n
        λsl = vcat(λsl, [x[crt_idx + i + 1:crt_idx + i + n]])
        i += n
        λsu = vcat(λsu, [x[crt_idx + i + 1:crt_idx + i + n]])
        i += n
        λchl = vcat(λchl, [x[crt_idx + i + 1:crt_idx + i + n]])
        i += n
        λchu = vcat(λchu, [x[crt_idx + i + 1:crt_idx + i + n]])
        i += n
        λdisl = vcat(λdisl, [x[crt_idx + i + 1:crt_idx + i + n]])
        i += n
        λdisu = vcat(λdisu, [x[crt_idx + i + 1:crt_idx + i + n]])
        i += n
        νs = vcat(νs, [x[crt_idx + i + 1:crt_idx + i + n]])
        i += n
    end
    @assert crt_idx + i == T * (static_length + storage_length)

    # sanity check for the sizes
    for t in 1:T
        @assert length(g[t]) == l
        @assert length(p[t]) == m
        @assert length(s[t]) == n
        @assert length(ch[t]) == n
        @assert length(dis[t]) == n
        @assert length(λpl[t]) == m
        @assert length(λpu[t]) == m
        @assert length(λgl[t]) == l
        @assert length(λgu[t]) == l
        @assert length(ν[t]) == n
        @assert length(λsl[t]) == n
        @assert length(λsu[t]) == n
        @assert length(λchl[t]) == n
        @assert length(λchu[t]) == n
        @assert length(λdisl[t]) == n
        @assert length(λdisu[t]) == n
        @assert length(νs[t]) == n
    end

    return g, p, s, ch, dis, λpl, λpu, λgl, λgu, ν, λsl, λsu, λchl, λchu, λdisl, λdisu, νs
end 

