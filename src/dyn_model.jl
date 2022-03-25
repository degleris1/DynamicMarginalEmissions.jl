# Code for building and running a dynamic power management problem

# ===
# POWER MANAGEMENT PROBLEM
# ===

# TODO: incorporate as kwarg? 
INIT_COND = .5
FINAL_COND = .5

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
η_c: discharge inefficiency
η_d: charge inefficiency
ρ: max ramping rate
"""
mutable struct DynamicPowerNetwork
    fq
    fl
    pmax
    gmax
    A
    B
    F
    S
    P
    C
    T
    τ
    η_c
    η_d
    ρ
end

function DynamicPowerNetwork(
    fq, fl, pmax, gmax, A, B, F, S, P, C, T; 
    τ=TAU, η_c=1.0, η_d=1.0, ρ=nothing
)
    # The generator will never change by more than `gmax` MW per timestep
    # Therefore, if no ramping constraint is specified, we can set the
    # ramping limit to ρ > gmax, and the constraint will never be binding
    ρ = something(ρ, 2*[maximum(x -> x[i], gmax) for i in 1:length(gmax[1])])
    return DynamicPowerNetwork(fq, fl, pmax, gmax, A, B, F, S, P, C, T, τ, η_c, η_d, ρ)
end


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
charge/discharge power for batteries and the maximum SOC, respectively. `S` maps storage to associated nodes.

Note:
-----
- Nodal storage variables are constrained to start at a value given by INIT_COND. 
- There is currently no final condition on the storage. 
"""
function DynamicPowerManagementProblem(
    fq, fl, d, pmax, gmax, A, B, F, S, P, C; 
    τ=TAU, η_c=1.0, η_d=1.0, ρ=nothing
)
    ρ = something(ρ, 2*[maximum(x -> x[i], gmax) for i in 1:length(gmax[1])])

    T = length(fq)
    n, _ = size(A)
    _, ns = size(S) #ns indicates the number of storage nodes

    # Define a storage variable for n nodes over T timesteps
    s = [Variable(ns) for _ in 1:T-1]
    # Define a charge and discharge variable for n nodes over T timesteps
    ch = [Variable(ns) for _ in 1:T]
    dis = [Variable(ns) for _ in 1:T]

    subproblems = vcat(
        [PowerManagementProblem(
            fq[t], fl[t], d[t], pmax[t], gmax[t], A, B, F; 
            τ=τ, ch=ch[t], dis=dis[t], S=S, η_c=η_c, η_d=η_d
        ) for t in 1:T]
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
        g_T[1] >= -ρ,  # λrampl
        g_T[1] <= ρ,  # λrampu
        0 == - s[1] + INIT_COND.*C + ch[1] * η_c - dis[1]/η_d, #νs (ν for storage)
    ])
    # running condition
    for t in 2:T-1
        add_constraints!(dynProblem,[
            0 <= s[t], # λsl
            s[t] <= C, # λsu
            ch[t] >= 0, #λchl
            ch[t] <= P, # λchu
            dis[t] >= 0, #λdisl
            dis[t] <= P, # λdisu
            g_T[t] >= g_T[t-1] - ρ,  # λrampl
            g_T[t] <= g_T[t-1] + ρ,  # λrampu
            0 == - s[t] + s[t-1] + ch[t] * η_c - dis[t]/η_d, #νs (ν for storage)
        ])
    end
    # final condition
    add_constraints!(dynProblem,[
            # the first two constraints have to be eliminated
            # because s[T] is not instantiated (dim(s) = T-1)
            # 0 <= s[t], # λsl
            # s[t] <= C, # λsu
            ch[T] >= 0, #λchl
            ch[T] <= P, # λchu
            dis[T] >= 0, #λdisl
            dis[T] <= P, # λdisu
            g_T[T] >= g_T[T-1] - ρ,  # λrampl
            g_T[T] <= g_T[T-1] + ρ,  # λrampu
            0 == - FINAL_COND.*C + s[T-1] + ch[T] * η_c - dis[T]/η_d, #νs (ν for storage)
        ])

    params = (
        fq = fq, fl = fl, d = d, pmax = pmax, gmax = gmax, 
        A = A, B = B, F = F, S = S, P = P, C = C, τ = τ
    )

    return PowerManagementProblem(dynProblem, p_T, g_T, s, ch, dis, params)
end

DynamicPowerManagementProblem(net::DynamicPowerNetwork, d) =
    DynamicPowerManagementProblem(
        net.fq, net.fl, d, net.pmax, net.gmax, net.A, net.B, net.F, net.S, net.P, net.C; 
        τ=net.τ, η_c=net.η_c, η_d=net.η_d, ρ=net.ρ
    )

# ===
# KKT OPERATOR for the DynPMP
# ===

"""
    kkt_dims_dyn(n, m, l, ns, T)

Compute the dimensions of the input / output of the KKT operator for
a network with `n` nodes, `m` edges, `l` generators for `T` timesteps.
`ns` is the number of storage nodes.

Details: 
--------
- size of variables in static problem: T*(l + m + n)
- static subproblem constraints: (n + 2l + 2m)*T
- storage constraints: T*(4ns)
"""
kkt_dims_dyn(n, m, l, ns, T) = T * (kkt_dims(n, m, l) + storage_kkt_dims(ns, l))

"""
    storage_kkt_dims(n, l)

Compute the dimensions of the KKT operator associated with storage constraints.
"""
storage_kkt_dims(ns, l) = 3ns + 7ns + 2l # 3ns for stationarity, 7ns+2l for complementary slackness

"""
    kkt_dyn(x, fq, fl, d, pmax, gmax, A, B, P, C, η_c, η_d, ρ; τ=TAU)

Compute the KKT operator applied to `x`, with parameters given by `fq`,
`fl`, `d`, `pmax`, `gmax`, `A`, `B`, `P`, `C`, and `τ`.
"""
function kkt_dyn(x, fq, fl, d, pmax, gmax, A, B, F, S, P, C, η_c, η_d, ρ; τ=TAU)

    n, m = size(A)
    _, l = size(B)
    _, ns = size(S)
    T = length(fq)

    # decompose `x` in arrays of T variables
    g, p, s, ch, dis, λpl, λpu, λgl, λgu, ν, νE, λsl, λsu, λchl, λchu, λdisl, λdisu, λrampl, λrampu, νs = 
        unflatten_variables_dyn(x, n, m, l, ns, T)

    # compute the KKTs for individual problems
    KKT_static_tot = Float64[]
    KKT_storage_tot = Float64[]
    for t in 1:T

        # extract primal/dual variables for the constraints of instance at time t
        x = [g[t]; p[t]; λpl[t]; λpu[t]; λgl[t]; λgu[t]; ν[t]; νE[t]]

        # handle edge cases
        # ?? TODO: figure out if there is an edge case for ch and dis? 
        # are they actually variables in n x T or n x (T-1)
        t == 1 ? s_prev = INIT_COND.*C : s_prev = s[t-1]
        t == T ? s_crt = FINAL_COND.*C : s_crt = s[t]
        t < T ? νs_next = νs[t + 1] : νs_next = zeros(ns)

        # compute the KKTs for the static subproblem
        λrampl_next, λrampu_next = (t == T) ? (zeros(l), zeros(l)) : (λrampl[t+1], λrampu[t+1])
        KKT = (
            kkt(
                x, fq[t], fl[t], d[t], pmax[t], gmax[t], A, B, F; 
                τ=TAU, ch=ch[t], dis=dis[t], S=S
            )
            + kkt_ramp(n, m, l, λrampl_next, λrampu_next, λrampl[t], λrampu[t])
        )

        # add the KKts for the storage
        gt_prev = (t == 1) ? zeros(l) : g[t-1]
        KKT_s = kkt_storage(
            s_crt, s_prev, ch[t], dis[t], λsu[t], λsl[t], λchu[t], λchl[t],
            λdisu[t], λdisl[t], λrampl[t], λrampu[t], ν[t], νE[t], νs[t], νs_next, F, S, P, C, 
            η_c, η_d, ρ, g[t], gt_prev,
        )

        # check the sizes of both matrices computed above
        @assert length(KKT) == kkt_dims(n, m, l)
        @assert length(KKT_s) == storage_kkt_dims(ns, l)

        # append both matrices to the general KKT operator
        KKT_static_tot = vcat(KKT_static_tot, KKT) 
        KKT_storage_tot = vcat(KKT_storage_tot, KKT_s) 
    end

    KKT_tot = vcat(KKT_static_tot, KKT_storage_tot)
    @assert length(KKT_tot) == kkt_dims_dyn(n, m, l, ns,T)

    return KKT_tot
end

"""
    kkt_dyn(x, net, d)
"""
kkt_dyn(x, net::DynamicPowerNetwork, d) = kkt_dyn(
    x, net.fq, net.fl, d, net.pmax, net.gmax, net.A, net.B, net.F, net.S, 
    net.P, net.C, net.η_c, net.η_d, net.ρ; τ=net.τ
)

"""
Compute the terms in the kkt matrix that are only related to the storage

Notes:
------
- Dimensions should be 10n, as there are 3 variables of size n + 7 constraints of size n
- The variables that are not in the diagonal block for the Jacobian are νs_next and ν
"""
function kkt_storage(
    s, s_prev, ch, dis, 
    λsu, λsl, λchu, λchl, λdisu, λdisl, λrampl, λrampu, ν, νE, νs_t, νs_next, 
    F, S, P, C, η_c, η_d, ρ, gt, gt_prev
)
    return [
        (λsu - λsl) + (νs_next - νs_t);  # ∇_s L, dim = ns
        (λchu - λchl) - (F*S)'*ν .+ νE + νs_t * η_c ;  # ∇_ch L, dim = ns
        (λdisu - λdisl) + (F*S)'*ν .- νE - νs_t/η_d;  # ∇_dis L, dim = ns
        λsl .* (-s); # dim = ns
        λsu .* (s - C); # dim = ns
        λchl .* (-ch); # dim = ns
        λchu .* (ch-P); # dim = ns
        λdisl .* (-dis); # dim = ns
        λdisu .* (dis-P); # dim = ns
        λrampl .* (gt_prev - ρ - gt); # dim = l 
        λrampu .* (gt - gt_prev - ρ); # dim = l
        - s .+ s_prev + ch * η_c - dis/η_d; # dim = ns
    ]
end

function kkt_ramp(n, m, l, λrampl_next, λrampu_next, λrampl, λrampu)
    # gt_prev - ρ - gt <= 0,  # λrampl
    # gt - gt_prev - ρ <= 0,  # λrampu
    K = [-λrampl + λrampu + λrampl_next - λrampu_next; spzeros(kkt_dims(n, m, l)-l)]
end

"""
    extract_vars_t(P::PowerManagementProblem, t)

Extracts both the primal and dual variables for a dynamic problem at timestep `t`.
"""
function extract_vars_t(P::PowerManagementProblem, t)

    n, m = size(P.params.A)
    _, l = size(P.params.B)
    _, ns = size(P.params.S)
    T = length(P.g)
    
    n_constraints_static = 6
    n_constraints_storage = 9

    g = P.g[t].value[:]
    p = evaluate(P.p[t]) 

    if t==T
        s = FINAL_COND.*P.params.C
    else
        s = evaluate(P.s[t]) 
    end

    ch = evaluate(P.ch[t]) 
    dis = evaluate(P.dis[t])

    make_array = x -> (typeof(x) <: Array) ? x : [x]

    p = make_array(p)
    s = make_array(s)
    ch = make_array(ch)
    dis = make_array(dis)

    start_index = (t - 1) * n_constraints_static
    λpl = P.problem.constraints[start_index + 1].dual  
    λpu = P.problem.constraints[start_index + 2].dual 
    λgl = P.problem.constraints[start_index + 3].dual
    λgu = P.problem.constraints[start_index + 4].dual
    ν = P.problem.constraints[start_index + 5].dual
    νE = P.problem.constraints[start_index + 6].dual

    # checking size consistency of primal variables
    @assert length(λpl) == m
    @assert length(λpu) == m
    @assert length(λgl) == l
    @assert length(λgu) == l
    @assert length(ν) == m
    @assert length(νE) == 1

    storage_index = n_constraints_static * T + (t - 1) * n_constraints_storage
    # two constraints are eliminated for the final step
    if t==T
        storage_index = storage_index - 2
    end
    λsl = (t==T) ? zeros(ns) : P.problem.constraints[storage_index + 1].dual 
    λsu = (t==T) ? zeros(ns) : P.problem.constraints[storage_index + 2].dual 
    λchl = P.problem.constraints[storage_index + 3].dual
    λchu = P.problem.constraints[storage_index + 4].dual
    λdisl = P.problem.constraints[storage_index + 5].dual
    λdisu = P.problem.constraints[storage_index + 6].dual
    λrampl = P.problem.constraints[storage_index + 7].dual
    λrampu = P.problem.constraints[storage_index + 8].dual
    νs = P.problem.constraints[storage_index + 9].dual

    # checking size consistency of dual variables
    @assert length(λsl) == ns
    @assert length(λsu) == ns
    @assert length(λchl) == ns
    @assert length(λchu) == ns
    @assert length(λdisl) == ns
    @assert length(λdisu) == ns
    @assert length(λrampl) == l
    @assert length(λrampu) == l
    @assert length(νs) == ns

    return g, p, s, ch, dis, λpl, λpu, λgl, λgu, ν, νE, λsl, λsu, λchl, λchu, λdisl, λdisu, λrampl, λrampu, νs
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
    _, ns = size(P.params.S)
    _, l = size(P.params.B)

    x_static = Float64[]
    x_storage = Float64[]

    # extract the variables at each timestep
    for t in 1:T
        g, p, s, ch, dis, λpl, λpu, λgl, λgu, ν, νE, λsl, λsu, λchl, λchu, λdisl, λdisu, λrampl, λrampu, νs = extract_vars_t(P, t)
        x_static = vcat(x_static, g, p, λpl, λpu, λgl, λgu, ν, νE)
        x_storage = vcat(x_storage, s, ch, dis, λsl, λsu, λchl, λchu, λdisl, λdisu, λrampl, λrampu, νs)
    end
    vars = vcat(x_static, x_storage)

    @assert length(x_static) == T*kkt_dims(n, m, l)
    @assert length(x_storage) == T*storage_kkt_dims(ns, l)

    # checking that the size is consistent with expectations
    @assert length(vars) == kkt_dims_dyn(n, m, l, ns, T)

    return vars
end

"""
    unflatten_variables_dyn(x, n, m, l, T)

Extract the primal and dual variables from `x`, an array containing them all. 
`n`, `m`, `l`, `T` are the dimensions of the problem. 
"""
function unflatten_variables_dyn(x, n, m, l, ns, T)

    static_length = kkt_dims(n, m, l)
    storage_length = storage_kkt_dims(ns, l)

    # variables
    g, p, s, ch, dis = [], [], [], [], []
    # static constriaints
    λpl, λpu, λgl, λgu, ν, νE = [], [], [], [], [], []
    # dynamic constraints
    λsl, λsu, λchl, λchu, λdisl, λdisu, λrampl, λrampu, νs = [], [], [], [], [], [], [], [], []

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
        ν = vcat(ν, [x[crt_idx + i + 1:crt_idx + i + m]])
        i += m
        νE = vcat(νE, [x[crt_idx + i + 1:crt_idx + i + 1]])
        i += 1
    end
    # check that we have gone through all the static subproblems
    @assert crt_idx + i == T * static_length

    # Step 2: going throught the storage variables
    for t in 1:T
        crt_idx = T * static_length + (t - 1) * storage_length
        i = 0
        s = vcat(s, [x[crt_idx + i + 1:crt_idx + i + ns]])
        i += ns
        ch = vcat(ch, [x[crt_idx + i + 1:crt_idx + i + ns]])
        i += ns
        dis = vcat(dis, [x[crt_idx + i + 1:crt_idx + i + ns]])
        i += ns
        λsl = vcat(λsl, [x[crt_idx + i + 1:crt_idx + i + ns]])
        i += ns
        λsu = vcat(λsu, [x[crt_idx + i + 1:crt_idx + i + ns]])
        i += ns
        λchl = vcat(λchl, [x[crt_idx + i + 1:crt_idx + i + ns]])
        i += ns
        λchu = vcat(λchu, [x[crt_idx + i + 1:crt_idx + i + ns]])
        i += ns
        λdisl = vcat(λdisl, [x[crt_idx + i + 1:crt_idx + i + ns]])
        i += ns
        λdisu = vcat(λdisu, [x[crt_idx + i + 1:crt_idx + i + ns]])
        i += ns
        λrampl = vcat(λrampl, [x[crt_idx + i + 1:crt_idx + i + l]])
        i += l
        λrampu = vcat(λrampu, [x[crt_idx + i + 1:crt_idx + i + l]])
        i += l
        νs = vcat(νs, [x[crt_idx + i + 1:crt_idx + i + ns]])
        i += ns
    end
    @assert crt_idx + i == T * (static_length + storage_length)

    # sanity check for the sizes
    for t in 1:T
        @assert length(g[t]) == l
        @assert length(p[t]) == m
        @assert length(s[t]) == ns
        @assert length(ch[t]) == ns
        @assert length(dis[t]) == ns
        @assert length(λpl[t]) == m
        @assert length(λpu[t]) == m
        @assert length(λgl[t]) == l
        @assert length(λgu[t]) == l
        @assert length(ν[t]) == m
        @assert length(νE[t]) == 1
        @assert length(λsl[t]) == ns
        @assert length(λsu[t]) == ns
        @assert length(λchl[t]) == ns
        @assert length(λchu[t]) == ns
        @assert length(λdisl[t]) == ns
        @assert length(λdisu[t]) == ns
        @assert length(λrampl[t]) == l
        @assert length(λrampu[t]) == l
        @assert length(νs[t]) == ns
    end

    return g, p, s, ch, dis, λpl, λpu, λgl, λgu, ν, νE, λsl, λsu, λchl, λchu, λdisl, λdisu, λrampl, λrampu, νs
end 

