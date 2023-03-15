"""
    dispatch_and_compute_emissions_rates(; kwargs...)

Set up and solve the electricity dispatch problem, then compute locational
marginal emissions rates. Below, we refer to the following variables.

- `T` : The number of time periods, same as `time_horizon`
- `n` : The number of nodes
- `m` : The number of transmission lines
- `k` : The number of generators
- `s` : The number of batteries

# Return Value

Return a `NamedTuple` with keys `(lmes, g, meta)`.

- `lmes::Matrix` with size `(n, T)`

    The entry `lmes[i, t]` is the marginal emissions rate at node `i` and 
    time `t`.

- `g::Matrix` with size `(n, T)`

    The entry `g[k, t]` is the power produced by generator `k` at time `t`.

- `meta::Dict`

    Contains the `DynamicPowerProblem` object used to solve the dispatch 
    problem and various debugging stats (e.g., solve time).

# Arguments

- `time_horizon::Int`
    
    The number of time periods in the problem.

- `demand_schedule::Vector{Vector}` with size `(T,)`

    Each entry `demand_schedule[t]` contains the `n` nodal demands for time `t`.

- `node_line_incidence_matrix::Matrix` with size `(n, m)`

    The incidence matrix mapping transmission lines to nodes.

- `node_generator_incidence_matrix::Matrix` with size `(n, k)`

    The incidence matrix mapping generators to nodes.

- `node_battery_incidence_matrix::Matrix` with size `(n, s)`

    The incidence matrix mapping batteries to nodes.

- `generator_costs_quadratic::Vector{Vector}` with size `(T,)`

    Each entry `generator_costs_quadratic[t]` contains the `k` generator quadratic cost 
    coefficients for time `t`.

- `generator_costs_linear::Vector{Vector}` with size `(T,)`

    Each entry `generator_costs_linear[t]` contains the `k` generator linear cost 
    coefficients for time `t`.

- `generator_capacities::Vector{Vector}` with size `(T,)`

    Each entry `generator_capacities[t]` contains the `k` generator capacities
    for time `t`.

- `generator_emissions_rates::Vector{<:Real}` with size `(k,)`

    The `k` generator emissions rates.

- `generator_ramping_rates::Union{Vector{<:Real}, Nothing} = nothing` with size `(k,)` (optional)

    The `k` generator ramping rates.

- `line_susceptances::Vector{<:Real}` with size `(m,)`

    The `m` line electrical susceptances (inverse of impedance).

- `line_capacities::Vector{Vector}` with size `(T,)`

    Each entry `line_capacities[t]` contains the `m` line capacities
    for time `t`.

- `battery_capacities::Vector` with size `(s,)`

    The maximum capacity of each battery.

- `battery_max_powers::Vector` with size `(s,)`

    The maximum charging or discharging rate of each battery.

- `battery_charge_efficiency::Real = 1.0` (optional)

    The efficiency of battery charging.

- `battery_discharge_efficiency::Real = 1.0` (optional)

    The efficiency of battery discharging.

- `solver=ECOS.Optimizer` (optional)

    The optimization solver used to solve the dispatch problem.  

"""
function dispatch_and_compute_emissions_rates(
    ;
    time_horizon::Int,
    demand_schedule,

    node_line_incidence_matrix,
    node_generator_incidence_matrix,
    node_battery_incidence_matrix,

    generator_costs_quadratic,
    generator_costs_linear,
    generator_capacities,
    generator_emissions_rates,
    generator_ramping_rates=nothing,

    line_susceptances,
    line_capacities,

    battery_capacities,
    battery_max_powers,
    battery_charge_efficiency::Real=1.0,
    battery_discharge_efficiency::Real=1.0,

    solver=ECOS_OPT,
)

    # Set up network and problem
    pfdf_matrix = make_pfdf_matrix(node_line_incidence_matrix, line_susceptances)

    network = DynamicPowerNetwork(
        generator_costs_quadratic, 
        generator_costs_linear, 
        line_capacities, 
        generator_capacities, 
        node_line_incidence_matrix, 
        node_generator_incidence_matrix, 
        pfdf_matrix, 
        node_battery_incidence_matrix,
        battery_max_powers, 
        battery_capacities, 
        time_horizon; 
        η_c=battery_charge_efficiency, 
        η_d=battery_discharge_efficiency, 
        ρ=generator_ramping_rates,
    )

    power_problem = DynamicPowerManagementProblem(network, demand_schedule)

    # Solve problem
    solve_time = @elapsed solve!(power_problem, solver)

    # Evaluate solution
    g = CarbonNetworks.evaluate.(power_problem.g)

    # Compute locational marginal emissions rates
    differentiation_time = @elapsed λ = compute_mefs(
        power_problem,
        network, 
        demand_schedule, 
        generator_emissions_rates
    )

    # Fill LMEs in a tensor
    n = size(node_line_incidence_matrix, 1)
    T = time_horizon

    # `dynamic_lmes[i, t_affected, t]` defines how a change in demand at node 
    # `i` and time `t` will affect emissions at time `t_affected`. Summing 
    # the middle dimensions should produce the total emissions affect of 
    # demand at time `t`.
    dynamic_lmes = zeros(n, T, T)

    for t in 1:T
        dynamic_lmes[:, :, t] .= λ[t];
    end

    # Aggregrate LMEs to find total emissions affect
    lmes = zeros(n, T)
    
    for (i, t) in Iterators.product(1:n, 1:T)
        lmes[i, t] = sum(dynamic_lmes[i, :, t])
    end

    # Populate metadata
    meta = Dict(
        :network => network,
        :power_problem => power_problem,
        :solve_time => solve_time,
        :differentiation_time => differentiation_time,
        :full_lmes => dynamic_lmes,
    )

    return (lmes=lmes, g=g, meta=meta)
end