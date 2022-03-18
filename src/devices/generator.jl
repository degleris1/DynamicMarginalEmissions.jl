abstract type AbstractGenerator <: AbstractDevice end

struct Generator <: AbstractGenerator
    pmin::Vector
    pmax::Vector
    cost::Vector
end

function StaticGenerator(pmin::Real, pmax::Real, cost::Real, T::Int)
    return Generator(pmin*ones(T), pmax*ones(T), cost*ones(T))
end

function Demand(d::Vector, curtailment_cost::Real)
    T = length(d)
    return Generator(-d, zeros(T), curtailment_cost*ones(T))
end

# Device Data

get_cost(d::AbstractGenerator, p) = d.cost' * p

make_aux_vars(d::AbstractGenerator) = nothing

make_constraints(d::AbstractGenerator, p) = [d.pmin <= p, p <= d.pmax]

# Dimensions, indices

get_time_horizon(d::AbstractGenerator) = length(d.pmin)

get_dims(d::AbstractGenerator)::Int = 3*get_time_horizon(d)

# KKTs

function get_device_dL_dx(d::AbstractGenerator, p, duals)
    c = d.cost
    μ_lower, μ_upper = duals
    
    return c - μ_lower + μ_upper
end

function get_device_comp_slack(d::AbstractGenerator, p, duals)
    (; pmin, pmax) = d
    μ_lower, μ_upper = duals
    
    return [
        μ_lower .* (pmin - p);
        μ_upper .* (p - pmax)
    ]
end

# Jacobian

function jacobian_kkt_z_device(d::AbstractGenerator, p, duals)
    (; pmin, pmax) = d
    μ_lower, μ_upper = duals
    T  = get_time_horizon(d)
    
    return [
        zeros(T, T) -I(T) I(T);
        -diagm(μ_lower[:, 1]) diagm(pmin - p) zeros(T, T);
        diagm(μ_upper[:, 1]) zeros(T, T) diagm(p - pmax);
    ]
end
