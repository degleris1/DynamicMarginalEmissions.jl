struct Generator <: AbstractDevice
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

get_cost(d::Generator, p) = d.cost' * p

get_constraints(d::Generator, p) = [d.pmin <= p, p <= d.pmax]

# Dimensions, indices

get_time_horizon(d::Generator) = length(d.pmin)

get_dims(d::Generator)::Int = 3*get_time_horizon(d)

# KKTs

function get_device_dL_dx(d::Generator, p, duals)
    c = d.cost
    μ_lower, μ_upper = duals
    
    return c - μ_lower + μ_upper
end

function get_device_comp_slack(d::Generator, p, duals)
    (; pmin, pmax) = d
    μ_lower, μ_upper = duals
    
    return [
        μ_lower .* (pmin - p);
        μ_upper .* (p - pmax)
    ]
end

# Jacobian

function jacobian_kkt_z_device(d::Generator, p, duals)
    μ_lower, μ_upper = duals
    T  = get_time_horizon(d)
    
    return [
        zeros(T, T) -I(T) I(T);
        -diagm(μ_lower[:, 1]) -diagm(p) zeros(T, T);
        diagm(μ_upper[:, 1]) zeros(T, T) diagm(p);
    ]
end
