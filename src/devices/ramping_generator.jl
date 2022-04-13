# Note: some behavior is inherited from the AbstractGenerator type

struct RampGenerator <: AbstractGenerator
    pmin::Vector
    pmax::Vector
    cost::Vector
    ramp
end

function with_ramping(d::Generator, ramp)
    return RampGenerator(d.pmin, d.pmax, d.cost, ramp)
end

# Device data

make_constraints(d::RampGenerator, p) = [
    d.pmin <= p, 
    p <= d.pmax, 
    p[2:end] - p[1:end-1] <= d.ramp,
    p[2:end] - p[1:end-1] >= -d.ramp
]

get_dims(d::RampGenerator)::Int = 3*get_time_horizon(d) + 2*(get_time_horizon(d)-1)

# KKTs

function get_device_dL_dx(d::RampGenerator, p, duals)
    c = d.cost
    μ_lower, μ_upper, ψ_upper, ψ_lower = duals

    μ = μ_upper - μ_lower

    # TODO: perhaps the `diff` operator would come in handy here?
    ψ = [0; ψ_upper] - [ψ_upper; 0] - [0; ψ_lower] + [ψ_lower; 0]
    
    return c + μ + ψ
end

function get_device_comp_slack(d::RampGenerator, p, duals)
    (; pmin, pmax, ramp) = d
    μ_lower, μ_upper, ψ_upper, ψ_lower = duals
    T = get_time_horizon(d)

    fx = [
        pmin - p;
        p - pmax;
        p[2:T] - p[1:T-1] .- ramp;
        -p[2:T] + p[1:T-1] .- ramp;
    ]
    
    return reduce(vcat, duals) .* fx
end

# Jacobian

function jacobian_kkt_z_device(d::RampGenerator, p, duals)
    (; pmin, pmax, ramp) = d
    μ_lower, μ_upper, ψ_upper, ψ_lower = duals
    T  = get_time_horizon(d)

    I_lag = I(T)[2:T, 1:T]
    I_lead = I(T)[1:T-1, 1:T]

    ψ_upper, ψ_lower = ψ_upper[:, 1], ψ_lower[:, 1]

    Dfx = [
        -I(T)
        I(T)
        I_lag - I_lead
        I_lead - I_lag
    ]

    fx = [
        pmin - p;
        p - pmax;
        p[2:T] - p[1:T-1] .- ramp;
        -p[2:T] + p[1:T-1] .- ramp;
    ]

    Dλ = Diagonal(reduce(vcat, duals)[:, 1])

    return [
        zeros(T, T) Dfx'
        Dλ*Dfx Diagonal(fx)
    ]
end
