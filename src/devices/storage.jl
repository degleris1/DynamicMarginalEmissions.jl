
struct Battery <: AbstractDevice
    T::Int  # Time horizon
    C::Real  # Capacity
    R::Real  # Power
    s0::Real  # Initial SOC
    sT::Real  # Final SOC
    η::Real  # Efficiency
end


# Device Data

get_cost(d::Battery, p) = 0

function make_aux_vars(d::Battery)
    (; T, C, R, s0, sT, η) = d

    # Charge and discharge
    ch = Variable(T)
    ds = Variable(T)

    # Storage at the end of each time period
    s = Variable(T)

    return (ch=ch, ds=ds, s=s)
end

function make_constraints(d::Battery, p, aux_vars)
    (; T, C, R, s0, sT, η) = d
    (; ch, ds, s) = aux_vars

    return [
        ch >= 0,
        ch <= R,
        ds >= 0,
        ds <= R,
        s >= 0,
        s <= C,
        0 == -p + sqrt(η)*ds - (1/sqrt(η))*ch,
        0 == -s - ds + ch + [s0; s[1:T-1]],
        0 == -s[T] + sT
    ]
end

# Dimensions, indices

get_time_horizon(d::Battery) = d.T

function get_dims(d::Battery)::Int
    T = get_time_horizon(d)

    num_vars = 4T
    num_ineq = 6T
    num_eq = 2T + 1

    return num_vars + num_ineq + num_eq
end

# KKTs

function get_device_dL_dx(d::Battery, p, duals, aux_vars)
    νp = duals[7]  
    return -νp
end

function get_device_dL_daux(d::Battery, p, duals, aux_vars)
    λcl, λcu, λdl, λdu, λsl, λsu, νp, νs, νT = duals
    (; T, C, R, s0, sT, η) = d
    (; ch, ds, s) = aux_vars

    dL_s = λsu - λsl - νs + [νs[2:T]; 0]
    dL_s[T] -= νT[1]

    return [
        λcu - λcl - (1/sqrt(η))*νp + νs;
        λdu - λdl + sqrt(η)*νp - νs;
        dL_s
    ]
end

function get_device_comp_slack(d::Battery, p, duals, aux_vars)
    (; T, C, R, s0, sT, η) = d
    (; ch, ds, s) = aux_vars

    constr = [
        -ch
        ch .- R
        -ds
        ds .- R
        -s
        s .- C
        -p + sqrt(η)*ds - (1/sqrt(η))*ch
        -s - ds + ch + [s0; s[1:T-1]]
        -s[T] + sT
    ]

    return reduce(vcat, duals) .* constr
end

# Jacobian

function jacobian_kkt_z_device(d::Battery, p, duals, aux_vars)
    λcl, λcu, λdl, λdu, λsl, λsu, νp, νs, νT = duals
    (; T, C, R, s0, sT, η) = d
    (; ch, ds, s) = aux_vars

    k = get_dims(d)

    fx = [
        -ch
        ch .- R
        -ds
        ds .- R
        -s
        s .- C
    ]

    Dλ = Diagonal(reduce(vcat, duals[1:6])[:, 1])

    o = spzeros(T, T)

    HLx = spzeros(4T, 4T)
    Dfx = [
        o -I(T) o o;
        o I(T) o o;
        o o -I(T) o;
        o o I(T) o;
        o o o -I(T);
        o o o I(T)
    ]
    Dhx = [
        -I(T)   -(1/sqrt(η))*I(T)   sqrt(η)*I(T)    o
        o       -I(T)               I(T)            -I(T)+[zeros(1, T); I(T)[1:T-1, :]]
        spzeros(1, 4T - 1) -1
    ]

    n_ineq = length(fx)
    n_eq = 2T+1
    
    return [
        HLx         Dfx'            Dhx';
        Dλ*Dfx      Diagonal(fx)    spzeros(6T, 2T+1);
        Dhx         spzeros(2T+1, 6T) spzeros(2T+1, 2T+1)
    ]
end
