# Make code backwards compatible with previous (non-modular code)
mutable struct OldPMP
    pmp::DynamicPowerManagementProblem
    net
    d
    generator_devices
    demand_devices
    pmp_result
end

# When creating a power network, don't do anything: just wrap everything in a NamedTuple!
PowerNetwork(fq, fl, pmax, gmax, A, B, F) = (fq=fq, fl=fl, pmax=pmax, gmax=gmax, A=A, B=B, F=F)

function PowerManagementProblem(net::NamedTuple, d)
    # Create T
    T = 1

    # Create devices
    # We should have `k+n` generators
    k = size(net.B, 2)
    n = size(net.A, 1)

    # The first `k` devices will be the generators
    gen_nodes = [findfirst(==(1), net.B[:, j]) for j in 1:k]
    gen_devices = [StaticGenerator(0, net.gmax[j], net.fl[j], T) for j in 1:k]

    # The next `n` devices will be the loads
    curtailment_cost = maximum(net.fl) * 2
    load_nodes = collect(1:n)
    load_devices = [Demand([d[i]], curtailment_cost) for i in 1:n]

    nodes = append!(gen_nodes, load_nodes)
    devices = append!(gen_devices, load_devices)

    # Create device PFDF matrix
    F = net.F[:, nodes]

    # Create fmax
    m = length(net.pmax)
    fmax = net.pmax * ones(1, T)
    fmax = [fmax; fmax]

    # Wrap problem in struct
    pmp = DynamicPowerManagementProblem(devices, F, fmax, T)
    return OldPMP(pmp, net, d, collect(1:k), collect(k+1:k+n), nothing)
end

function Base.getproperty(obj::OldPMP, sym::Symbol)
    if sym == :g
        return obj.pmp_result.p[obj.generator_devices, :]
    elseif sym == :problem
        return obj.pmp_result.problem
    else
        return getfield(obj, sym)
    end
end

function solve!(old_pmp::OldPMP, solver)
    pmp_result = solve!(old_pmp.pmp, solver)
    old_pmp.pmp_result = pmp_result
    return nothing
end

function compute_mefs(old_pmp::OldPMP, net, d, co2_rates)
    (; pmp, pmp_result, generator_devices, demand_devices) = old_pmp
    return get_lmes(pmp, pmp_result, demand_devices, generator_devices, co2_rates)
end