using Distributions
using LightGraphs
using Random

using Convex, ECOS

ECOS_QUIET = () -> ECOS.Optimizer(verbose=false)

function make_synthetic_network(n)
    # Make graph
    G = watts_strogatz(n, 4, 0.2)

    # Convert to incidence matrix
    A = incidence_matrix(G, oriented=true)
    m = size(A, 2)

    # Generate costs
    f = rand(Exponential(5), n) .+ 2

    # Generate generation and flow capacities
    gmax = rand(Gamma(5.0, 5.0), n)
    pmax = rand(Gamma(1.0, 1.0), m)

    B = I;
    
    return f, pmax, gmax, A, B
end

function make_cases(num_cases, f, pmax, gmax, A)
    cases = []
    for _ in 1:num_cases
        d = rand(Uniform(0.0, 1.05), n) .* gmax

        opf = PowerManagementProblem(f, d, pmax, gmax, A)
        solve!(opf, ECOS_QUIET)

        # println(opf.problem.status)
        if Int(opf.problem.status) == 2
            continue
        end

        push!(cases, (d=d, g=evaluate(opf.g)))
    end
    
    return cases
end