include("util.jl")

dim1, dim2 = size(PARETO_CONFIG)
for (i, j) in Iterators.product(1:dim1, 1:dim2)
    println(i, j)
    run_expansion_planning(PARETO_CONFIG[i, j], "pareto$(i)_$(j)")
end