using CarbonNetworks
using Random

# script to test the dynamic networks and the computation of KKT conditions
# to be erased when the code works
Random.seed!(2)
n = 10
l = 100

# Make graph
G = watts_strogatz(n, 4, 0.2)

# Convert to incidence matrix
A = incidence_matrix(G, oriented=true)
m = size(A, 2)

# Create generator-node mapping
node_map = vcat(collect(1:n), rand(1:n, l-n))
B = spzeros(n, l)
for (gen, node) in enumerate(node_map)
    B[node, gen] = 1
end

# Generate carbon costs
cq = zeros(l)
cl = rand(Exponential(2), l)

# Generate demands
d = rand(Bernoulli(0.8), n) .* rand(Gamma(3.0, 3.0), n)

# Generate generation and flow capacities
gmax = rand(Gamma(4.0, 3.0), l) + (B'd)  # This is just to make the problem easier
pmax = rand(Gamma(1.0, 0.1), m)