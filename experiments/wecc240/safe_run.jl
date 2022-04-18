using TOML

# Find script directory
println("Script executed from: $(pwd())")
println("Script location: $(@__DIR__)")

# Get script and argument names
command = joinpath(@__DIR__, "run.jl")
setup = ARGS[1]

# Initialize
run(`julia $command $setup 0`)

# Run all the days separately
num_days = TOML.parsefile(setup)["NUM_DAYS"]
horizon = TOML.parsefile(setup)["TIME_HORIZON"]

num_runs = (num_days * 24) ÷ horizon
@show num_runs

for i in 1:num_runs
    run(`julia $command $setup $i`)
end
