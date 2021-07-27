module CarbonNetworks

# Module Imports
using CSV
using Convex
using DataFrames
using ECOS
using LightGraphs
using LinearAlgebra
using SimpleWeightedGraphs
using Zygote
using RecipesBase

using Base.Iterators: product
using SparseArrays

# Exports
export open_datasets, parse_network_data, load_case, create_generation_map

export PowerNetwork, PowerManagementProblem, solve!
export get_lmps, kkt_dims, flatten_variables, unflatten_variables

export DynamicPowerNetwork, DynamicPowerManagementProblem
export flatten_variables_dyn, unflatten_variables_dyn
export kkt_dyn, sensitivity_demand_dyn, kkt_dims_dyn

export sensitivity_price, sensitivity_demand
export loss_and_grad, stochastic_loss_and_grad
export compute_mefs
export plot_sensitivity_check

# Files
include("model.jl")
include("dyn_model.jl")

include("parse_data.jl")
include("sensitivity.jl")
include("dyn_sensitivity.jl")
include("descent.jl")

OPT = () -> ECOS.Optimizer(verbose=false)

end
