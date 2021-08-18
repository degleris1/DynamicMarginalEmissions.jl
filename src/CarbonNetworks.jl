module CarbonNetworks

# Module Imports
using CSV
using Convex
using DataFrames
using ECOS
using LightGraphs
using LinearAlgebra
using SimpleWeightedGraphs
using SparseArrays
using Zygote

using Base.Iterators: product
using PowerModels: parse_file, make_basic_network, make_per_unit!,
    calc_basic_incidence_matrix


# Exports
export open_datasets, parse_network_data, load_case, create_generation_map
export load_synthetic_network, load_demand_data, load_renewable_data

export PowerNetwork, PowerManagementProblem, solve!
export get_lmps, kkt_dims, flatten_variables, unflatten_variables

export DynamicPowerNetwork, DynamicPowerManagementProblem
export flatten_variables_dyn, unflatten_variables_dyn
export kkt_dyn, sensitivity_demand_dyn, kkt_dims_dyn, storage_kkt_dims

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
