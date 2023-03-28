module DynamicMarginalEmissions

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
import PowerModels


# Exports
export open_datasets, parse_network_data, load_case, create_generation_map
export load_synthetic_network, load_demand_data, load_renewable_data

export PowerNetwork, PowerManagementProblem, solve!
export get_lmps, kkt_dims, flatten_variables, unflatten_variables, kkt

export DynamicPowerNetwork, DynamicPowerManagementProblem
export flatten_variables_dyn, unflatten_variables_dyn
export kkt_dyn, sensitivity_demand_dyn, kkt_dims_dyn, storage_kkt_dims

export sensitivity_price, sensitivity_demand
export compute_mefs
export plot_sensitivity_check
export make_dynamic, generate_random_data

export compute_jacobian_kkt_dyn, get_problem_dims
export make_pfdf_matrix

export generate_network

export INIT_COND, FINAL_COND

export evaluate, dispatch_and_compute_emissions_rates

# Files
include("model.jl")
include("dyn_model.jl")

include("sensitivity.jl")
include("dyn_sensitivity.jl")

include("parse_data.jl")
include("utils.jl")
include("simple_api.jl")

ECOS_OPT = () -> ECOS.Optimizer()
OPT = ECOS_OPT

end
