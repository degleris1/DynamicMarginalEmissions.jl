module CarbonNetworks

# Module Imports
using Convex
using ECOS

using SparseArrays: sparse, spzeros, blockdiag

import Convex: solve!


# using CSV
# using DataFrames
# using LightGraphs
# using LinearAlgebra
# using SimpleWeightedGraphs
# using SparseArrays
# using Zygote

# using Base.Iterators: product
# using PowerModels: parse_file, make_basic_network, make_per_unit!,
#     calc_basic_incidence_matrix


# Exports - Devices
export get_cost, get_constraints
export Generator, StaticGenerator, Demand

# Exports - Power Management Problem
export DynamicPowerManagementProblem
export solve!, kkt
export flatten_pmp_result, unflatten_pmp_result
export get_num_lines, get_num_devices, get_time_horizon, get_dims

# Exports - Utilities
export make_pfdf_matrix, make_device_pfdf_matrix

# export open_datasets, parse_network_data, load_case, create_generation_map
# export load_synthetic_network, load_demand_data, load_renewable_data

# export get_lmps

# export DynamicPowerNetwork, DynamicPowerManagementProblem
# export flatten_variables_dyn, unflatten_variables_dyn
# export kkt_dyn, sensitivity_demand_dyn, kkt_dims_dyn, storage_kkt_dims

# export sensitivity_price, sensitivity_demand
# export compute_mefs
# export generate_random_data

# export compute_jacobian_kkt_dyn, get_problem_dims
# export make_pfdf_matrix

# export generate_network


# Files
include("devices/template.jl")
include("devices/generator.jl")

include("pmp.jl")
include("jacobian.jl")

include("network_utils.jl")


#include("parse_data.jl")
#include("utils.jl")
end
