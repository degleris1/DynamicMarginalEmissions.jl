module CarbonNetworks

# Module Imports
using Convex
using ECOS
using LinearAlgebra

using SparseArrays: sparse, spzeros, blockdiag
using Zygote: forward_jacobian, gradient

# Imports for redefinition
import Convex: solve!, evaluate
import Base: reshape

OPT = ECOS.Optimizer

# using CSV
# using DataFrames
# using LightGraphs
# using SimpleWeightedGraphs
# using SparseArrays
# using Zygote

# using Base.Iterators: product
# using PowerModels: parse_file, make_basic_network, make_per_unit!,
#     calc_basic_incidence_matrix


# Exports - Devices
export get_cost
export Generator, StaticGenerator, Demand
export RampGenerator, with_ramping
export Battery

# Exports - Power Management Problem
export DynamicPowerManagementProblem
export solve!, kkt
export flatten_pmp_result, unflatten_pmp_result
export get_num_lines, get_num_devices, get_time_horizon, get_dims

# Exports - Network Utilities
export make_pfdf_matrix, make_device_pfdf_matrix

# Exports - Emissions Utilities
export get_total_emissions, get_lmes

# Exports - Backwards Compatibility
export PowerNetwork, PowerManagementProblem, compute_mefs


# export load_synthetic_network, load_demand_data, load_renewable_data
# export DynamicPowerNetwork, DynamicPowerManagementProblem
# export generate_random_data
# export compute_jacobian_kkt_dyn, get_problem_dims
# export generate_network


# Files
include("devices/template.jl")
include("devices/generator.jl")
include("devices/ramping_generator.jl")
include("devices/storage.jl")

include("pmp.jl")
include("jacobian.jl")

include("network_utils.jl")
include("emissions_utils.jl")
include("backwards_compat.jl")

#include("parse_data.jl")
#include("utils.jl")
end
