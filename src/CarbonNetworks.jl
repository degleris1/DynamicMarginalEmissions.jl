module CarbonNetworks

using CSV
using DataFrames
using LightGraphs
using SimpleWeightedGraphs

DATA_PATH = ENV["CARBON_NETWORKS_DATA"]

include("./parse_data.jl")

export open_datasets, parse_network_data

end
