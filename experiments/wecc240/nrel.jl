
NREL_DIR = joinpath(config["data"]["DATA_DIR"], "nrel")
PATH_NETWORK = joinpath(NREL_DIR, "WECC240_2018_Generation_scheduling.xlsx")
PATH_GIS = joinpath(NREL_DIR, "Bus GIS.xlsx")
PATH_LOAD = joinpath(NREL_DIR, "Time_Series/TIMESERIES_DAUC_LOAD_V1.csv")
PATH_GEN = joinpath(NREL_DIR, "Time_Series/TIMESERIES_DAUC_GEN_V1.csv")

NREL_FUEL_MAP = Dict("Steam" => 'C', "Gas" => 'G')

function get_nrel_data(date)
    node = get_node_params()
    gen = get_gen_params(node.id_map)
    storage = get_storage_params(node.id_map)
    line = get_line_params(node.id_map)

    dt = get_dynamic_demand(date, node.id_map)
    gt = get_dynamic_gmax(date, gen)

    return (node=node, gen=gen, storage=storage, line=line, dt=dt, gt=gt)
end


function get_node_params()
    df = DataFrame(XLSX.readtable(PATH_NETWORK, "Bus")...)
    df_gis = DataFrame(XLSX.readtable(PATH_GIS, "Test1")...)
    
    busnames = df.busname
    regions = df.Region
    id_map = Dict(name => i for (i, name) in enumerate(busnames))

    nickname = df_gis[:, "Bus  Name"]
    lat = df_gis.Lat
    lon = df_gis.Long

    return (name=busnames, region=regions, id_map=id_map, nickname=nickname, lat=lat, lon=lon)
end

function get_gen_params(id_map)
    # How are we going to get heat rates?
    # One option: divide cost by price of fuel
    
    df = DataFrame(XLSX.readtable(PATH_NETWORK, "Generator")...)

    @assert all(df.Pmin .< df.MW1)
    @assert minimum(df.MW1) > 0
    
    gen_buses = Int[]
    gen_names = String[]
    gmin = Float64[]
    gmax = Float64[]
    ramp = Float64[]
    fuel = String[]
    cost = Float64[]

    # Handle generators with multiple costs
    for r in eachrow(df)
        b = r.busname
        append!(gen_buses, fill(typeof(b) <: Number ? b : parse(Int, b), 4))

        append!(gen_names, fill(r.genname, 4))
        
        append!(gmin, [r.Pmin, 0, 0, 0])

        capacities = [r.MW1, r.MW2, r.MW3, r.MW4]
        append!(gmax, [r.MW1; diff(capacities)])

        append!(ramp, fill(r.Ramp_Rate, 4))

        append!(fuel, fill(r.Gen_Type, 4))

        append!(cost, [r.Cost1, r.Cost2, r.Cost3, r.Cost4])
    end

    params = (bus=gen_buses, name=gen_names, gmin=gmin, gmax=gmax, ramp=ramp, fuel=fuel, cost=cost)

    # Remove generators without any capacity
    valid_gens = gmax .> 0
    params = map(p -> p[valid_gens], params)

    # Create node-generator map
    B = zeros(length(id_map), length(params.bus))
    for (k, b) in enumerate(params.bus)
        B[id_map[b], k] = 1
    end

    return (params..., B=B)
end

function get_storage_params(id_map)
    df = DataFrame(XLSX.readtable(PATH_NETWORK, "ESS")...)

    bus = df.busname
    C = df.Capacity .* (df.SOC_Max - df.SOC_Min)
    P = df.CH_Max
    ηc = df.ETA_ch
    ηd = df.ETA_dis

    S = zeros(length(id_map), length(bus))
    for (k, b) in enumerate(bus)
        S[id_map[b], k] = 1
    end
    
    return (bus=bus, C=C, P=P, ηc=ηc, ηd=ηd, S=S)
end

function get_line_params(id_map)
    df = DataFrame(XLSX.readtable(PATH_NETWORK, "Line")...)

    src = df.StartBusName
    snk = df.EndBusName
    β = 1 ./ df.X
    fmax = df.FlowLim

    m = length(src)
    A = zeros(length(id_map), m)
    for j in 1:m
        A[id_map[src[j]], j] = -1
        A[id_map[snk[j]], j] = 1
    end

    return (src=src, snk=snk, β=β, fmax=fmax, A=A)
end

function get_dynamic_demand(date, id_map)
    @assert year(date) == 2018

    # Find appropriate row
    df_load = DataFrame(CSV.File(PATH_LOAD))
    d = dayofyear(date)
    row = (d-1)*24 + hour(date) + 1
    #_f = r -> DateTime(2018, r.Month, r.Day, r.Period-1) == date
    #row = findfirst(_f, eachrow(df_load))

    # Create demand vector
    d = zeros(length(id_map))
    for (name, ind) in id_map
        d[ind] = df_load[row, string(name)]
    end

    return d
end

function get_dynamic_gmax(date, gen_params)
    @assert year(date) == 2018
    gmax = deepcopy(gen_params.gmax)

    # Find appropriate row
    df_gen = DataFrame(CSV.File(PATH_GEN))
    d = dayofyear(date)
    row = (d-1)*24 + hour(date) + 1
    #_f = r -> DateTime(2018, r.Month, r.Day, r.Period-1) == date
    #row = findfirst(_f, eachrow(df_gen))

    # Update max capacities
    gens = names(df_gen[row, 5:end])
    for g in gens
        @assert sum(gen_params.name .== g) == 1
    
        gen_ind = findfirst(gen_params.name .== g)
        gmax[gen_ind] = df_gen[row, g]
    end

    return gmax
end