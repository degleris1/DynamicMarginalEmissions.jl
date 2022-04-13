include("util.jl")

SETUP_FILE = (length(ARGS) > 0) ? ARGS[1] : nothing

function run(setup_file::String)
    # Constants
    HOURS_PER_DAY = 24
    RESULTS_SUBFOLDER = "results240"
    DATE_FORMAT = "yy-mm-dd-HH"

    # Load run configuration
    setup_file = joinpath(@__DIR__, setup_file)
    setup = TOML.parsefile(setup_file)

    run_id = setup["RUN_ID"]
    start_date = DateTime(setup["START_DATETIME"], DATE_FORMAT)
    num_days = setup["NUM_DAYS"]
    T = setup["TIME_HORIZON"]

    # Verify inputs
    @assert 1 <= num_days <= 364  # NREL is missing 7 hours of data, so we can run at most 364 days
    @assert 1 <= T

    # Decide whether or not to use static or dynamic runs
    load_case = (T == 1) ? make_static_case : (d -> make_dynamic_case(d, T))
    run_model = (T == 1) ? formulate_and_solve_static : (d -> formulate_and_solve_dynamic(d, T))

    # Create dates
    dates = start_date .+ Hour.(0 : T : (HOURS_PER_DAY*num_days-1))

    @assert all(d -> year(d) in [2004, 2018], dates)
    @assert all(d -> dayofyear(d) <= 364, dates)

    # Set up directory
    results_path = joinpath(SAVE_DIR, RESULTS_SUBFOLDER, run_id)
    if ispath(results_path)
        @warn "Path $results_path already exists: you may be overwriting data."
    end
    mkpath(results_path)
    rp = f -> joinpath(results_path, f)

    # Put config file in directory for posterity
    cp(setup_file, rp("setup.toml"), force=true)

    # Save a single case for later analyses
    case = load_case(first(dates))
    bson(rp("case.bson"); case...)

    for d in dates
        results_d = run_model(d)
        name_d = Dates.format(d, DATE_FORMAT) * ".bson"

        bson(rp(name_d); results_d...)
    end
end

run(SETUP_FILE)
