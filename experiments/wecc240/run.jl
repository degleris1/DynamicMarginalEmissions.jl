include("util.jl")

SETUP_FILE = (length(ARGS) >= 1) ? ARGS[1] : nothing
RUN_NUMBER = (length(ARGS) >= 2) ? parse(Int, ARGS[2]) : "all"

function run(setup_file::String, run_number)
    
    # Constants
    HOURS_PER_DAY = 24
    RESULTS_SUBFOLDER = "results240"
    DATE_FORMAT = "yy-mm-dd-HH"




    # Load run configuration
    setup = TOML.parsefile(setup_file)

    run_id = setup["RUN_ID"]
    flag = setup["FLAG"]

    start_date = DateTime(setup["START_DATETIME"], DATE_FORMAT)
    num_days = setup["NUM_DAYS"]
    T = setup["TIME_HORIZON"]
    added_storage = get(setup, "ADDED_STORAGE", 0.0)




    # Decide whether or not to use static or dynamic runs
    load_case = (T == 1) ? make_static_case : (d -> make_dynamic_case(d, T))
    run_model = (T == 1) ? formulate_and_solve_static : (d -> formulate_and_solve_dynamic(d, T; added_storage=added_storage))

    # Create dates
    dates = start_date .+ Hour.(0 : T : (HOURS_PER_DAY*num_days-1))




    # Verify inputs
    @assert flag in ["", "clean"]
    if flag == "clean"
        @assert run == "all"
    end
    @assert 1 <= num_days <= 364  # NREL is missing 7 hours of data, so we can run at most 364 days
    @assert 1 <= T

    @assert all(d -> year(d) in [2004, 2018], dates)
    @assert all(d -> dayofyear(d) <= 364, dates)
    @assert (run_number == "all") || (0 <= run_number <= length(dates))




    # Set up directory
    results_path = joinpath(SAVE_DIR, RESULTS_SUBFOLDER, run_id)
    rp = f -> joinpath(results_path, f)

    if (run_number == "all") || (run_number == 0)
        # Check if we are overwriting data
        if ispath(results_path)
            @warn "Path $results_path already exists: you may be overwriting data."

            # If clean flag is passed, remove the directory and start from scratch
            if flag == "clean"
                rm(results_path, recursive=true)
            end
        end

        # Create directory if it doesn't exist already
        mkpath(results_path)

        # Put config file in directory for posterity
        cp(setup_file, rp("setup.toml"), force=true)

        # Save a single case for later analyses
        case = load_case(first(dates))
        bson(rp("case.bson"); case...)
    end



    # Now actually run something!
    run_and_save = d -> bson(rp(Dates.format(d, DATE_FORMAT) * ".bson"); run_model(d)...)

    if run_number == "all"  # Run all cases
        run_and_save.(dates)
        return
    elseif run_number == 0  # Just initialize the directory
        return
    else  # run_number > 0
        run_and_save(dates[run_number])
    end
end

run(SETUP_FILE, RUN_NUMBER)
