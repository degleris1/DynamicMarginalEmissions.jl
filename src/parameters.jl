# From deChalendar - Tracking emissions codebase
# UNK is 2017 average US power grid intensity according to Schivley 2018
# unit is kg / MWh
EMISSIONS_FACTORS = Dict(
		"WAT" => 4,
        "NUC" => 16,
        "SUN" => 46,
        "NG" => 469,
        "WND" => 12,
        "COL" => 1000,
        "OIL" => 840,
        "OTH" => 439,
        "UNK" => 439,
        "BIO" => 230,
        "GEO" => 42,
	)