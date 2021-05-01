# CarbonNetworks.jl



## Data

To generate datasets, run `python make_files.py` in `src/`.

`branch_data.csv`
- source_node
- sink_node
- max_exchanges
- min_exchanges

`node_data.csv`
- id
- name
- demand_mef
- demand_mef_regression_score
- gen_max_1
- gen_max_2
- gen_max_3
- ...

`resource_data.csv`
- id
- name
- emission_factor

`case_X.csv`
- node_id
- demand
- gen1
- gen2
- gen3
- ...

TODO
- Make `gen_max` variable (if desired)
