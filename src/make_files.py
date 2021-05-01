import numpy as np
import pandas as pd
import os

DATA_PATH = os.getenv('CARBON_NETWORKS_DATA')


def main():
    # load
    df_co2, df_elec = load_datasets()

    # build datasets
    df_branch = build_branch_data(df_elec)
    df_node = build_node_data(df_co2, df_elec)
    df_emissions = build_resource_data(df_elec)

    # write datasets
    df_branch.to_csv(os.path.join(DATA_PATH, 'branch_data.csv'))
    df_node.to_csv(os.path.join(DATA_PATH, 'node_data.csv'))
    df_emissions.to_csv(os.path.join(DATA_PATH, 'resource_data.csv'))
    return


def load_datasets():
    fnm_co2 = os.path.join(DATA_PATH, 'EBA_co2.csv')
    fnm_elec = os.path.join(DATA_PATH, 'EBA_elec.csv')
    df_co2 = pd.read_csv(fnm_co2, index_col=0, parse_dates=True)
    df_elec = pd.read_csv(fnm_elec, index_col=0, parse_dates=True)
    return df_co2, df_elec


def extract_BAs(df_co2):
    # extracting the names of the BAs present in the dataset
    nms = []
    for c in df_co2.columns:
        nms.append(c.split('_')[1])
    BAs = []
    for nm in nms:
        if '-' in nm:
            pass
        else:
            BAs.append(nm)
    BAs = np.unique(BAs)
    return BAs

def extract_resources(df_elec):
    # resources
    resources = []
    for c in df_elec.columns:
        if '-ALL.NG.' in c:
            resource = c.split('-ALL.NG.')[-1]
            if resource == 'H':
                continue
            else:
                resources.append(resource.split('.')[0])

    resources = np.unique(resources)
    return resources


def build_branch_data(df_elec):
    df = pd.DataFrame(
        columns=['source_node', 'sink_node', 'max_exchanges', 'min_exchanges'])
    sinks = []
    sources = []
    min_exchanges = []
    max_exchanges = []
    for c in df_elec.columns:
        if 'ID' in c:
            source_ba = c.split('.')[1].split('-')[0]
            sink_ba = c.split('.')[1].split('-')[1].split('.')[0]
            sources.append(source_ba)
            sinks.append(sink_ba)
            min_, max_ = df_elec[c].abs().min(), df_elec[c].abs().max()
            min_exchanges.append(min_)
            max_exchanges.append(max_)

    df.source_node = sources
    df.sink_node = sinks
    df.max_exchanges = max_exchanges
    df.min_exchanges = min_exchanges
    return df


def build_node_data(df_co2, df_elec):
    
    resources = extract_resources(df_elec)
    columns_node = ['id', 'name', 'demand_mef', 'demand_mef_regression_score']
    columns_res = [res+'_max' for res in resources]

    df_node = pd.DataFrame(columns=columns_node+columns_res)
    df_node = df_node.set_index('id')

    BAs = extract_BAs(df_co2)
    for res in resources:
        for ba in BAs:
            col = f'EBA.{ba}-ALL.NG.{res}.H'
            if col in df_elec.columns:
                df_node.loc[ba, f'{res}_max'] = df_elec[col].max()

    return df_node


def build_resource_data(df_elec):
    # UNK is 2017 average US power grid intensity according to Schivley 2018
    # unit is kg / MWh
    EMISSIONS_FACTORS = {
        "CO2": {
            "WAT": 4,
            "NUC": 16,
            "SUN": 46,
            "NG": 469,
            "WND": 12,
            "COL": 1000,
            "OIL": 840,
            "OTH": 439,
            "UNK": 439,
            "BIO": 230,
            "GEO": 42,
        }
    }

    resources = extract_resources(df_elec)
    df_emissions = pd.DataFrame(index=resources, columns=[
                                'name', 'emission_factor'])
    for res in resources:
        df_emissions.loc[res,
                         'emission_factor'] = EMISSIONS_FACTORS['CO2'][res]

    return df_emissions

if __name__ == "__main__":
    main()