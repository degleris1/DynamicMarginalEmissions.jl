import numpy as np
import pandas as pd
import os
from sklearn.linear_model import LinearRegression
from open_evil_pickle import open_evil_pickle, build_cost_df

DATA_PATH = os.getenv('CARBON_NETWORKS_DATA')

CASES_DATES = [
    '2020-01-30 12:00:00+0000',  # lowest sun/wnd generation in 2020
    '2020-12-23 19:00:00+0000',  # highest sun/wnd generation in 2020
    '2020-09-11 08:00:00+0000',  # lowest sun/wnd generation in summer 2020
    '2020-09-26 21:00:00+0000',  # highest sun/wnd generation in summer 2020
    '2020-09-20 10:00:00+0000',  # lowest total demand for 2020
    '2020-08-24 22:00:00+0000',  # highest total demand for 2020
]


def main():
    # load
    df_co2, df_elec = load_datasets()

    # build datasets
    df_branch = build_branch_data(df_elec)
    df_node = build_node_data(df_co2, df_elec)
    df_emissions = build_resource_data(df_elec)
    df_cost = build_cost_df(open_evil_pickle())

    # write datasets
    df_branch.to_csv(os.path.join(DATA_PATH, 'branch_data.csv'))
    df_node.to_csv(os.path.join(DATA_PATH, 'node_data.csv'))
    df_emissions.to_csv(os.path.join(DATA_PATH, 'resource_data.csv'))
    df_cost.to_csv(os.path.join(DATA_PATH, 'monetary_costs.csv'))

    # build cases:
    for i, date in enumerate(CASES_DATES):
        df_case = build_case(df_elec, df_co2, date)
        df_case.to_csv(os.path.join(DATA_PATH, f'case_{i}.csv'))
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
    columns_node = ['id', 'name', 'demand_mef', 'demand_mef_r2', 'net_demand_mef', 'net_demand_mef_r2',
                    'generation_mef', 'generation_mef_r2', 'net_generation_mef', 'net_generation_mef_r2']
    columns_res = [res+'_max' for res in resources]

    df_node = pd.DataFrame(columns=columns_node+columns_res)
    df_node = df_node.set_index('id')

    BAs = extract_BAs(df_co2)
    for res in resources:
        for ba in BAs:
            col = f'EBA.{ba}-ALL.NG.{res}.H'
            if col in df_elec.columns:
                df_node.loc[ba, f'{res}_max'] = df_elec[col].max()

    # compute all mefs
    for ba in BAs:
        for which in ['generation', 'net_generation', 'demand', 'net_demand']:
            (ba_, ba_co2), _ = extract_cols(ba, df_elec, df_co2, which=which)
            _, mef, r2 = compute_mef(ba_, ba_co2)
            df_node.loc[ba, f'{which}_mef'] = mef
            df_node.loc[ba, f'{which}_mef_r2'] = r2

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


def compute_mef(ba_, ba_co2):

    X, y = ba_.values.reshape(-1, 1), ba_co2.values
    reg = LinearRegression(fit_intercept=True).fit(X, y)

    return reg.predict(X), reg.coef_[0], reg.score(X, y)


def extract_cols(ba, df_elec, df_co2, which='generation'):

    D_elec_col = f'EBA.{ba}-ALL.D.H'
    D_co2_col = f'CO2_{ba}_D'
    # TI_col = f'EBA.{ba}-ALL.TI.H'
    NG_elec_col = f'EBA.{ba}-ALL.NG.H'
    NG_co2_col = f'CO2_{ba}_NG'
    WND_col = f'EBA.{ba}-ALL.NG.WND.H'
    SUN_col = f'EBA.{ba}-ALL.NG.SUN.H'

    if which == 'generation':
        ba_ = df_elec[NG_elec_col]
        ba_co2 = df_co2[NG_co2_col]
        xlabel = 'Delta in Generation'
    elif which == 'net_generation':
        ba_ = df_elec[NG_elec_col]
        ba_co2 = df_co2[NG_co2_col]
        if WND_col in df_elec.columns:
            ba_ = ba_-df_elec[WND_col]
        if SUN_col in df_elec.columns:
            ba_ = ba_-df_elec[SUN_col]
        xlabel = 'Delta in Net Generation'
    elif which == 'demand':
        ba_ = df_elec[D_elec_col]
        ba_co2 = df_co2[D_co2_col]
        xlabel = 'Delta in Demand'
    elif which == 'net_demand':
        ba_ = df_elec[D_elec_col]
        ba_co2 = df_co2[D_co2_col]
        if WND_col in df_elec.columns:
            ba_ = ba_-df_elec[WND_col]
        if SUN_col in df_elec.columns:
            ba_ = ba_-df_elec[SUN_col]
        xlabel = 'Delta in Net Generation'

    idx = ba_.index.intersection(ba_co2.index)

    ba_ = ba_.loc[idx]
    ba_co2 = ba_co2.loc[idx]
    ba_ = ba_.diff()
    ba_co2 = ba_co2.diff()

    ba_ = ba_.loc[~ba_.isna()]
    ba_co2 = ba_co2.loc[~ba_co2.isna()]

    return (ba_, ba_co2), xlabel


def build_case(df_elec, df_co2, date, save=False):

    resources = extract_resources(df_elec)
    BAs = extract_BAs(df_co2)
    df_case = pd.DataFrame(
        columns=['demand', 'net_demand'] + resources.tolist(), index=BAs)

    df_crt = df_elec.loc[date].copy()

    for ba in BAs:

        D_elec_col = f'EBA.{ba}-ALL.D.H'
        WND_col = f'EBA.{ba}-ALL.NG.WND.H'
        SUN_col = f'EBA.{ba}-ALL.NG.SUN.H'

        df_case.loc[ba, 'demand'] = df_crt[D_elec_col]
        ND = df_crt[D_elec_col]
        if WND_col in df_elec.columns:
            ND -= df_crt[WND_col]
        if SUN_col in df_elec.columns:
            ND -= df_crt[SUN_col]

        df_case.loc[ba, 'net_demand'] = ND

        for res in resources:
            res_col = f'EBA.{ba}-ALL.NG.{res}.H'
            if res_col in df_elec.columns:
                df_case.loc[ba, res] = df_crt[res_col]

    if save:
        df_case.to_csv()

    return df_case


if __name__ == "__main__":
    main()
