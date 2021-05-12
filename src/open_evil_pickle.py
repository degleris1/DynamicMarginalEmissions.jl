# From https://stackoverflow.com/questions/45368255/error-in-loading-pickle
# https://stackoverflow.com/questions/2613800/how-to-convert-dos-windows-newline-crlf-to-unix-newline-lf-in-a-bash-script/19702943#19702943

import pickle
import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np

DATA_PATH = Path(os.getenv('CARBON_NETWORKS_DATA'))
PICKLE_PATH = DATA_PATH / 'SimpleDispatchData'
NERC_REGIONS = ['TRE', 'MRO', 'WECC', 'SPP', 'SERC', 'RFC', 'FRCC', 'NPCC']
YEAR = 2017

def cleaned(path):
    return path.parent / (path.name[:-4] + '_clean.obj')

def open_evil_pickle():
    gd_short = {}
    for region in NERC_REGIONS:
        try:
            path = PICKLE_PATH / f'generator_data_short_{region}_{YEAR}.obj'

            content = ''
            outsize = 0
            with open(path, 'rb') as infile:
                content = infile.read()

            with open(cleaned(path), 'wb') as output:
                for line in content.splitlines():
                    outsize += len(line) + 1
                    output.write(line + b'\n')

            # print(f"Done. Saved {(len(content)-outsize)} bytes for {region}.")

            with open(cleaned(path), 'rb') as f:
                gd_short[region] = pickle.load(f, encoding="latin1")

        except EOFError:
            print(f"Error -- No data for {region}")
    
    return gd_short

def compute_resource_cost(df, ba, resource):
    weekly_cost = pd.DataFrame(columns=np.arange(1, 53))
    for i in df.index:
        if df.loc[i, f'is_{resource}'] and df.loc[i, 'ba']==ba:
            for wk in np.arange(1, 53):
                weekly_cost.loc[i, wk] = df.loc[i, f'heat_rate{wk}']*df.loc[i, f'fuel_price{wk}']
        weekly_cost = weekly_cost.astype(float)
    return weekly_cost


def build_cost_df(gd_short):
    print("Building cost DF - this takes a bit of time.")
    resources = ['gas', 'coal', 'oil', 'nuclear', 'hydro', 'geothermal', 'biomass']
    df_cost = pd.DataFrame(columns = resources)
    for region in NERC_REGIONS:
        df = gd_short[region]['df']
        bas = df.ba.unique()
        for ba in bas:
            for resource in resources:
                weekly_cost = compute_resource_cost(df, ba, resource)
                if len(weekly_cost)>0:
                    df_cost.loc[ba, resource] = weekly_cost.median(axis=0).median()
    return df_cost
        
                