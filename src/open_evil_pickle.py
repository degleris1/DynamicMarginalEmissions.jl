# From https://stackoverflow.com/questions/45368255/error-in-loading-pickle
# https://stackoverflow.com/questions/2613800/how-to-convert-dos-windows-newline-crlf-to-unix-newline-lf-in-a-bash-script/19702943#19702943

import pickle
import sys
import os
from pathlib import Path

DATA_PATH = Path(os.getenv('CARBON_NETWORKS_DATA'))
PICKLE_PATH = DATA_PATH / 'SimpleDispatchData'
NERC_REGIONS = ['TRE', 'MRO', 'WECC', 'SPP', 'SERC', 'RFC', 'FRCC', 'NPCC']
YEAR = 2017

def cleaned(path):
    return path.parent / (path.name[:-4] + '_clean.obj')

def main():
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

            print(f"Done. Saved {(len(content)-outsize)} bytes for {region}.")

            with open(cleaned(path), 'rb') as f:
                gd_short[region] = pickle.load(f, encoding="latin1")

        except EOFError:
            print(f"Error -- No data for {region}")
    
    return gd_short