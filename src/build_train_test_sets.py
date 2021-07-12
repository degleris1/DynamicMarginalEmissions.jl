from make_files import load_datasets, build_case
import argparse
import os
import numpy as np

DATA_PATH = os.getenv('CARBON_NETWORKS_DATA')


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", help="Number of timepoints to extract", type=int)
    args = parser.parse_args()
    n = args.n

    df_co2, df_elec = load_datasets()

    #DATES = np.random.choice(df_elec.index, n, replace=False)
    DATES = df_elec.index

    # create repository
    TRAIN_PATH = os.path.join(
        DATA_PATH, 'TRAIN'
    )
    os.makedirs(TRAIN_PATH, exist_ok=True)

    # build cases:
    for i, date in enumerate(DATES):
        nm = date.strftime("%Y%m%d_%H")
        df_case = build_case(df_elec, df_co2, date)
        df_case.to_csv(os.path.join(TRAIN_PATH, f'{nm}.csv'))

    return


if __name__ == "__main__":
    main()
