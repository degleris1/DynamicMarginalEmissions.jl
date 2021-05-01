import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import numpy as np

def plot_mef(ba, df_elec, df_co2, which='generation'):

    (ba_, ba_co2), xlabel = extract_cols(ba, df_elec, df_co2, which=which)    

    #compute MEF
    preds, mef, r2 = compute_mef(ba_, ba_co2)

    #plot
    _, ax = plt.subplots()
    plt.scatter(ba_, ba_co2, s=1)
    plt.plot(ba_.values, preds, color='darkorange')
    plt.xlabel(xlabel)
    plt.ylabel('Delta in Emissions')
    plt.title(f'R2: {np.around(r2, 2)} - MEF: {np.around(mef, 2)}')
    plt.grid()

    return (ba_, ba_co2)

def compute_mef(ba_, ba_co2):

    X, y = ba_.values.reshape(-1, 1), ba_co2.values
    reg = LinearRegression(fit_intercept=True).fit(X, y)

    return reg.predict(X), reg.coef_[0], reg.score(X,y)

def extract_cols(ba, df_elec, df_co2, which='generation'):

    D_elec_col = f'EBA.{ba}-ALL.D.H'
    D_co2_col = f'CO2_{ba}_D'
    # TI_col = f'EBA.{ba}-ALL.TI.H'
    NG_elec_col = f'EBA.{ba}-ALL.NG.H'
    NG_co2_col = f'CO2_{ba}_NG'
    WND_col = f'EBA.{ba}-ALL.NG.WND.H'
    SUN_col = f'EBA.{ba}-ALL.NG.SUN.H'

    if which=='generation':
        ba_ = df_elec[NG_elec_col]
        ba_co2 = df_co2[NG_co2_col]
        xlabel='Delta in Generation'
    elif which=='net_generation':
        ba_ = df_elec[NG_elec_col]
        ba_co2 = df_co2[NG_co2_col]
        if WND_col in df_elec.columns:
            ba_ = ba_-df_elec[WND_col]
        if SUN_col in df_elec.columns:
            ba_ = ba_-df_elec[SUN_col]
        xlabel='Delta in Net Generation'
    elif which=='demand':
        ba_ = df_elec[D_elec_col]
        ba_co2 = df_co2[D_co2_col]
        xlabel = 'Delta in Demand'
    elif which=='net_demand':
        ba_ = df_elec[D_elec_col]
        ba_co2 = df_co2[D_co2_col]
        if WND_col in df_elec.columns:
            ba_ = ba_-df_elec[WND_col]
        if SUN_col in df_elec.columns:
            ba_ = ba_-df_elec[SUN_col]
        xlabel='Delta in Net Generation'

    idx = ba_.index.intersection(ba_co2.index)

    ba_ = ba_.loc[idx]
    ba_co2 = ba_co2.loc[idx]
    ba_ = ba_.diff()
    ba_co2 = ba_co2.diff()

    ba_ = ba_.loc[~ba_.isna()]
    ba_co2 = ba_co2.loc[~ba_co2.isna()]

    return (ba_, ba_co2), xlabel

