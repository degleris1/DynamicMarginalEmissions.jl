import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import numpy as np

from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Times New Roman'], 'size':10})
rc('text', usetex=True)

FIGSIZE = (3.5, 2.5)

props = dict( facecolor='white', alpha=1)



def plot_mef(ba, df_elec, df_co2, which='generation'):

    (ba_, ba_co2), xlabel = extract_cols(ba, df_elec, df_co2, which=which)    

    #compute MEF
    preds, mef, r2 = compute_mef(ba_, ba_co2)

    #plot
    _, ax = plt.subplots(figsize=FIGSIZE, tight_layout=True)
    ax.scatter(ba_, ba_co2, s=1)
    ax.plot(ba_.values, preds, color='darkorange')
    ax.set_xlabel(r'$\Delta$'+f'{xlabel} [MWh/h]')
    ax.set_ylabel(r'$\Delta\ CO_2$ [kg/h]')
    # place a text box in upper left in axes coords
    x_text = ba_.quantile(.0)
    y_text = ba_co2.quantile(.99)
    textstr=rf'$R^2$: {np.around(r2, 2)}' +'\n' +rf'MEF: {np.around(mef, 2)} kg/MWh'
    ax.text(x_text, y_text, textstr, bbox=props, fontsize=8)
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
        xlabel='Generation'
    elif which=='net_generation':
        ba_ = df_elec[NG_elec_col]
        ba_co2 = df_co2[NG_co2_col]
        if WND_col in df_elec.columns:
            ba_ = ba_-df_elec[WND_col]
        if SUN_col in df_elec.columns:
            ba_ = ba_-df_elec[SUN_col]
        xlabel='Net Generation'
    elif which=='demand':
        ba_ = df_elec[D_elec_col]
        ba_co2 = df_co2[D_co2_col]
        xlabel = 'Demand'
    elif which=='net_demand':
        ba_ = df_elec[D_elec_col]
        ba_co2 = df_co2[D_co2_col]
        if WND_col in df_elec.columns:
            ba_ = ba_-df_elec[WND_col]
        if SUN_col in df_elec.columns:
            ba_ = ba_-df_elec[SUN_col]
        xlabel='Net Demand'

    idx = ba_.index.intersection(ba_co2.index)

    ba_ = ba_.loc[idx]
    ba_co2 = ba_co2.loc[idx]
    ba_ = ba_.diff()
    ba_co2 = ba_co2.diff()

    ba_ = ba_.loc[~ba_.isna()]
    ba_co2 = ba_co2.loc[~ba_co2.isna()]

    return (ba_, ba_co2), xlabel

