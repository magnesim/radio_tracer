import numpy as np
import matplotlib.pyplot as plt 


def plot_scatter_obsmodel(obs, mod, isotope, folder, box ):

    # observation dates
    odate = obs['Date'].to_numpy()
    if len(obs['Date'])==0:
        print('No obs available, ')
        return

    # observation values
    if isotope=='ratio':
        oval = obs['129I Concentration (at/L)'].astype(float) / obs['236U Concentration (at/L)'].astype(float)
    else:
        oval = obs['{} Concentration (at/L)'.format(isotope)] .astype(float)

    # model dates
    mdate = np.array([item.values for item in mod['time'][:]]) 
    # index of corrsponding dates 
    tdate = np.array([ np.argmin( np.abs(mdate - item) ) for item in odate ])

    # model values at obs dates
    modSF = mod.isel(time=tdate, origin_marker=0)
    modLH = mod.isel(time=tdate, origin_marker=1)
    modT  = mod.isel(time=tdate).sum(dim='origin_marker')#.astype(float)

    # plotting
    figs, axs =plt.subplots()
    maxcal = np.maximum( np.max(modT ), np.max(oval) ) * 1.1 
    axs.plot(oval, modSF, marker='s', markeredgecolor='k', markerfacecolor='none', ls='none', label='Sellfield')
    axs.plot(oval, modLH, marker='*', markeredgecolor='k', markerfacecolor='none', ls='none', label='LaHague')
    axs.plot(oval, modT,  marker='^', markeredgecolor='k', markerfacecolor='none', ls='none', label='Total')
    axs.plot([0,maxcal], [0,maxcal], ls='--', lw=1.4, alpha=.6, c='k')
    axs.set_xlim([0,maxcal])
    axs.set_ylim([0,maxcal])
    axs.set_xlabel('obs')
    axs.set_ylabel('model')
    axs.legend()
    axs.grid()
    axs.set_title(box+' '+isotope)
    figs.savefig('{}/scatter_{}_{}.png'.format(folder, box, isotope ) )
    plt.close()

    return

