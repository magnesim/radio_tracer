import numpy as np
import matplotlib.pyplot as plt 


def plot_vertdistr(oa, boxes, vint):
    import pandas as pd

    timefromfile  = np.array([pd.to_datetime(item).to_pydatetime() for item in oa.ds['time'].values])
    timedifffromfile = (timefromfile[-1] - timefromfile[0])
    nt = len(timefromfile)



    for box in boxes:
        print('Vertical profile: ', box['text'], box['lon'], box['lat'])
        min_lon=box['lon'][0]
        max_lon=box['lon'][1]
        min_lat=box['lat'][0]
        max_lat=box['lat'][1]
        mask_lon = (oa.ds.lon >= min_lon) & (oa.ds.lon <= max_lon)
        mask_lat = (oa.ds.lat >= min_lat) & (oa.ds.lat <= max_lat)

        r1 = oa.ds.where(mask_lon & mask_lat)
        z = r1.z


        fig=plt.figure()
        ax=plt.subplot()

        for ii in range(0, nt, nt//vint):
            tt = timefromfile[ii]
            print(ii, tt)
            zu=z.sel(time=tt )
            hist, bin_edges = np.histogram(zu,bins=np.arange(-100,10,4))
            bin_centr = (bin_edges[1:]+bin_edges[:-1]) /2
            ax.plot(hist, bin_centr, label=tt)
        ax.legend()
        ax.grid()
        fn = '../plots/vertical_distribution_{}.png'.format(box['text'].replace(' ',''))
        plt.savefig(fn, dpi=200, bbox_inches='tight')
        plt.close()

    return







def plot_scatter_obsmodel(obs, mod, isotope, folder, box , printtoscreen=False, ratiosum=False):
    import pandas as pd

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
    if isotope=='ratio':
        modT = ratiosum.isel(time=tdate)
    else:
        modT  = mod.isel(time=tdate).sum(dim='origin_marker')#.astype(float)


    obslines= []
    if printtoscreen:
        for ii in range(len(tdate)):
            oline =  '{}; {}; {}; {:.3E}; {:.3E}; {:.3E}; {:.3E}'.format(box, pd.to_datetime(odate[ii]).strftime("%Y-%m-%d"), isotope, oval.values[ii], modT.values[ii], modSF.values[ii], modLH.values[ii])
            obslines.append(oline)

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

    return obslines

