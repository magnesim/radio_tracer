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







def plot_scatter_obsmodel(obs, mod, isotope, folder, box , printtoscreen=False, ratiosum=False, tag=''):
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
    figs.savefig('{}/scatter_{}_{}{}.png'.format(folder, box, isotope, tag ) )
    plt.close()

    return obslines






def plot_diffratio(data, isotop=None, sources=None, isofmt=None, tag='', vmaxd=None, vmaxr=None, datestr=''):
    import cartopy.crs as ccrs
    from cartopy import feature as cfeature
    
    map_proj = ccrs.Orthographic(2, 65)
    proj_pp=ccrs.PlateCarree()


    if vmaxd==None:
        vmaxdiff = 1e8
        vmindiff = -vmaxdiff
    else:
        vmindiff = -vmaxd
        vmaxdiff = vmaxd
    if vmaxr==None:
        vmaxrat = 4
        vminrat = -vmaxrat
    else:
        vminrat = -vmaxr
        vmaxrat = vmaxr

    print('plotting diff and ratio for ',isotop, isofmt[isotop])
    print('Diff: {} - {} '.format(sources[0], sources[1]))
    print('Ratio: {} / {} '.format(sources[0], sources[1]))

    fig=plt.figure(figsize=[16,7])
    ax1 = plt.subplot(1,2,1, projection=map_proj)
    ax2 = plt.subplot(1,2,2, projection=map_proj)
    #b=np.log10(b)
    LONS, LATS = np.meshgrid(data['lon_bin'], data['lat_bin'])

    # diff plot
    diff = data.isel(origin_marker=0) - data.isel(origin_marker=1)
    m1 = ax1.pcolormesh(LONS,LATS, diff.transpose(), vmin=vmindiff, vmax=vmaxdiff,  cmap='RdBu' , shading='nearest', transform=proj_pp, zorder=4)
    cb=plt.colorbar(m1, label='Diff {} Concentration (at/L)'.format(isofmt[isotop]))
    ax1.set_title(isofmt[isotop]+' Diff '+sources[0]+' - '+sources[1])

    # ratio plot
    ratio = np.log10( data.isel(origin_marker=0) / data.isel(origin_marker=1))
    m2 = ax2.pcolormesh(LONS,LATS, ratio.transpose(), vmin=vminrat, vmax=vmaxrat,  cmap='RdBu' , shading='nearest', transform=proj_pp, zorder=4)
    cb=plt.colorbar(m2, label='log10(Ratio) {} Concentration (at/L)'.format(isofmt[isotop]))
    ax2.set_title(isofmt[isotop]+' Ratio '+sources[0]+' / '+sources[1])
    plt.suptitle(datestr)
    for ax in [ax1,ax2]:
        ax.coastlines(zorder=6)
        ax.gridlines(zorder=7)
        ax.add_feature(cfeature.LAND, zorder=5)
    
    fn = '../plots/tracerdiffratio_{}{}.png'.format( isotop,tag)
    fig.savefig(fn, dpi=290)
    plt.close(fig)


    return



def running_mean_over_time(dataset, window):
    import xarray as xr
    """
    Compute a running mean over the time dimension of a given xarray dataset.

    Parameters:
    - dataset (xarray.Dataset or xarray.DataArray): Input dataset.
    - window (int): Size of the rolling window.

    Returns:
    - xarray.Dataset or xarray.DataArray: Dataset or DataArray with running mean applied.
    """
    # Check if dataset is a DataArray
    if isinstance(dataset, xr.DataArray):
        return dataset.rolling(time=window, min_periods=1, center=True).mean()
    # Check if dataset is a Dataset
    elif isinstance(dataset, xr.Dataset):
        return dataset.apply(lambda x: x.rolling(time=window, min_periods=1, center=True).mean(), keep_attrs=True)
    else:
        raise TypeError("Input must be an xarray.Dataset or xarray.DataArray")

# Example usage:
# Assuming 'ds' is your xarray dataset
# running_mean_ds = running_mean_over_time(ds, window=3)



def xarray_datasets_to_csv(datasets, output_file):
    import pandas as pd
    """
    Write a number of xarray datasets, all at the same format as columns, to a CSV file.

    Parameters:
    - datasets (dict): Dictionary of xarray datasets with keys as column names.
    - output_file (str): Output CSV file path.
    """
    # Convert xarray datasets to pandas DataFrames
    dataframes = {name: ds.to_dataframe() for name, ds in datasets.items()}
    
    # Merge DataFrames on the index (time dimension assumed)
    merged_dataframe = pd.concat(dataframes, axis=1)
    
    # Write merged DataFrame to CSV
    merged_dataframe.to_csv(output_file)

# Example usage:
# Assuming 'datasets' is a dictionary containing xarray datasets
# xarray_datasets_to_csv(datasets, "output.csv")
