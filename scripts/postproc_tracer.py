import numpy as np
from datetime import datetime, timedelta
import opendrift 
#from opendrift.models.oceandrift import OceanDrift
#from opendrift.readers import reader_netCDF_CF_generic 
from release_descriptions import monthly_release, get_daily_weights, get_seed_date
import matplotlib.pyplot as plt 

import cartopy.crs as ccrs
import pandas as pd 




# ###########################################################################
# Set the input file

#infn = '/home/magnes/projects/CERAD/RadioTracer/model_output/opendrift_tracers.nc'
#infn = '/home/magnes/projects/CERAD/RadioTracer/model_output/opendrift_tracers_2000-2003.nc'
#infn = '/home/magnes/projects/CERAD/RadioTracer/model_output/opendrift_tracers_1993-1997_svim.nc'
infn = '/home/magnes/projects/CERAD/RadioTracer/model_output/opendrift_tracers_1993-1997_cmems.nc'


# Select isotopes
isotops  = [
            '129I', 
            '236U',
            ]

isotops_fmt = {'129I': r'$^{129}$I',
               '236U': r'$^{236}$U'}

print('{} {}'.format( isotops_fmt[isotops[0]], isotops_fmt[isotops[1]] ))

# Size (m) of horisontal pixels for concentration 
psize   = 40000


# Filter out particles deeper than zmin
# this number will also be included in the file names 
zmin = 100
tag = '_{}m'.format(zmin)
tag = tag+''   # here, you can add specific name tag that will appear in all file names




# ######################################################################
# O B S E R V A T I O N S
obs_compare = False
obs_compare = True
obsfile = '../data/obs_template.csv'


# ###########################################################################
# Define/edit boxes for analysis (time series etc) here

box1_lon = [5, 7]        # Denmark west coast
box1_lat = [53, 55]
box2_lon = [7, 10]       # Norwegian coast
box2_lat = [64, 66]
box3_lon = [8.8, 11.1]   # Skagerrak
box3_lat = [57.5, 59.1]
box4_lon = [1.3, 3.3]    # North Sea
box4_lat = [53.9, 56.9] 
box5_lon = [16.2, 19.2]  # Tromsø
box5_lat = [70.1, 72.0]


boxes_org = [
    {'lon': box1_lon, 'lat': box1_lat, 'text': 'Dan W Coast', 'fc': 'none', 'alpha': 0.8, 'lw': 1, 'ec': 'k'},
    {'lon': box2_lon, 'lat': box2_lat, 'text': 'Norw Coast',  'fc': 'none', 'alpha': 0.8, 'lw': 1, 'ec': 'k'},
    {'lon': box3_lon, 'lat': box3_lat, 'text': 'Skagerrak',   'fc': 'none', 'alpha': 0.8, 'lw': 1, 'ec': 'k'},
    {'lon': box4_lon, 'lat': box4_lat, 'text': 'North Sea',   'fc': 'none', 'alpha': 0.8, 'lw': 1, 'ec': 'k'},
    {'lon': box5_lon, 'lat': box5_lat, 'text': 'Tromsø',      'fc': 'none', 'alpha': 0.8, 'lw': 1, 'ec': 'k'},
        ]


# if False:
#     o.io_import_file(infn)

#     fn = '../plots/tracer_particles.mp4'
#     o.animation(fast=False, #clabel='Ocean current [m/s]',
#                 #color='z',
#     #            background=['x_sea_water_velocity', 'y_sea_water_velocity'],
#                 filename=fn,
#                 )




# plot the boxes on a map
fig = plt.figure()
ax=plt.subplot(projection=ccrs.Orthographic(5,68))
proj_pp=ccrs.PlateCarree()
for ibox in boxes_org:
    [lon1, lon2] = ibox['lon']
    [lat1, lat2] = ibox['lat']
    box_coords = [(lon1,lat1), (lon2,lat1), (lon2,lat2), (lon1,lat2), (lon1,lat1)]
    ax.plot(*zip(*box_coords), transform=proj_pp, color='red', linewidth=2 )
    ax.text(lon2, lat2, ibox['text'], horizontalalignment = 'right', transform=proj_pp, zorder=5)
ax.set_extent([-10.0, 20, 50, 80], proj_pp)
ax.coastlines()
ax.gridlines()
ax.stock_img()

fn='../plots/map_locations{}.png'.format(tag)
plt.savefig(fn)



# ###########################################################################
# Read Opendrift output file as xarray

oa = opendrift.open_xarray(infn)
#print(oa.ds.trajectory, len(oa.ds.trajectory))
#ntra = len(oa.ds.trajectory)



# Mask on depth 
# deeper than zmin will be masked out
oa.ds = oa.ds.where(oa.ds.z > -zmin)



# ###########################################################################
# Get the time period (datesfromfile), number of trajectories (ntraj)
# and their corresponding seeding date (seed_dates) from opendrft output file (infn)
# Note: First all Sellafield trajectories, then all LaHague trajectories

[seed_dates, datesfromfile, ntraj] = get_seed_date(infn)
[d0,d1] = datesfromfile

# Compute number of days and the number of released particles each day
# Assume equal number for each location 
ndays = int((d1 - d0).total_seconds() / 86400 + 1) 
print(f'ndays: {ndays}, first day: {d0}, last day: {d1}',ndays)
date_arr = [d0 + timedelta(days=item) for item in range(ndays)] 
ntrajperday = np.zeros(ndays)
for ii,dd in enumerate(date_arr):
    ntrajperday[ii] = np.sum(seed_dates[:ntraj//2]==dd)






# ###############################################################################
# Store the histograms 
h_save = []
hage_save = []

# Create figure for release data
fig0=plt.figure(figsize=[10,7])
ax0=plt.subplot()





# ##############################################################################
# Run through loop over the isotopes, 
# Plot figures and store necessary data 

for isotop in isotops:
    boxes=boxes_org.copy()

    # Extract the mothly releases from SF and LH
    [rel_sf_t, rel_sf_y] = monthly_release(isotop, dates=[d0,d1], location='Sellafield') 
    [rel_lh_t, rel_lh_y] = monthly_release(isotop, dates=[d0,d1], location='LaHague') 



    # Get the daily weighted release from SF and LH
    [date_arrSF, weightsSF, total_atomsSF] = get_daily_weights( release=[rel_sf_t, rel_sf_y], dates=[d0,d1] )
    [date_arrLH, weightsLH, total_atomsLH] = get_daily_weights( release=[rel_lh_t, rel_lh_y], dates=[d0,d1] )


    # Get the daily weights,
    # scaled by number of trajectories released each day in opendrift simulation
    weightsSF = weightsSF / ntrajperday
    weightsLH = weightsLH / ntrajperday

    # Get the weights per trajectory,
    # scaled by the total release from SF and LH
    trajweightsSF = np.zeros(ntraj//2)
    trajweightsLH = np.zeros(ntraj//2)
    for ii in range(ntraj//2):
        dd = seed_dates[ii]
        idx = date_arr.index(dd)
        trajweightsSF[ii] = weightsSF[idx] * total_atomsSF 
        trajweightsLH[ii] = weightsLH[idx] * total_atomsLH

    # Concatenate the releases 
    # and trajectory weights
    # trajweights is a list with length of the total number of trajectories,
    # containing the number of atoms corresponding to each trajectory
    release = np.concatenate((rel_sf_y, rel_lh_y))
    dates = np.concatenate((rel_sf_t, rel_lh_t))
    trajweights = np.concatenate((trajweightsSF, trajweightsLH))





    # #######################
    # plot number of atoms per trajectory, 
    # first SF trajectories, the LH trajectories
    ax0.plot(trajweights,label=isotop)
    ax0.legend()
    ax0.set_yscale('log')
    ax0.grid()
    fn = '../plots/releases_ts{}.png'.format(tag)
    fig0.savefig(fname=fn)





    # #########################
    # Use opendrift function to compute horizontal histograms 
    h  = oa.get_histogram(pixelsize_m=psize, weights=trajweights)
    # Scale by pixel volume to get concentration (atoms/m3)
    h = h / (psize*psize*zmin)
    h = h / 1000.  # Convert to atoms/L


    # Filter on time
    h = h.sel( time=slice(d0, d1) )
#    h = h.sel( time=slice(d1-timedelta(days=365), d1) )




# """
# text = [{'s': o.origin_marker[0], 'x': 8.55, 'y': 58.56, 'fontsize': 20, 'color': 'g',
#          'backgroundcolor': 'white', 'bbox': dict(facecolor='white', alpha=0.8), 'zorder': 1000},
#         {'s': o.origin_marker[1], 'x': 8.35, 'y': 58.42, 'fontsize': 20, 'color': 'g',
#          'backgroundcolor': 'white', 'bbox': dict(facecolor='white', alpha=0.8), 'zorder': 1000},
#         {'s': '* Station', 'x': station_lon, 'y': station_lat, 'fontsize': 20, 'color': 'k',
#          'backgroundcolor': 'white', 'bbox': dict(facecolor='none', edgecolor='none', alpha=0.4), 'zorder': 1000}]
# """




    # ######################
    # Plot maps of tracer concentration integrated over time
    vminconc=5
    vmaxconc=13
    map_proj = ccrs.Orthographic(-25, 65)

    sources = ['Sellafield','LaHague','Total']
    for ii,b in enumerate([h.isel(origin_marker=0).sum(dim='time'), 
                           h.isel(origin_marker=1).sum(dim='time'), 
                           h.sum(dim='origin_marker').sum(dim='time') ]):
        fig=plt.figure(figsize=[12,7])
        ax = plt.subplot(projection=map_proj)
        b=np.log10(b)
        LONS, LATS = np.meshgrid(b['lon_bin'], b['lat_bin'])
        m1 = ax.pcolormesh(LONS,LATS, b.transpose(), vmin=vminconc, vmax=vmaxconc,  cmap='plasma' , shading='nearest', transform=proj_pp)
        ax.coastlines()
        ax.gridlines()
        cb=plt.colorbar(m1, label='log10 {} Concentration (at/L)'.format(isotops_fmt[isotop]))
        ax.set_title(isotops_fmt[isotop]+' '+sources[ii])
        fn = '../plots/tracer_{}_{}{}.png'.format(sources[ii], isotop,tag)
        fig.savefig(fn)
        #oa.plot(transform=proj_pp, projection=map_proj, background=b.where(b>0), fast=False, show_elements=False, vmin=vminconc, vmax=vmaxconc, clabel='log10 Concentration '+isotops_fmt[isotop], filename=fn, subplot_kws={'projection': map_proj, 'transform':proj_pp})

    # # La Hague releases
    # b=h.isel(origin_marker=1).sum(dim='time')
    # b=np.log10(b)
    # fn = '../plots/tracer_lahague_{}{}.png'.format(isotop,tag)
    # oa.plot(background=b.where(b>0), fast=True, show_elements=False, vmin=vminconc, vmax=vmaxconc, clabel='log10 Concentration '+isotops_fmt[isotop], filename=fn)
    # # Total releases
    # b=h.sum(dim='origin_marker').sum(dim='time')
    # b=np.log10(b)
    # fn = '../plots/tracer_total_{}{}.png'.format(isotop,tag)
    # oa.plot(background=b.where(b>0), fast=True, show_elements=False, vmin=vminconc, vmax=vmaxconc, clabel='log10 Concentration '+isotops_fmt[isotop], filename=fn)



    if isotop == isotops[0]:
        # ######################
        # Plot maps of tracer age integrated over time
        # This is similar for both radionuclides, and 
        # is only nesecary to do for the first radionuclide
        hage  = oa.get_histogram(pixelsize_m=psize, weights=oa.ds['age_seconds'], density=False)
        num   = oa.get_histogram(pixelsize_m=psize, weights=None, density=False)
        hage = hage / (86400*356)    # years
        hage = hage / num
#        hage = hage.sel( time=slice(d1-timedelta(days=365), d1))
        maxage=5

        for ii,b in enumerate([hage.isel(origin_marker=0).mean(dim='time'), 
                            hage.isel(origin_marker=1).mean(dim='time'), 
                            hage.mean(dim='origin_marker').mean(dim='time') ]):
            fig=plt.figure(figsize=[12,7])
            ax = plt.subplot(projection=map_proj)
            LONS, LATS = np.meshgrid(b['lon_bin'], b['lat_bin'])
            m1 = ax.pcolormesh(LONS,LATS, b.transpose(), vmin=0, vmax=maxage, cmap='rainbow', shading='nearest', transform=proj_pp)
            ax.coastlines()
            ax.gridlines()
            cb=plt.colorbar(m1, label='Age {} (years)'.format(isotops_fmt[isotop]))
            ax.set_title(isotops_fmt[isotop]+' '+sources[ii])
            fn = '../plots/tracerage_{}_{}{}.png'.format(sources[ii], isotop,tag)
            fig.savefig(fn)




        # # Sellafield releases
        # b=hage.isel(origin_marker=0).mean(dim='time')
        # fn = '../plots/tracerage_sellaf_{}{}.png'.format(isotop,tag)
        # oa.plot(background=b.where(b>0), fast=True, show_elements=False, vmin=0, vmax=maxage, clabel='Age '+isotops_fmt[isotop], filename=fn)
        # # La Hague releases
        # b=hage.isel(origin_marker=1).mean(dim='time')
        # fn = '../plots/tracerage_lahague_{}{}.png'.format(isotop,tag)
        # oa.plot(background=b.where(b>0), fast=True, show_elements=False, vmin=0, vmax=maxage, clabel='Age '+isotops_fmt[isotop], filename=fn)
        # # Total releases
        # b=hage.mean(dim='time').mean(dim='origin_marker')
        # fn = '../plots/tracerage_total_{}{}.png'.format(isotop,tag)
        # oa.plot(background=b.where(b>0), fast=True, show_elements=False, vmin=0, vmax=maxage, clabel='Age '+isotops_fmt[isotop], filename=fn)




    if False:
        # Make animations
        bbox= boxes.copy()
        print('animation',bbox)

        rw = h.sum(dim='origin_marker')

        fn = '../plots/tracer_density_{}.mp4'.format(isotop)
        oa.animation(background=rw.where(rw>0), bgalpha=1, fast=False,
                show_elements=False, vmin=0, vmax=2.e5, filename=fn) #, box=bbox, text=text)

    if len(isotops)==2:
        h_save.append(h)
        hage_save.append(hage)



# #################################################################
# Compute and plot ratio time series between the isotopes
if not len(h_save)==0:
    ratstr = isotops_fmt[isotops[0]]+'/'+isotops_fmt[isotops[1]]
    print('Compute isotope ratio '+ratstr)

    # Isotope 1
    r1 = h_save[0]
    # Isotope 2
    r2 = h_save[1]

    # Loop over the selected boxes (locations)
    # Compute time series for each location for each isotope for each source (origin_marker)
    # Plot time series, including total (sum of sources)
    for ibox in boxes:
        fig=plt.figure(figsize=[9,12])
        ax1=plt.subplot(4,1,1)
        ax2=plt.subplot(4,1,2)
        ax3=plt.subplot(4,1,3)
        ax4=plt.subplot(4,1,4)
        # Isotope 1
        t1 = r1.sel(lon_bin=slice(ibox['lon'][0], ibox['lon'][1]), lat_bin=slice(ibox['lat'][0], ibox['lat'][1])).sum(('lon_bin','lat_bin'))
        # Isotope 2
        t2 = r2.sel(lon_bin=slice(ibox['lon'][0], ibox['lon'][1]), lat_bin=slice(ibox['lat'][0], ibox['lat'][1])).sum(('lon_bin','lat_bin'))
        # Isotope ratio
        t3 = t1 / t2 

        # Isotope 1
        t1.isel(origin_marker=0).plot(label=isotops_fmt[isotops[0]]+' SF',ax=ax1)
        t1.isel(origin_marker=1).plot(label=isotops_fmt[isotops[0]]+' LH',ax=ax1)
        t1.sum(dim='origin_marker').plot(label=isotops_fmt[isotops[0]]+' total',ax=ax1)
        # Isotope 2
        t2.isel(origin_marker=0).plot(label=isotops_fmt[isotops[1]]+' SF',ax=ax2)
        t2.isel(origin_marker=1).plot(label=isotops_fmt[isotops[1]]+' LH',ax=ax2)
        t2.sum(dim='origin_marker').plot(label=isotops_fmt[isotops[1]]+' total',ax=ax2)
        # Isotope ratio
        t3.isel(origin_marker=0).plot(label=ratstr+' SF',ax=ax3)
        t3.isel(origin_marker=1).plot(label=ratstr+' LH',ax=ax3)
        t3.sum(dim='origin_marker').plot(label=ratstr+' total',ax=ax3)
        # Source contribution 129I
        t4SF = t1.isel(origin_marker=0) / t1.sum(dim='origin_marker')*100.
        t4LH = t1.isel(origin_marker=1) / t1.sum(dim='origin_marker')*100.
        t4SF.plot(c='C0', label='SF contribution {}'.format(isotops_fmt[isotops[0]]),ax=ax4)
        t4LH.plot(c='C1', label='LH contribution {}'.format(isotops_fmt[isotops[0]]),ax=ax4)
        # Source contribution 236U
        t5SF = t2.isel(origin_marker=0) / t2.sum(dim='origin_marker')*100.
        t5LH = t2.isel(origin_marker=1) / t2.sum(dim='origin_marker')*100.
        t5SF.plot(ls=':', c='C0', label='SF contribution {}'.format(isotops_fmt[isotops[1]]),ax=ax4)
        t5LH.plot(ls=':', c='C1', label='LH contribution {}'.format(isotops_fmt[isotops[1]]),ax=ax4)

        ax1.set_title(isotops_fmt[isotops[0]]+' concentration')
        ax2.set_title(isotops_fmt[isotops[1]]+' concentration')
        ax3.set_title('Isotope ratio '+ratstr)
        ax4.set_title('Source contribution')




        if obs_compare:
            obsdata = pd.read_csv(obsfile, header=0)

            # Convert to datetime object
            obsdata['Date'] = np.array([datetime.strptime(str(item), '%Y%m%d') for item in obsdata["Date"] ])

            # Find the relevant observations, within the actual box
            obsdata = obsdata[((obsdata['Longitude']>=ibox['lon'][0]) & (obsdata['Longitude']<=ibox['lon'][1]) &
                               (obsdata['Latitude']>=ibox['lat'][0]) & (obsdata['Latitude']<=ibox['lat'][1]))]
            
            # ...and time
            obsdata = obsdata[ ((obsdata['Date']>=d0) & (obsdata['Date']<=d1)) ]

            # ... and depth
            obsdata = obsdata[ np.abs(obsdata['Depth'])<=np.abs(zmin) ]
            
            # Pick out only iodine / uranium observations, and compute ratio where both nuclides are present
            obs_iodine = obsdata[obsdata['129I Concentration (at/L)']!='Nan']
            obs_uran = obsdata[obsdata['236U Concentration (at/L)']!='Nan']
            obs_ratio = obsdata[ (obsdata['129I Concentration (at/L)']!='Nan') & (obsdata['236U Concentration (at/L)']!='Nan') ]
            ratio = obs_ratio['129I Concentration (at/L)'].astype(float) /  obs_ratio['236U Concentration (at/L)'].astype(float)


            # Plot with the time series
            print('{:22} {:6}; N obs: {}'.format(ibox['text'], isotops[0], len(obs_iodine['129I Concentration (at/L)'])))
            p1= ax1.plot(obs_iodine['Date'], obs_iodine['129I Concentration (at/L)'].astype(float), ls='none', marker='o', markeredgecolor='k', markerfacecolor='none', label='obs {}'.format(isotops_fmt[isotops[0]]))
            print('{:22} {:6}; N obs: {}'.format(ibox['text'], isotops[1], len(obs_uran['236U Concentration (at/L)'])))
            p2= ax2.plot(obs_uran['Date'], obs_uran['236U Concentration (at/L)'].astype(float), ls='none', marker='o', markeredgecolor='k', markerfacecolor='none', label='obs {}'.format(isotops_fmt[isotops[1]]))
            print('{:22} ratio ; N obs: {}'.format(ibox['text'], len(ratio)))
            p3= ax3.plot(obs_ratio['Date'], ratio, ls='none', marker='o', markeredgecolor='k', markerfacecolor='none', label='obs ratio')

        for ax in [ax1,ax2,ax3,ax4]:
            ax.legend()
            ax.grid()
            ax.set_xlim([d0,d1])
        ax1.set_ylabel('{} concentration (at/L)'.format(isotops_fmt[isotops[0]]))
        ax2.set_ylabel('{} concentration (at/L)'.format(isotops_fmt[isotops[1]]))
        ax3.set_ylabel('ratio')
        ax4.set_ylabel('percent')
        ax3.set_yscale('log')
        plt.suptitle(ibox['text'])
        fn = '../plots/location_ts_{}{}.png'.format(ibox['text'].replace(' ',''), tag)
        plt.savefig(fn)
        

        if obs_compare:
            # Compare obs with nearest model value
            from plotting_tools import plot_scatter_obsmodel
            plot_scatter_obsmodel(obs_iodine, t1, isotope='129I', folder='../plots', box=ibox['text'].replace(' ',''))
            plot_scatter_obsmodel(obs_uran, t2, isotope='236U', folder='../plots', box=ibox['text'].replace(' ',''))
            plot_scatter_obsmodel(obs_ratio, t3, isotope='ratio', folder='../plots', box=ibox['text'].replace(' ',''))








# #################################################################
# Compute and plot tracer age time series (and age ratios) between the isotopes
if not len(h_save)==0:
    ratstr = isotops_fmt[isotops[0]]+'/'+isotops_fmt[isotops[1]]
    print('Compute mean age '+ratstr)

    # Isotope 1
    r1 = hage_save[0]
    # Isotope 2
    r2 = hage_save[1]

    # Loop over the selected boxes (locations)
    # Compute time series for each location for each isotope for each source (origin_marker)
    # Plot time series, including total (sum of sources)
    for ibox in boxes:
        fig=plt.figure(figsize=[9,7])
        ax1=plt.subplot(1,1,1)
#        ax2=plt.subplot(2,1,2)
#        ax3=plt.subplot(3,1,3)
        # Isotope 1
        t1 = r1.sel(lon_bin=slice(ibox['lon'][0], ibox['lon'][1]), lat_bin=slice(ibox['lat'][0], ibox['lat'][1])).mean(('lon_bin','lat_bin'))
#        t1std = r1.sel(lon_bin=slice(ibox['lon'][0], ibox['lon'][1]), lat_bin=slice(ibox['lat'][0], ibox['lat'][1])).std(('lon_bin','lat_bin'))
#        t1std = t1std+t1
        # Isotope 2
        t2 = r2.sel(lon_bin=slice(ibox['lon'][0], ibox['lon'][1]), lat_bin=slice(ibox['lat'][0], ibox['lat'][1])).mean(('lon_bin','lat_bin'))
        t2std = r2.sel(lon_bin=slice(ibox['lon'][0], ibox['lon'][1]), lat_bin=slice(ibox['lat'][0], ibox['lat'][1])).std(('lon_bin','lat_bin'))
        # age ratio
#        t3 = t1 / t2 

        # Isotope 1
        t1.isel(origin_marker=0).plot(label=isotops_fmt[isotops[0]]+' SF',ax=ax1)
        #t1std.isel(origin_marker=0).plot(label=isotops[0]+' std SF',ax=ax1)
        t1.isel(origin_marker=1).plot(label=isotops_fmt[isotops[0]]+' LH',ax=ax1)
        t1.mean(dim='origin_marker').plot(label=isotops_fmt[isotops[0]]+' total',ax=ax1)
        # # Isotope 2
        # t2.isel(origin_marker=0).plot(label=isotops[1]+' SF',ax=ax2)
        # t2.isel(origin_marker=1).plot(label=isotops[1]+' LH',ax=ax2)
        # t2.mean(dim='origin_marker').plot(label=isotops[1]+' total',ax=ax2)
        # Isotope ratio
#        t3.isel(origin_marker=0).plot(label=ratstr+' SF',ax=ax3)
#        t3.isel(origin_marker=1).plot(label=ratstr+' LH',ax=ax3)
#        t3.mean(dim='origin_marker').plot(label=ratstr+' total',ax=ax3)

        ax1.set_title(isotops_fmt[isotops[0]]+' age')
#        ax2.set_title(isotops[1]+' age')
#        ax3.set_title('age ratio '+ratstr)

        for ax in [ax1]:
            ax.legend()
            ax.grid()
        plt.suptitle(ibox['text'])
        fn = '../plots/location_agets_{}{}.png'.format(ibox['text'].replace(' ',''), tag)
        plt.savefig(fn)
        






# #################################################################
# Compute and plot ratio maps between the isotopes
        
if not len(h_save)==0:
    ratstr = isotops_fmt[isotops[0]]+'/'+isotops_fmt[isotops[1]]
    print('Compute isotope ratio '+ratstr)

    r1 = h_save[0]
    r1 = r1.sum(dim='origin_marker')

    r2 = h_save[1]
    r2 = r2.sum(dim='origin_marker')

    r1=r1.isel(time=-1)
    r2=r2.isel(time=-1)

    ratio = r1 / r2 
#    ratio = ratio.isel(time=-1)
    ratio = np.log10(ratio)

    fn = '../plots/map_isotope1_{}{}.png'.format(isotops[0], tag)
    oa.plot(background=r1.where(r1>0), fast=True, show_elements=False, vmin=0, vmax=20000, clabel = 'r1 '+isotops[0], filename=fn)
    fn = '../plots/map_isotope2_{}{}.png'.format(isotops[1], tag)
    oa.plot(background=r2.where(r2>0), fast=True, show_elements=False, vmin=0, vmax=6000, clabel = 'r2 '+isotops[1], filename=fn)

    fn = '../plots/map_ratio{}.png'.format(tag)
    oa.plot(background=ratio.where(ratio>-19), fast=True, show_elements=False, vmin=-2, vmax=3, clabel = 'log10 ratio '+ratstr, filename=fn)



# for om in [0, 1]:
#     background=h.isel(origin_marker=om)
#     fn = '../plots/tracer_density_{:02d}.mp4'.format(om)
#     oa.animation(background=background.where(background>0), bgalpha=1,
#                 corners=[-10.0, 30., 48., 80.], fast=False, show_elements=False, vmin=0, vmax=200, filename=fn)



