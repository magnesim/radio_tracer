import numpy as np
from datetime import datetime, timedelta
import opendrift 
#from opendrift.models.oceandrift import OceanDrift
#from opendrift.readers import reader_netCDF_CF_generic 
from release_descriptions import monthly_release, get_daily_weights, get_seed_date
import matplotlib.pyplot as plt 






# ###########################################################################
# Set the input file

#infn = '/home/magnes/projects/CERAD/RadioTracer/model_output/opendrift_tracers.nc'
#infn = '/home/magnes/projects/CERAD/RadioTracer/model_output/opendrift_tracers_2000-2003.nc'
infn = '/home/magnes/projects/CERAD/RadioTracer/model_output/opendrift_tracers_1993-1997_svim.nc'


# Select isotopes
isotops  = [
            '129I', 
            '236U',
            ]

# Size (m) of horisontal pixels for concentration 
psize   = 40000


# Filter out particles deeper than zmin
# this number will also be included in the file names 
zmin = 100
tag = '_{}m'.format(zmin)



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

# Create figure for release data
fig0=plt.figure(figsize=[10,7])
ax=plt.subplot()





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
    ax.plot(trajweights,label=isotop)
    ax.legend()
    ax.grid()
    ax.set_yscale('log')
    fn = '../plots/timeseries_releases.png'
    fig0.savefig(fname=fn)





    # #########################
    # Use opendrift function to compute horizontal histograms 
    h  = oa.get_histogram(pixelsize_m=psize, weights=trajweights)
    # Scale by pixel volume to get concentration (atoms/m3)
    h = h / (psize*psize*zmin)


    # Filter on time
    h = h.sel( time=slice(d0, d1) )






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
    # Sellafield releases
    b=h.isel(origin_marker=0).sum(dim='time')
    fn = '../plots/tracer_sellaf_{}{}.png'.format(isotop,tag)
    oa.plot(background=b.where(b>0), fast=True, show_elements=False, vmin=0, vmax=8.e14, clabel='Concentration '+isotop, filename=fn)
    # La Hague releases
    b=h.isel(origin_marker=1).sum(dim='time')
    fn = '../plots/tracer_lahague_{}{}.png'.format(isotop,tag)
    oa.plot(background=b.where(b>0), fast=True, show_elements=False, vmin=0, vmax=8.e14, clabel='Concentration '+isotop, filename=fn)
    # Total releases
    b=h.sum(dim='origin_marker').sum(dim='time')
    fn = '../plots/tracer_total_{}{}.png'.format(isotop,tag)
    oa.plot(background=b.where(b>0), fast=True, show_elements=False, vmin=0, vmax=8.e14, clabel='Concentration '+isotop, filename=fn)


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



# #################################################################
# Compute and plot ratio time series between the isotopes
if not len(h_save)==0:
    ratstr = isotops[0]+'/'+isotops[1]
    print('Compute isotope ratio '+ratstr)

    # Isotope 1
    r1 = h_save[0]
    # Isotope 2
    r2 = h_save[1]

    ii=0
    # Loop over the selected boxes (locations)
    # Compute time series for each location for each isotope for each source (origin_marker)
    # Plot time series, including total (sum of sources)
    for ibox in boxes:
        ii+=1
        fig=plt.figure(figsize=[9,9])
        ax1=plt.subplot(3,1,1)
        ax2=plt.subplot(3,1,2)
        ax3=plt.subplot(3,1,3)
        # Isotope 1
        t1 = r1.sel(lon_bin=slice(ibox['lon'][0], ibox['lon'][1]), lat_bin=slice(ibox['lat'][0], ibox['lat'][1])).sum(('lon_bin','lat_bin'))
        # Isotope 2
        t2 = r2.sel(lon_bin=slice(ibox['lon'][0], ibox['lon'][1]), lat_bin=slice(ibox['lat'][0], ibox['lat'][1])).sum(('lon_bin','lat_bin'))
        # Isotope ratio
        t3 = t1 / t2 

        # Isotope 1
        t1.isel(origin_marker=0).plot(label=isotops[0]+' SF',ax=ax1)
        t1.isel(origin_marker=1).plot(label=isotops[0]+' LH',ax=ax1)
        t1.sum(dim='origin_marker').plot(label=isotops[0]+' total',ax=ax1)
        # Isotope 2
        t2.isel(origin_marker=0).plot(label=isotops[1]+' SF',ax=ax2)
        t2.isel(origin_marker=1).plot(label=isotops[1]+' LH',ax=ax2)
        t2.sum(dim='origin_marker').plot(label=isotops[1]+' total',ax=ax2)
        # Isotope ratio
        t3.isel(origin_marker=0).plot(label=ratstr+' SF',ax=ax3)
        t3.isel(origin_marker=1).plot(label=ratstr+' LH',ax=ax3)
        t3.sum(dim='origin_marker').plot(label=ratstr+' total',ax=ax3)

        ax1.set_title(isotops[0]+' concentration')
        ax2.set_title(isotops[1]+' concentration')
        ax3.set_title('Isotope ratio '+ratstr)

        for ax in [ax1,ax2,ax3]:
            ax.legend()
            ax.grid()
        plt.suptitle(ibox['text'])
        fn = '../plots/timeseries_'+ibox['text']+'.png'
        plt.savefig(fn)
        



# #################################################################
# Compute and plot ratio maps between the isotopes
        
if not len(h_save)==0:
    ratstr = isotops[0]+'/'+isotops[1]
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

    fn = '../plots/r1.png'
    oa.plot(background=r1.where(r1>0), fast=True, show_elements=False, vmin=0, vmax=20000, clabel = 'r1 '+isotops[0], filename=fn)
    fn = '../plots/r2.png'
    oa.plot(background=r2.where(r2>0), fast=True, show_elements=False, vmin=0, vmax=6000, clabel = 'r2 '+isotops[1], filename=fn)

    fn = '../plots/ratio.png'
    oa.plot(background=ratio.where(ratio>-19), fast=True, show_elements=False, vmin=-2, vmax=3, clabel = 'log10 ratio '+ratstr, filename=fn)



# for om in [0, 1]:
#     background=h.isel(origin_marker=om)
#     fn = '../plots/tracer_density_{:02d}.mp4'.format(om)
#     oa.animation(background=background.where(background>0), bgalpha=1,
#                 corners=[-10.0, 30., 48., 80.], fast=False, show_elements=False, vmin=0, vmax=200, filename=fn)



