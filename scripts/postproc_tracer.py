import numpy as np
from datetime import datetime, timedelta
import opendrift 
#from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic 
from release_descriptions import releases_sellafield, releases_lahague
import matplotlib.pyplot as plt 





#o = OceanDrift(loglevel=20, seed=0)
#


infn = '/home/magnes/projects/CERAD/RadioTracer/model_output/opendrift_tracers_2000-2003.nc'
#infn = '/home/magnes/projects/CERAD/RadioTracer/model_output/opendrift_tracers.nc'

d0 = datetime(2000,1,1)
d1 = datetime(2003,12,31)


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



isotops  = [
            '129I', 
            '236U',
            ]

# if False:
#     o.io_import_file(infn)

#     fn = '../plots/tracer_particles.mp4'
#     o.animation(fast=False, #clabel='Ocean current [m/s]',
#                 #color='z',
#     #            background=['x_sea_water_velocity', 'y_sea_water_velocity'],
#                 filename=fn,
#                 )




oa = opendrift.open_xarray(infn)
print(oa.ds.trajectory, len(oa.ds.trajectory))
ntra = len(oa.ds.trajectory)


h_save = []

# figure for release data
fig0=plt.figure(figsize=[10,7])
ax=plt.subplot()

for isotop in isotops:
    boxes=boxes_org.copy()
    [rel_sf_t, rel_sf_y] = releases_sellafield(isotop=isotop, dates=[d0,d1], number=ntra/2)
    [rel_lh_t, rel_lh_y] = releases_lahague(isotop=isotop, dates=[d0,d1], number=ntra/2)
    release = np.concatenate((rel_sf_y, rel_lh_y))

    # plot release data
    ax.plot(release,label=isotop)
    ax.legend()
    fn = '../plots/timeseries_releases.png'
    fig0.savefig(fname=fn)

    h  = oa.get_histogram(pixelsize_m=40000, weights=release)

#    print (h)
    h=h.sel( time=slice(d0, d1) )
#    print(h)
    
#    rw = h.sum(dim='origin_marker')


    # Plot time series in boxes
    if len(boxes)>3:
        figsize=[9,15]
    else:
        figsize=[9,10]
    fig1=plt.figure(figsize=figsize)
    ii=0
    for ibox in boxes:
        ii+=1
        print(ii,ibox)
        ax1 = plt.subplot(len(boxes),1,ii)
        t1 = h.sel(lon_bin=slice(ibox['lon'][0], ibox['lon'][1]), lat_bin=slice(ibox['lat'][0], ibox['lat'][1]))
        t1 = t1.sum(('lon_bin', 'lat_bin'))
        t1.isel(origin_marker=0).plot(label=isotop+' SF', ax=ax1)
        t1.isel(origin_marker=1).plot(label=isotop+' LH', ax=ax1)
        t1.sum(dim='origin_marker').plot(label='Total '+isotop, linestyle='--', ax=ax1)

        ax1.legend()
        ax1.set_title(ibox['text'])

    fn1 = '../plots/timeseries_{}.png'.format(isotop)
    fig1.savefig(fname=fn1)


# """
# text = [{'s': o.origin_marker[0], 'x': 8.55, 'y': 58.56, 'fontsize': 20, 'color': 'g',
#          'backgroundcolor': 'white', 'bbox': dict(facecolor='white', alpha=0.8), 'zorder': 1000},
#         {'s': o.origin_marker[1], 'x': 8.35, 'y': 58.42, 'fontsize': 20, 'color': 'g',
#          'backgroundcolor': 'white', 'bbox': dict(facecolor='white', alpha=0.8), 'zorder': 1000},
#         {'s': '* Station', 'x': station_lon, 'y': station_lat, 'fontsize': 20, 'color': 'k',
#          'backgroundcolor': 'white', 'bbox': dict(facecolor='none', edgecolor='none', alpha=0.4), 'zorder': 1000}]
# """


    b=h.isel(origin_marker=0).sum(dim='time')
    fn = '../plots/tracer_sellaf_{}.png'.format(isotop)
    oa.plot(background=b.where(b>0), fast=True, show_elements=False, vmin=0, vmax=8.e5, clabel='Concentration '+isotop, filename=fn)

    b=h.isel(origin_marker=1).sum(dim='time')
    fn = '../plots/tracer_lahague_{}.png'.format(isotop)
    oa.plot(background=b.where(b>0), fast=True, show_elements=False, vmin=0, vmax=8.e5, clabel='Concentration '+isotop, filename=fn)

    b=h.sum(dim='origin_marker').sum(dim='time')
    fn = '../plots/tracer_total_{}.png'.format(isotop)
    oa.plot(background=b.where(b>0), fast=True, show_elements=False, vmin=0, vmax=8.e5, clabel='Concentration '+isotop, filename=fn)


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



