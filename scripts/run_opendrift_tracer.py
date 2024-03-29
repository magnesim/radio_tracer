import numpy as np
from datetime import datetime, timedelta
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import  reader_ROMS_native, reader_netCDF_CF_generic 

import copernicus_marine_client as copernicusmarine
#from opendrift.readers.reader_netCDF_CF_generic import Reader



# Approx positions Sellafield
#lon_SF = -3.710
#lat_SF = 54.357
# A bit further away from shore
lon_SF = -3.9232
lat_SF = 54.1750

# Approx position LaHague
#lon_LH = -1.918
#lat_LH = 49.656
# A bit further away from shore
lon_LH = -2.05
lat_LH = 49.744



# Input parameters
# 
outdir = '/home/magnes/projects/CERAD/RadioTracer/model_output'
time_start = datetime(1993,1,1,0)
time_end   = datetime(1997,1,1,0)
#time_end   = datetime(1993,7,10,0)

simlen_h   = (time_end - time_start).total_seconds() / 3600  # total time in hours
time_step  = 3600*4     # in seconds
ntraj      = 26000      # total number of trajectories







# Initialize OpenDrift object and add reader data

o = OceanDrift(loglevel=20, seed=0)
#reader_cmems = reader_netCDF_CF_generic.Reader('https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_my_0.083_P1D-m')








# #################
# R E A D E R S
# ocean current data


#o.add_readers_from_list([
    #'https://nrt.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i',
                         #'https://nrt.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_anfc_merged-uv_PT1H-i',
                       #'https://nrt.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_anfc_merged-uv_PT1H-i'  # funker
                        #'https://nrt.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_anfc_0.083deg_PT1H-m'    # funker
#                        'https://thredds.met.no/thredds/dodsC/nansen-legacy-ocean/svim_daily_agg'
#                         ])

#o.add_reader([reader_cmems])

r_svim  = reader_ROMS_native.Reader('https://thredds.met.no/thredds/dodsC/nansen-legacy-ocean/svim_daily_agg')

ds_cmems   = copernicusmarine.open_dataset(
#    dataset_id='cmems_mod_glo_phy_anfc_merged-uv_PT1H-i',
            dataset_id='cmems_mod_glo_phy_my_0.083deg_P1D-m',
            username='<your cmems username>', password='<your cmems password>'
            )

ds_cmems_int   = copernicusmarine.open_dataset(
#    dataset_id='cmems_mod_glo_phy_anfc_merged-uv_PT1H-i',
            dataset_id='cmems_mod_glo_phy_myint_0.083deg_P1D-m',
            username='<your cmems username>', password='<your cmems password>'
            )

r_cmems = reader_netCDF_CF_generic.Reader(ds_cmems)
r_cmems_int = reader_netCDF_CF_generic.Reader(ds_cmems_int)




# Add the selected readers in prioritized order
o.add_reader([
              r_svim, 
              r_cmems,
              #r_cmems_int,
              ])


# Adjusting some configuration
o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:diffusivitymodel','constant') 
o.set_config('environment:fallback:ocean_vertical_diffusivity', 1.e-6)

# By default, radionuclides do not strand towards coastline
o.set_config('general:coastline_action', 'previous')
# Vertical mixing requires fast time step
o.set_config('vertical_mixing:timestep', 900.) # seconds
o.set_config('drift:horizontal_diffusivity', 30.)

o.set_config('drift:deactivate_north_of', 80.1)
o.set_config('drift:deactivate_east_of', 26.)

o.list_configspec()







# Seed the elements


iniz=np.random.rand(ntraj//2) * -10. # seeding the radionuclides in the upper 10m


o.seed_elements(lon=lon_SF, lat=lat_SF, number=int(ntraj/2), radius=1000, 
                time = [time_start, time_end],
                z=iniz,
                #origin_marker=0
                )
o.seed_elements(lon=lon_LH, lat=lat_LH, number=int(ntraj/2), radius=1000, 
                time = [time_start, time_end],
                z=iniz,
                #origin_marker=1
                )







# Run the model

fn = outdir + '/opendrift_tracers.nc'
o.run( steps=simlen_h*3600/time_step+1,  time_step=time_step,  
      time_step_output=3600*24, outfile=fn )








# Take care of the output

print('make animation..')
fn = outdir + '/animation_tracers.mp4'
o.animation(fast=False, #clabel='Ocean current [m/s]',
            #color='z',
#            background=['x_sea_water_velocity', 'y_sea_water_velocity'],
            filename=fn,
            )




