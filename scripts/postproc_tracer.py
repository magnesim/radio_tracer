import numpy as np
from datetime import datetime, timedelta
import opendrift 
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic 






o = OceanDrift(loglevel=20, seed=0)
#


infn = '/home/magnes/projects/CERAD/RadioTracer/model_output/opendrift_tracers.nc'



o.io_import_file(infn)

fn = '../plots/tracer_particles.mp4'
o.animation(fast=False, #clabel='Ocean current [m/s]',
            #color='z',
#            background=['x_sea_water_velocity', 'y_sea_water_velocity'],
            filename=fn,
            )




oa = opendrift.open_xarray(infn)

h = oa.get_histogram(pixelsize_m=10000)

b=h.isel(origin_marker=0).sum(dim='time')
fn = '../plots/tracer_sellaf.png'
oa.plot(background=b.where(b>0), fast=True, show_elements=False, vmin=0, vmax=1000, clabel='First seeding', filename=fn)

b=h.isel(origin_marker=1).sum(dim='time')
fn = '../plots/tracer_lahague.png'
oa.plot(background=b.where(b>0), fast=True, show_elements=False, vmin=0, vmax=1000, clabel='First seeding', filename=fn)

for om in [0, 1]:
    background=h.isel(origin_marker=om)
    fn = '../plots/tracer_density_{:02d}.mp4'.format(om)
    oa.animation(background=background.where(background>0), bgalpha=1,
                corners=[-10.0, 30., 48., 80.], fast=False, show_elements=False, vmin=0, vmax=200, filename=fn)
