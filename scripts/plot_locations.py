
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs

saveplots=False
saveplots=True


outfolder  = '../plots'

lats  =  np.array([ 
62.7718, 62.7718333333333, 66.9385, 66.9385, 66.8345, 66.8345, 68.315, 68.315, 70.6392, 70.997
])


lons =  np.array([
4.3663, 4.3663, 1.73, 1.73, 9.3355, 9.3355, 10.7918, 10.7918, 0.0097, 9.3575,
])

sample_depths =  np.array([
5, 319, 5, 400, 5, 326, 5, 646, 5, 5
])


station_number = np.array([
'CTD 261 ',
'CTD 261 ',
'CTD 290 ',
'CTD 290 ',
'CTD 294 ',
'CTD 294 ',
'CTD 301 ',
'CTD 301 ',
'CTD 310 ',
'CTD 314 ',
])








proj_pp=ccrs.PlateCarree()
ax = plt.axes(projection=ccrs.Orthographic(5,68))
ax.coastlines()
ax.stock_img()


latsSurf = lats[np.where(sample_depths<10)]
lonsSurf = lons[np.where(sample_depths<10)]
stnoSurf = station_number[np.where(sample_depths<10)]

latsDeep = lats[np.where(sample_depths>10)]
lonsDeep = lons[np.where(sample_depths>10)]
stnoDeep = station_number[np.where(sample_depths>10)]


ax.plot(lonsSurf, latsSurf, marker='o', markersize=8, c='darkgreen', ls='none', transform=proj_pp, zorder=5, label='surface sample')
ax.plot(lonsDeep, latsDeep, marker='o', markersize=12, c='red', ls='none', transform=proj_pp, zorder=4, label='deep water sample')
for ii in range(len(lonsSurf)):
    ax.text(lonsSurf[ii], latsSurf[ii], stnoSurf[ii], horizontalalignment = 'right', transform=proj_pp, zorder=5)

ax.set_extent([-12.0, 20, 58, 73], proj_pp)
ax.legend()
if saveplots:
    plt.savefig(outfolder+'/locations.png',dpi=190, bbox_inches='tight')

if not saveplots:
    plt.show()
