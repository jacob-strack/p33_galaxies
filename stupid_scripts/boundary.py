from astropy_healpix import HEALPix

hp = HEALPix(nside=32, order='nested')

print(hp.boundaries_lonlat([120], step=1) )
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from astropy_healpix.core import boundaries_lonlat

ax = plt.subplot(1, 1, 1)
plt.clf()
for step, color in [(1,'red')]:#[(1, 'red'), (100, 'black')]:
    lon, lat = boundaries_lonlat([5952], nside=1, step=step)
    lon = lon.to(u.deg).value
    lat = lat.to(u.deg).value
    vertices = np.vstack([lon.ravel(), lat.ravel()]).transpose()
    p = Polygon(vertices, closed=True, edgecolor=color, facecolor='none')
    ax.add_patch(p)

plt.xlim(210, 330)
plt.ylim(-50, 50)
plt.savefig('%s/test'%plot_dir)
