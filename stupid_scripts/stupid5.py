

NSIDE = 4
NPIX = hp.nside2npix(NSIDE)
m = np.arange(NPIX)

import dtools.davetools as dt
lon_ext = dt.extents()
lat_ext = dt.extents()
for ipix in m:
    xyz = hp.boundaries(NSIDE, ipix, step=1)
    theta, phi = hp.vec2ang(xyz.T)
    lon_ext(theta)
    lat_ext(phi)
print(lon_ext)
print(lat_ext)


