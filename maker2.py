
import yt
from starter2 import *
import pdb
import projector.s2p as s2p
import projector.proj as proj
#import healpy as hp
reload(s2p)
reload(proj)

if 1:
    cube, xyz, dxyz = proj.make_cube_full(10)
    dxyz/=2
    proj_center = nar([0.5,0.5,0.5])
#proj_center = nar([0,0,0])
    proj_axis   = nar([0,1,0],dtype='float')
    corners=s2p.project(cube,xyz,dxyz,proj_center,proj_axis)

if 0:
    wmap_map_I = hp.read_map("wmap_band_iqumap_r9_7yr_W_v4.fits")
    hp.mollview(
            wmap_map_I,
            coord=["G", "E"],
            title="Histogram equalized Ecliptic",
            unit="mK",
            norm="hist",
            min=-1,
            max=1,

    )
    plt.savefig('derp')



