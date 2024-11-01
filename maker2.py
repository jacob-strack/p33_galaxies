
import yt
from starter2 import *
import pdb
import projector.s2p as s2p
import projector.proj as proj
import healpy as hp
reload(s2p)
reload(proj)

if 1:
    cube, xyz, dxyz = proj.make_cube_full(1)
    dxyz/=128
    new_center = nar([0.15,0.0,-0.3])
    new_center.shape=new_center.size,1
    xyz += new_center
    proj_center = nar([.6, -1.5, 1.5])
#proj_center = nar([0,0,0])
    #proj_axis   = nar([0,1,0],dtype='float')
    dcenter = nar([0.5]*3)-proj_center
    proj_axis = -dcenter/(dcenter**2).sum()
    #proj_axis=nar([0,1,1],dtype='float')
    #proj_axis/=(proj_axis**2).sum()

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



