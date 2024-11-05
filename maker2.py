
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
    dxyz/=8
    old_center = nar([0.5]*3)
    old_center.shape = old_center.size,1
    new_center = nar([0.0,0.0,0.0])
    new_center.shape=new_center.size,1
    xyz += new_center-old_center
    proj_center = nar([0.5,0.7,0.6])
    if 0:
        #proj_center = nar([0.5]*3)
        #proj_center = nar([0,0,0])
        #proj_axis   = nar([0,1,0],dtype='float')
        dcenter = proj_center-new_center.flatten()
        proj_axis_tmp= -dcenter/(dcenter**2).sum()
        proj_axis = np.cross(proj_axis_tmp,[0,0,1])
        proj_axis=proj_axis_tmp

    #proj_axis=nar([0,1,1],dtype='float')
    #proj_axis/=(proj_axis**2).sum()

    if 1:
        #slew the camera around the zone
        theta = np.linspace(np.pi/7,np.pi*5/6,10)-np.pi/2
        phi   = np.linspace(-np.pi,np.pi,10)
        theta = theta[5:6]
        phi = phi[1:2]
        r=2
        for nph,ph in enumerate(phi):
            for nth,th in enumerate(theta):
                print(nph,nth)
                x=r*np.sin(th)*np.cos(ph)
                y=r*np.sin(th)*np.sin(ph)
                z=r*np.cos(th)
                proj_center=nar([x,y,z])
                dcenter = proj_center-new_center.flatten()
                proj_axis = -dcenter/(dcenter**2).sum()
                bucket={'theta':th,'phi':ph}
                corners=s2p.project(cube,xyz,dxyz,proj_center,proj_axis, bucket=bucket)



    if 0:
        #rotate the camera, fixed in space
        theta = np.linspace(np.pi/7,np.pi*5/6,10)-np.pi/2
        phi   = np.linspace(-np.pi,np.pi,10)
        for ph in phi:
            for th in theta:
                proj_axis=nar(astropy.coordinates.spherical_to_cartesian(1,th,ph))
                print(proj_axis)
                bucket={'theta':th,'phi':ph}
                corners=s2p.project(cube,xyz,dxyz,proj_center,proj_axis, bucket=bucket)

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



