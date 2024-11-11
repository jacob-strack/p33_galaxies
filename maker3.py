
import yt
from starter2 import *
import pdb
import projector.s2p as s2p
import projector.proj as proj
import healpy as hp
reload(s2p)
reload(proj)

if 1:
    #strafe the camera around the zone
    cube, xyz, dxyz = proj.make_cube_full(32, stick_or_sphere=1)
    molplot=False
    theta = np.linspace(np.pi/7,np.pi*5/6,10)-np.pi/2
    phi   = np.linspace(-np.pi,np.pi,10)
    #theta = theta[9:10] #broken
    theta = theta[0:1]
    phi = phi[0:1]
    cube_center = xyz.mean(axis=1)

    r=2
    for nph,ph in enumerate(phi):
        for nth,th in enumerate(theta):
            print("Nth, Nph",nph,nth)
            x=r*np.sin(th)*np.cos(ph)
            y=r*np.sin(th)*np.sin(ph)
            z=r*np.cos(th)
            proj_center=nar([x,y,z])
            dcenter = proj_center-cube_center.flatten()
            proj_axis = -dcenter/(dcenter**2).sum()
            bucket={'theta':th,'phi':ph}

            image=s2p.project(cube,xyz,dxyz,proj_center,proj_axis, bucket=bucket,molplot=molplot, NSIDE=64)

            plt.clf()
            hp.mollview(image, title="Two Sticks")
            prefix='%s/proj_sticks'%plot_dir
            nplots = len(glob.glob(prefix+"*"))
            plt.savefig(prefix+"%03d"%nplots)



if 0:
    cube, xyz, dxyz = proj.make_cube_full(2)
    #dxyz/=8
    old_center = nar([0.5]*3)
    old_center.shape = old_center.size,1
    new_center = nar([0.0,0.0,0.0])
    new_center.shape=new_center.size,1
    #xyz += new_center-old_center
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

