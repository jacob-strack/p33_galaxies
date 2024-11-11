
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
    cube, xyz, dxyz = proj.make_cube_full(32, stick_or_sphere=0)
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

            image=s2p.project(cube,xyz,dxyz,proj_center,proj_axis, bucket=bucket,molplot=molplot, NSIDE=128)

            plt.clf()
            hp.mollview(image, title="Two Sticks")
            prefix='%s/proj_sticks'%plot_dir
            nplots = len(glob.glob(prefix+"*"))
            plt.savefig(prefix+"%03d"%nplots)


if 0:
    ds = yt.load('datasets/IsolatedGalaxy/galaxy0030/galaxy0030')
    #ad = ds.all_data()
    c = ds.arr([0.5]*3, 'code_length')
    r = ds.quan(1./16, 'code_length')
    L = ds.arr([0.5-1./32]*3,'code_length')
    R = ds.arr([0.5+1./32]*3,'code_length')
    dx = ds.index.get_smallest_dx()*16
    Nz = (R-L)/dx
    cg = ds.covering_grid(4,L,Nz)
    ad=cg
    #sphere = ds.sphere(c, r)
    cube = ad['density'].in_units('code_density').flatten().v
    space_units='code_length'
    xyz = np.stack([ad['x'].in_units(space_units).v.flatten(),
                    ad['y'].in_units(space_units).v.flatten(),
                    ad['z'].in_units(space_units).v.flatten()])
    dxyz = np.stack([ad['dx'].in_units(space_units).flatten().v,
                     ad['dy'].in_units(space_units).flatten().v,
                     ad['dz'].in_units(space_units).flatten().v])
    if 0:
        density=cg['density']
        print("Nz %0.2e"%density.size)
        fig,ax=plt.subplots(1,1)
        ax.imshow( np.log(density.sum(axis=2)))
        fig.savefig("%s/galaxy_cg"%plot_dir)

    if 0:
        proj=ds.proj('density',2, data_source=sphere)
        pw=proj.to_pw()
        #pw.zoom(16)
        pw.save('%s/galaxy'%plot_dir)

    dx_min = dxyz.min()

    proj_center = np.array([0.5]*3)+np.array([0,0,-1./32]) + dx_min/128
    proj_axis = np.array([0,0,1.])
    image=s2p.project(cube,xyz,dxyz,proj_center,proj_axis, bucket=bucket,molplot=molplot, NSIDE=128)

    plt.clf()
    hp.mollview(image, title="Isolated Galaxy")
    prefix='%s/galaxy'%plot_dir
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

