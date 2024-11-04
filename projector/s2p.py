
from starter2 import *
from scipy.ndimage import gaussian_filter
import pdb
import projector.proj as proj
import healpy as hp
import astropy


def project(cube, xyz, dxyz, proj_center, proj_axis):

    
    Nz = cube.size
    proj_center.shape=3,1
    xyz= xyz - proj_center
    if 1:
        #this shift puts them in an order that makes sense to draw lines
        shifter = np.array([[[-0.5, -0.5, +0.5, +0.5,  0.5, +0.5, -0.5, -0.5]],
                            [[-0.5, -0.5, -0.5, -0.5,  0.5, +0.5, +0.5, +0.5]],
                            [[-0.5,  0.5, +0.5, -0.5, -0.5, +0.5, +0.5, -0.5]]])

    dxyz.shape = 3,Nz,1
    xyz.shape = 3,Nz,1

    #get the corners of each zone.
    #Get projected theta and phi for each corner.
    #Get the 6 furthest from the center of projection
    #sort them in a clockwise fashion.
    corners = xyz+shifter*dxyz
    #pdb.set_trace()
    corners_new,phi_new,theta_new=proj.make_phi_theta(corners, proj_axis)
    phi_cen=phi_new.mean(axis=1)
    theta_cen=theta_new.mean(axis=1)
    phi_cen.shape=phi_cen.size,1
    theta_cen.shape=theta_cen.size,1
    phi = phi_new-phi_cen
    theta=theta_new-theta_cen
    distance = theta**2+phi**2
    asrt = np.argsort(distance,axis=1)
    sorted_theta= np.take_along_axis(theta_new,asrt,axis=1)
    sorted_phi= np.take_along_axis(phi_new,asrt,axis=1)
    ec_theta = sorted_theta[...,2:]
    ec_phi   = sorted_phi[...,2:]
    ec_theta_cen = ec_theta.mean(axis=1)
    ec_phi_cen = ec_phi.mean(axis=1)
    #ec_phi -= ec_phi_cen
    #ec_theta -= ec_theta_cen
    psi = np.arctan2(ec_theta-ec_theta_cen,ec_phi- ec_phi_cen)
    asrt = np.argsort(psi,axis=1)
    ec_theta = np.take_along_axis(ec_theta, asrt, axis=1)
    ec_phi =   np.take_along_axis(ec_phi, asrt, axis=1)
    psi2 = np.take_along_axis(psi,asrt,axis=1)


    NSIDE = 4
    for izone in range(ec_theta.shape[0]):
        #this works.
        edge_theta=ec_theta[izone]
        edge_phi  =ec_phi[izone]

        xyzpoly = astropy.coordinates.spherical_to_cartesian(1,np.pi/2-edge_theta,edge_phi)
        poly = np.array(xyzpoly).T

        my_pix = hp.query_polygon(NSIDE,poly, inclusive=True)
        #print(my_pix)

        if 1:
            #temp plot stuff
            plt.clf()
            NPIX = hp.nside2npix(NSIDE)
            m = np.arange(NPIX)
            m[my_pix]=m.max()
            hp.mollview(m, title="Mollview image RING")
            hp.projscatter(edge_theta,edge_phi, c='r')

        #code that works to find boundaries.
        for ipix in my_pix:
            ######yes this works
            xyz = hp.boundaries(NSIDE, ipix, step=1)
            theta, phi = hp.vec2ang(xyz.T)
            ######
            hp.projscatter(theta,phi)

        prefix='%s/moltest'%plot_dir
        nplots = len(glob.glob(prefix+"*"))
        plt.savefig(prefix+"%03d"%nplots)
        
        if 1:
            fig,axes=plt.subplots(1,2,figsize=(8,4))
            ax0=axes[0];ax1=axes[1]
            ax0.set_aspect('equal')
            ax1.set_aspect('equal')
            ax0.plot(corners[0][izone], corners[1][izone],c='r')
            ax0.plot(corners_new[0][izone], corners_new[1][izone],marker='*')
            ax1.plot(phi_new[izone],theta_new[izone])
            ax1.set(ylim=[np.pi,0], xlim=[np.pi, -np.pi])
            ax1.plot(edge_phi, edge_theta,c='orange')
            prefix = '%s/corners'%plot_dir
            nplot = len(glob.glob(prefix+"*"))
            fig.savefig(prefix+"%03d"%nplot)
            #pdb.set_trace()
