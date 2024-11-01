
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
    corners = xyz+shifter*dxyz
    corners_new,phi_new,theta_new=proj.make_phi_theta(corners, proj_axis)
    phi_cen=phi_new.mean(axis=1)
    theta_cen=theta_new.mean(axis=1)
    phi_cen.shape=phi_cen.size,1
    theta_cen.shape=theta_cen.size,1
    phi = phi_new-phi_cen
    theta=theta_new-theta_cen
    distance = theta**2+phi**2
    asrt = np.argsort(distance,axis=1)
    ec_theta = np.take_along_axis(theta_new,asrt,axis=1)[...,2:]
    ec_phi   = np.take_along_axis(phi_new,asrt,axis=1)[...,2:]
    ec_theta_cen = ec_theta.mean(axis=1)
    ec_phi_cen = ec_phi.mean(axis=1)
    ec_phi -= ec_phi_cen
    ec_theta -= ec_theta_cen
    psi = np.arctan2(ec_theta,ec_phi)
    asrt = np.argsort(psi,axis=1)
    ec_theta = np.take_along_axis(ec_theta, asrt, axis=1)
    ec_phi =   np.take_along_axis(ec_phi, asrt, axis=1)
    psi2 = np.take_along_axis(psi,asrt,axis=1)

    #KLUDGE here we break with the truly stride one operations.

        #plt.clf()
        #plt.scatter(phi_new,theta_new)
        #plt.scatter(ec_phi, ec_theta)
        #plt.plot(phi_new[0,:], theta_new[0,:])
        #plt.savefig('%s/theta_phi'%plot_dir)
        #plt.clf()
        #plt.scatter(corners_new[0].flatten(),corners_new[1].flatten(),c='r')
        #plt.plot(corners_new[0,0,:], corners_new[1,0,:])
        #plt.savefig('%s/corners'%plot_dir)
        #plt.clf()

    NSIDE = 32
    #https://github.com/healpy/healpy/issues/393
    #
    for izone in range(ec_theta.shape[0]):
        #this works.
        edge_theta=ec_theta[izone]
        edge_phi  =ec_phi[izone]
        xyzpoly = astropy.coordinates.spherical_to_cartesian(1,edge_theta,edge_phi)
        poly = np.array(xyzpoly).T

        ipix = hp.query_polygon(NSIDE,poly, inclusive=True)

        from astropy_healpix import HEALPix
        hp2 = HEALPix(nside=NSIDE, order='ring')
        n = 0
        for pix in ipix[n:n+1]:
            #THIS DOES NOT WORK.
            #edges = hp.boundaries(nside=NSIDE,pix=pix,step=1)
            ra, dec = hp2.boundaries_lonlat([pix],step=1)
            #print(edges.shape)
            #angles = astropy.coordinates.cartesian_to_spherical(*edges)
            #print(angles)
            #ra,dec=hp.pixelfunc.vec2ang(edges.T, lonlat=True)
            #rad,ra,dec=angles

            #print('edges',edges)
        #pdb.set_trace()
        #This works except for ra,dec.
        plt.clf()
        plt.plot(edge_phi,edge_theta,c='g')
        for i in range(len(edge_phi)):
            plt.text(edge_phi[i], edge_theta[i], "%f"%psi2[izone][i])
        plt.savefig('%s/ec_sort'%plot_dir)
        NPIX = hp.nside2npix(NSIDE)
        m = np.arange(NPIX)
        m[ipix]=m.max()
        hp.mollview(m, title="Mollview image RING")
        plt.plot(edge_theta,edge_phi,c='r')
        print('theta',edge_theta)
        print('phi  ',edge_phi)
        #plt.plot(angles[1],angles[2],c='k')
        plt.plot(ra,dec,c='orange',linewidth=11)
        print('dec ',dec)
        print('ra  ',ra)
        plt.savefig('%s/moltest'%plot_dir)
        #pdb.set_trace()


