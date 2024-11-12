
from starter2 import *
from scipy.ndimage import gaussian_filter
import pdb
import projector.proj as proj
import healpy as hp
import astropy
import dtools.davetools as dt
import shapely
from shapely.geometry import Polygon


def project(cube, xyz, dxyz, proj_center, proj_axis,bucket=None, molplot=False, moreplots=False, NSIDE = 4):

    
    verbose=True
    Nz = cube.size
    proj_center.shape=3,1
    xyz= xyz - proj_center

    #make the final map

    NPIX = hp.nside2npix(NSIDE)
    final_map = np.zeros(NPIX)
    

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
    if verbose: print('corners')
    corners = xyz+shifter*dxyz

    #we'll need this later.
    #Check for the zones whose projections will be squares.
    #Zones with only 4 corners in projection are those that contain the origin
    #along one coordiate. Thus, for zones whose unrotated coordinates have a mix of
    #positive and negative values will only have 4 corners as seen from the origin at 0.
    #four_corners = (np.abs(corners).sum(axis=2) - np.abs(corners.sum(axis=2)) >0).any(axis=0)
    four_corners = (np.abs(corners).sum(axis=2) - np.abs(corners.sum(axis=2)) >0).sum(axis=0) > 1

    #this can be streamlined, taken out of make_phi_theta
    if verbose: print('rotate')
    xyz_p = proj.rotate(xyz,proj_axis)

    #the thing to accumulate
    rrr2 =  (xyz_p**2).sum(axis=0).flatten()
    zone_volume = dxyz.prod(axis=0).flatten()
    zone_emission = cube/rrr2*zone_volume


    #the orthographic projection is used to determine the exterior corners.
    if verbose: print('work')
    #cor_p = proj.rotate(corners, proj_axis)
    #corners_oblique, phi_oblique, theta_oblique = proj.obliqueproj(xyz_p, cor_p)
    corners_persp,phi_persp,theta_persp=proj.make_phi_theta(corners, proj_axis)
    if 0:
        xyz_p.shape = 3,Nz
        phi_cen=phi_oblique.mean(axis=1)
        theta_cen=theta_oblique.mean(axis=1)
        phi_cen.shape=phi_cen.size,1
        theta_cen.shape=theta_cen.size,1

        #Decide if we're using perspective or orth projections
        if 0:
            theta_use=theta_oblique
            phi_use = phi_oblique
            distance_oblique = (phi_oblique-phi_cen)**2 + (theta_oblique-theta_cen)**2
            distance = distance_oblique
        else:
            theta_use=theta_persp
            phi_use = phi_persp
            distance_persp = (phi_persp-phi_cen)**2 + (theta_persp-theta_cen)**2
            distance = distance_persp
        #distance = np.maximum(distance_oblique, distance_persp)

        if verbose: print('more sort')
        asrt_distance = np.argsort(distance,axis=1)

        sorted_distance = np.take_along_axis(distance,asrt_distance,axis=1)
        sorted_theta= np.take_along_axis(theta_use,asrt_distance,axis=1)
        sorted_phi= np.take_along_axis(phi_use,asrt_distance,axis=1)
        ec_theta = sorted_theta[...,2:]
        ec_phi   = sorted_phi[...,2:]
        #if there's only 4 corners, don't keep two points
        keepers = np.ones_like(ec_phi,dtype='bool')
        keepers = keepers.T
        keepers[0,four_corners]=False
        keepers[1,four_corners]=False
        keepers = keepers.T

        #sort them clockwise.
        ec_theta_cen = ec_theta.mean(axis=1)
        ec_theta_cen.shape = Nz,1
        ec_phi_cen = ec_phi.mean(axis=1)
        ec_phi_cen.shape = Nz,1
        psi = np.arctan2(ec_theta-ec_theta_cen,ec_phi- ec_phi_cen)
        asrt_psi = np.argsort(psi,axis=1)
        ec_theta = np.take_along_axis(ec_theta, asrt_psi, axis=1)
        ec_phi =   np.take_along_axis(ec_phi, asrt_psi, axis=1)
        psi2 = np.take_along_axis(psi,asrt_psi,axis=1)
        keepers = np.take_along_axis(keepers,asrt_psi,axis=1)


    #From here, this needs to be in cython.
    import s2p_loop 
    return s2p_loop.zone_loop(Nz, zone_emission, theta_persp, phi_persp, NPIX,NSIDE, final_map)
