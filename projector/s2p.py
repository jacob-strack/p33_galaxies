
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
    for izone in range(Nz):
        print("Izone/Nzone %d/%d = %0.2f"%(izone,Nz, izone/Nz) )
        quantity = zone_emission[izone] #q*V/r^2
        all_theta = theta_persp[izone]
        all_phi   = phi_persp[izone]

        all_corners = np.stack([all_theta,all_phi])
        corner_poly = Polygon(all_corners.T)
        hull = shapely.convex_hull(corner_poly)
        edge_theta, edge_phi = hull.exterior.xy
        #the polygon comes out closed.
        edge_theta = np.array(edge_theta)[:-1]
        edge_phi = np.array(edge_phi)[:-1]
        #pdb.set_trace()
        if 0:
            edge_theta=ec_theta[izone]
            edge_phi  =ec_phi[izone]
            #sometimes 6, sometimes 4 corners.
            #This is a problem for truly stride-one vectorization.
            edge_theta = edge_theta[keepers[izone]]
            edge_phi = edge_phi[keepers[izone]]

        #healpy can't deal if there are colinear points.
        #Compute the area for every set of three points.  
        #If there are colinear points, reject the middle one.
        if 1:
            for nroll in range(len(edge_theta)):
                it=0
                area = edge_theta[it  ]*(edge_phi[it+1]-edge_phi[it+2])+\
                       edge_theta[it+1]*(edge_phi[it+2]-edge_phi[it+0])+\
                       edge_theta[it+2]*(edge_phi[it+0]-edge_phi[it+1])
                if np.abs(area) < 1e-6:
                    edge_theta = np.delete(edge_theta,1)
                    edge_phi = np.delete(edge_phi,1)
                else:
                    edge_theta = np.roll(edge_theta,1)
                    edge_phi   = np.roll(edge_phi,1)

            #check for convexity
            #https://stackoverflow.com/questions/471962/how-do-i-efficiently-determine-if-a-polygon-is-convex-non-convex-or-complex
            zxprod = np.zeros(edge_phi.size)
            x=edge_theta
            y=edge_phi
            #if zprod changes sign, its not convex
            for k in np.arange(zxprod.size)-2:
                dx1 = x[k+1]-x[k]
                dy1 = y[k+1]-y[k]
                dx2 = x[k+2]-x[k+1]
                dy2 = y[k+2]-y[k+1]
                zxprod[k] = dx1*dy2 - dy1*dx2
            if np.abs(zxprod).sum() - np.abs(zxprod.sum()) > 0:
                print("Convex Structure, FIX ME VERY BAD")
                continue

        if 1:
            #compute the polygon in 3d
            xyzpoly = astropy.coordinates.spherical_to_cartesian(1,np.pi/2-edge_theta,edge_phi)
            xyzpoly = np.array(xyzpoly)
            poly = np.array(xyzpoly).T

        #check for degenerate corners.  
        #for some reason I'm still getting a degenerate corner when I don't think I should.
        #https://healpix.sourceforge.io/html/Healpix_cxx/healpix__base_8cc_source.html line 1000
        if 1:
            degenerate=False
            for i in np.arange(poly.shape[1])-2:
                normal = np.cross(poly[i], poly[i+1])
                hnd = np.dot(normal, poly[i+2])
                if np.abs(hnd) < 1e-10:
                    print("Degenerate Corner FIX ME",i)
                    degenerate=True
            if degenerate:
                continue

        #the all important pixel query

        try:
            my_pix = hp.query_polygon(NSIDE,poly, inclusive=True)
        except:
            print("MISSED A BAD CORNER")
            continue

        #we'll also need this.  Can be streamlined.
        zone_poly = np.stack([edge_theta,edge_phi])
        q = Polygon(zone_poly.T)
        zone_area = q.area
        zone_column_density = quantity/zone_area

        #plot the points we'll use.
        if moreplots:
            fig,ax1=plt.subplots(1,1,figsize=(8,8))
            ax1.plot(phi_use[izone],theta_use[izone])
            ax1.plot(edge_phi, edge_theta,c='r')
            #ax1.set(ylim=[np.pi,0], xlim=[np.pi, -np.pi])
            #ax1.plot(edge_phi, edge_theta,c='orange')
            ax1.scatter(edge_phi, edge_theta,c='orange')
            ax1.set(title="th %0.2f ph %0.2f"%(bucket['theta'],bucket['phi']))
            ax1.set(xlabel='phi',ylabel='theta')

            ax1.scatter(phi_cen[izone],theta_cen[izone],c='r')
            d = distance[izone]
            for nd, dd in enumerate(d):
                ax1.text(phi_use[izone][nd], theta_use[izone][nd],"%0.3e"%dd)

            prefix = '%s/cubeproj'%plot_dir
            nplot = len(glob.glob(prefix+"*"))
            fig.savefig(prefix+"%03d"%nplot)
            plt.close(fig)
            #pdb.set_trace()
        
            #plot the polygon
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            ax.plot3D(*poly.T)
            ax.scatter3D(*poly.T)
            for p in poly:
                b = np.stack([[0,0,0],p])
                ax.plot3D(*b.T)
            for psi2 in np.linspace(0,90,5):
                ax.view_init(elev=27,azim=psi2)
                prefix = '%s/3d'%plot_dir
                nplot = len(glob.glob(prefix+"*"))
                fig.savefig(prefix+"%03d"%nplot)
            plt.close(fig)


        if molplot:
            #temp plot stuff
            plt.clf()
            m = np.arange(NPIX)
            m[my_pix]=m.max()
            hp.mollview(m, title="Mollview image RING")
            hp.projscatter(edge_theta,edge_phi, c='r')

        #code that works to find boundaries.
        for ipix in my_pix:
            ######yes this works
            xyz = hp.boundaries(NSIDE, ipix, step=1)
            theta, phi = hp.vec2ang(xyz.T)

            ray_poly = np.stack([theta,phi])
            p = Polygon(ray_poly.T)
            intersection = p.intersection(q)
            area = intersection.area

            #the final quantity = q*V/A*intersection/r^2
            net_light=area*zone_column_density
            final_map[ipix] += net_light
            #final_map[ipix] = 1

            #
            # plotting 
            #
            if molplot:
                hp.projscatter(theta,phi,c='k')

            if moreplots:
                fig,ax=plt.subplots(1,1)
                ax.plot(*q.exterior.xy)
                ax.plot(*p.exterior.xy)
                ax.plot(*intersection.exterior.xy)
                #ax.plot(theta,phi)
                ##ax.plot(edge_theta,edge_phi)
                #pdb.set_trace()
                prefix = '%s/intersector_'%plot_dir
                nplot = len(glob.glob(prefix+"*"))
                fig.savefig(prefix+"%03d"%nplot)

            #if I decide to roll my own:
            #2.) InZone
            #3.) CrossInPixel
            #4.) CrossInZone
            #5.) Sort both
            #6.) Intersection 1
            #7.) Intersection 2
            #8.) Collect interior points
            #9.) sort interior points clockwise
            #10.) area of interior points
            #11.) Add the right thing to the destination plot.
            ######


        if molplot:
            #finish plotting
            prefix='%s/moltest'%plot_dir
            nplots = len(glob.glob(prefix+"*"))
            plt.savefig(prefix+"%03d"%nplots)
    return final_map
        
