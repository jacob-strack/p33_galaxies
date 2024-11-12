from starter2 import *
from scipy.ndimage import gaussian_filter
import pdb
import projector.proj as proj
import healpy as hp
import astropy
import dtools.davetools as dt
import shapely
from shapely.geometry import Polygon
def zone_loop(Nz, zone_emission, theta_persp, phi_persp, NPIX, NSIDE,final_map):
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


    print("final_map", final_map)
    return final_map
