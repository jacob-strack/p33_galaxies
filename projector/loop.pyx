#Loop as is in projector code
import numpy as np
import healpy as hp
import yt
import matplotlib.pyplot as plt
def field_to_healpix(H, field_ind_mask,mask,healpix_ind,field):  
    print("H loop shape", np.shape(H))
    for i in field_ind_mask:
        H[healpix_ind[i]] += field[i] #fill healpix map with field values
    return H

def loop_H(H,mask_ind, mask, cube_ind, field): 
    for Hi, Ii in enumerate(cube_ind[mask]):
        H[mask_ind[Hi]] += field[Ii]
    return H 

def enzo_to_healpix(arr, theta, phi, nside):
    arr = arr.flatten()
    theta = theta.flatten()
    phi = phi.flatten()
    hp_inds = hp.ang2pix(nside, theta, phi)
    m = yt.YTArray(np.zeros(hp.nside2npix(nside)), arr.units)
    for i in range(len(arr)):
        if m[hp_inds[i]] == 0:
            m[hp_inds[i]] = arr[i]

        else: 
            m[hp_inds[i]] = np.mean(arr[i], m[hp_inds[i]])
    return m

def big_loop(H,F,field,Nphi_max,Ntheta_max,min_phi_bin,min_theta_bin,Ntheta_bins,Nphi,Ntheta,dtheta, dphi,max_theta_bin,max_phi_bin,min_theta,max_theta,min_phi,max_phi,verbose=False,method="projector"):
    r1 = range(Nphi_max)
    r2 = range(Ntheta_max)
    field_ind = np.arange(len(field))
    print(Ntheta_bins)
    m = yt.YTArray(np.zeros(hp.nside2npix(15)),field.units) #TODO: pass in nside as a parameter to big_loop
    for i in r1:
        for j in r2:
            if verbose:
                print('loop',i,j,Nphi_max,Ntheta_max)
            #which zones to take
            mask = ((Ntheta > j)*(Nphi > i)).flatten() #swapped i and j took out equals 
            #the position in the 2d histogram H of the 3d zone
            index  = ((min_theta_bin + j) + (Ntheta_bins*(min_phi_bin + i))).flatten() #swapped i and j 
            #pre compute some things
            #healpix index array goes here
            #theta and phi here use the already binned quantities, maybe loses information
            #theta = (min_theta_bin + max_theta_bin)*dtheta/2
            #phi = (min_phi_bin + max_phi_bin)*dphi/2
            #theta and phi using angles from zones (i think this is the right thing?)
            phi = (min_phi + max_phi)/2 
            theta = (min_theta + max_theta)/2
            print("min/max phi loop", np.min(phi), np.max(phi))
            print("min/max theta loop", np.min(theta), np.max(theta))
            healpix_ind = hp.ang2pix(15, theta.value, phi.value)
            if 1:
                import matplotlib.pyplot as plt
                plt.clf()
                plt.hist(healpix_ind, bins = 100)
                plt.savefig("healpixinds.png")
                plt.clf()
                plt.hist(theta.value, bins = 100)
                plt.savefig("thetas.png")
                plt.clf()
                plt.hist(phi.value, bins = 100)
                plt.savefig("phis.png")
                plt.clf()
            field_ind_mask = field_ind[mask]
            mask_ind=index[mask]
            #print(mask_ind.max(), np.shape(H))
            #Do the actual filling of the histogram.
            #Loops are to be avoided, this should be at least re done in cython.
            #I tried to do this with a mask, but it failed.  
            #Here is a place for speed improvements.
            import loop
            if method=="healpix":
                m = field_to_healpix(m, field_ind_mask, mask, healpix_ind, field) 
            if method=="projector": 
                H = loop_H(H, mask_ind, mask, field_ind, field) 
            #for Hi, Ii in enumerate(cube_ind[mask]):
                #H[mask_ind[Hi]] += cube_flat[Ii]
            #The mask of found pixels.
            F[index[mask]]=0

    print("map max", np.max(m))
    import matplotlib.pyplot as plt
    plt.clf()
    plt.hist(healpix_ind,bins = 100)
    plt.savefig("healpixinds.png")
    plt.clf()
    return H,F,m
