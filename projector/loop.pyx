#Loop as is in projector code
import numpy as np
import healpy as hp
import yt
def loop_H(H, field_ind_mask,mask,healpix_ind,field):  
  #  for Hi, Ii in enumerate(cube_ind[mask]):
  #      print(mask_ind[Hi], Ii)
  #      H[mask_ind[Hi]] += field[Ii]
    print("H loop shape", np.shape(H))
    for Hi,Fi in enumerate(field_ind_mask):
        H[healpix_ind[Hi]] += field[Fi] #fill healpix map with field values 
    return H

def big_loop(H,F,field,Nphi_max,Ntheta_max,min_phi_bin,min_theta_bin,Ntheta_bins,Nphi,Ntheta,theta_cen, phi_cen,verbose=False):
    r1 = range(Nphi_max)
    r2 = range(Ntheta_max)
    field_ind = np.arange(len(field))
    print(Ntheta_bins)
    m = yt.YTArray(np.zeros(hp.nside2npix(70)),field.units) #TODO: pass in nside as a parameter to big_loop
    healpix_ind = hp.ang2pix(70, theta_cen.value, phi_cen.value) #array which contains healpix map indices in input field ordering
    for i in r1:
        for j in r2:
            if verbose:
                print('loop',i,j,Nphi_max,Ntheta_max)
            #which zones to take
            mask = ((Ntheta > j)*(Nphi > i)).flatten() #swapped i and j took out equals 
            #the position in the 2d histogram H of the 3d zone
            index  = ((min_theta_bin + j) + (Ntheta_bins*(min_phi_bin + i))).flatten() #swapped i and j 
            #pre compute some things
            field_ind_mask = field_ind[mask]
            mask_ind=index[mask]
            hp_ind_mask = healpix_ind[mask]
            print("unique healpix inds", np.unique(hp_ind_mask))
            #print(mask_ind.max(), np.shape(H))
            #Do the actual filling of the histogram.
            #Loops are to be avoided, this should be at least re done in cython.
            #I tried to do this with a mask, but it failed.  
            #Here is a place for speed improvements.
            import loop
            print("H before loop shape", np.shape(H))
            m = loop_H(m, field_ind_mask, mask, healpix_ind, field) 
            #for Hi, Ii in enumerate(cube_ind[mask]):
                #H[mask_ind[Hi]] += cube_flat[Ii]
            #The mask of found pixels.
            F[index[mask]]=0
    return H,F,m
