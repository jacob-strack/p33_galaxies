#Loop as is in projector code
import numpy as np

def loop_H(H, mask_ind,mask,cube_ind,cube_flat): 
    for Hi, Ii in enumerate(cube_ind[mask]):
        H[mask_ind[Hi]] += cube_flat[Ii]
    return H

def big_loop(H,F,cube,Nphi_max,Ntheta_max,min_phi_bin,min_theta_bin,Ntheta_bins,Nphi,Ntheta,verbose=False):
    r1 = range(Nphi_max)
    r2 = range(Ntheta_max)
    cube_ind = np.arange(cube.size)
    cube_flat = cube.flatten()
    for i in r1:
        for j in r2:
            if verbose:
                print('loop',i,j,Nphi_max,Ntheta_max)
            #which zones to take
            mask = ((Ntheta > j)*(Nphi > i)).flatten() #swapped i and j took out equals 
            #the position in the 2d histogram H of the 3d zone
            index  = ((min_theta_bin + j) + (Ntheta_bins*(min_phi_bin + i))).flatten() #swapped i and j 
            #pre compute some things
            cube_ind_mask = cube_ind[mask]
            mask_ind=index[mask]
            #print(mask_ind.max(), np.shape(H))
            #Do the actual filling of the histogram.
            #Loops are to be avoided, this should be at least re done in cython.
            #I tried to do this with a mask, but it failed.  
            #Here is a place for speed improvements.
            import loop
            H = loop_H(H, mask_ind, mask, cube_ind, cube_flat)
            #for Hi, Ii in enumerate(cube_ind[mask]):
                #H[mask_ind[Hi]] += cube_flat[Ii]
            #The mask of found pixels.
            F[index[mask]]=0
    return H,F
