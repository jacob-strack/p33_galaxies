#A stride-one perspective projector
#that can deal with heavily oversampled rays.
#Somewhat approximate at the moment. 
#Not a general geometry.

#1.) Find the (phi, theta) bins for every corner of every zone.
#2.) Quantize the phi and theta bin arrays
#3.) Each zone is covered by Ntheta * Nphi bins.  Loop over Ntheta and Nphi, fill the hist.


from starter2 import *
from scipy.ndimage import gaussian_filter

def make_xyz(cube):
    """presently hard coded to project along Y"""
    N = cube.shape[0]
    dx_v = 1./N
    center = nar([0.5,-0.2,0.5])
    xyz = np.stack(np.mgrid[0:1:dx_v,0:1:dx_v,0:1:dx_v])
    return xyz, dx_v

def make_phi_theta(xyz,center,projax):
    """ presently not very general"""
    center.shape=(3,1,1,1)
    xyz_p = xyz-center
    r = np.sqrt((xyz_p**2).sum(axis=0))
    theta = np.arccos( xyz_p[2]/r )
    phi   = np.arctan2( xyz_p[1], xyz_p[0])
    return r, phi, theta

class s1p():
    def __init__(self,cube, center=nar([0.5,-0.2,0.5]), projax=[0,1,0], verbose=False, Nbins=None):
        self.cube=cube
        self.proj(center, projax, verbose, Nbins)
    def proj(self, center=nar([0.5,-0.2,0.5]),projax=[0,1,0], verbose=False, Nbins=None):
        cube=self.cube
        if Nbins is None:
            Nbins = cube.shape[0]

        xyz, dx_v = make_xyz(cube)
        #center works, projax is a dummy variable.  Presently hard coded along y.
        r,phi,theta = make_phi_theta(xyz,center,projax)
        cube=self.cube

        #1.) Find the (phi, theta) bins for every corner of every zone.

        got=False
        for sX in [-0.5,0.5]:
            for sY in [-0.5,0.5]:
                for sZ in [-0.5,0.5]:
                    shift = np.stack([sX,sY,sZ])*dx_v
                    shift.shape=(3,1,1,1)
                    xyz_shift = xyz+shift
                    r_shift, this_phi, this_theta = make_phi_theta(xyz_shift, center, projax)
                    if not got:
                        max_theta = this_theta
                        min_theta = this_theta
                        max_phi = this_phi
                        min_phi = this_phi
                        got=True
                    else:
                        max_theta = np.maximum(max_theta, this_theta)
                        min_theta = np.minimum(min_theta, this_theta)
                        max_phi   = np.maximum(max_phi, this_phi)
                        min_phi   = np.minimum(min_phi, this_phi)


        #2.) Quantize the phi and theta bin arrays
        #the eps is to make sure the last point isn't in the Nbins+1 bin.

        eps = 1e-15
        minmin_theta = min_theta.min()-eps
        minmin_phi = min_phi.min() - eps
        maxmax_theta = max_theta.max() +eps
        maxmax_phi = max_phi.max() + eps
        dtheta = (maxmax_theta-minmin_theta)/(Nbins)
        dphi   = (maxmax_phi - minmin_phi)/(Nbins)
        min_theta_bin =((min_theta - minmin_theta)//dtheta).astype('int')
        min_phi_bin =  ((min_phi - minmin_phi)//dphi).astype('int')
        max_theta_bin =((max_theta - minmin_theta)//dtheta).astype('int')
        max_phi_bin =  ((max_phi - minmin_phi)//dphi).astype('int')
        Ntheta = max_theta_bin-min_theta_bin
        Nphi   = max_phi_bin - min_phi_bin
        Ntheta_max = Ntheta.max()
        Nphi_max = Nphi.max()
        if verbose:
            print("Ntheta across, Npphi",Ntheta_max, Nphi_max)
            print("max theta bin",max_theta_bin.max())
        Ntheta_bins=Nbins
        Nphi_bins = Nbins


        #set up target array H and mask array F
        #and theta and phi bins
        H = np.zeros(Nphi_bins*Ntheta_bins)
        F = np.zeros(Nphi_bins*Ntheta_bins)
        F[:]=np.nan
        theta_bins = np.linspace( minmin_theta, maxmax_theta, Ntheta_bins)
        phi_bins = np.linspace( minmin_phi, maxmax_phi, Nphi_bins)
        theta_cen = 0.5*(theta_bins)
        phi_cen = 0.5*(phi_bins)
        coordPhi, coordTheta = np.meshgrid(phi_cen,theta_cen, indexing='ij')

        r1 = range(Nphi_max)
        r2 = range(Ntheta_max)
        
        #3.) Each zone is covered by Ntheta * Nphi bins.  Loop over Ntheta and Nphi, fill the hist.
        #the mask is to not overstep a zone's influence.
        cube_ind = np.arange(cube.size)
        cube_flat = cube.flatten()
        for i in r1:
            for j in r2:
                if verbose:
                    print('loop',i,j,Ntheta_max)
                #which zones to take
                mask = ((Ntheta >= i)*(Nphi >= j)).flatten()
                #the position in the 2d histogram H of the 3d zone
                index  = (min_theta_bin+i + Ntheta_bins*(min_phi_bin +j)).flatten()
                #pre compute some things
                cube_ind_mask = cube_ind[mask]
                mask_ind=index[mask]
                #Do the actual filling of the histogram.
                #Loops are to be avoided, this should be at least re done in cython.
                #I tried to do this with a mask, but it failed.  
                #Here is a place for speed improvements.
                for Hi, Ii in enumerate(cube_ind[mask]):
                    H[mask_ind[Hi]] += cube_flat[Ii]
                #The mask of found pixels.
                F[index[mask]]=0
        not_ok = np.isnan(F)
        #H[not_ok]=F[not_ok]
        H.shape = Nphi_bins,Ntheta_bins
        F.shape = Nphi_bins,Ntheta_bins

        self.coordPhi=coordPhi
        self.coordTheta=coordTheta
        self.H = H
        self.mask = F

def plot_image(coordPhi, coordTheta, Hin, fname, mask=None):
    H = Hin + 0 #make a copy
    #nrm = mpl.colors.Normalize(vmin=den[ok][den[ok]>0].min(),vmax=den[ok].max())
    ok = ~np.isnan(H)
    nrm = mpl.colors.Normalize( vmin=H[ok].min(), vmax=H[ok].max())
    #nrm = mpl.colors.Normalize(vmin=0,vmax=4)
    cmap = copy.copy(mpl.cm.get_cmap("Reds"))
    cmap.set_under('g')
    cmap.set_bad('k')
    
    if mask is not None:
        not_ok = np.isnan(mask)
        H[not_ok] = np.nan

    fig,axes=plt.subplots(1,1)
    ax0=axes
    p=ax0.pcolormesh(coordPhi, coordTheta, H, cmap=cmap,norm=nrm )
    fig.colorbar(p,ax=ax0)
    
    fig.savefig('%s/test2'%plot_dir)
