#A stride-one perspective projector
#that can deal with heavily oversampled rays.
#Somewhat approximate at the moment. 
#Not a general geometry.

#1.) Find the (phi, theta) bins for every corner of every zone.
#2.) Quantize the phi and theta bin arrays
#3.) Each zone is covered by Ntheta * Nphi bins.  Loop over Ntheta and Nphi, fill the hist.

from starter2 import *
from scipy.ndimage import gaussian_filter
import pdb

def make_xyz(cube):
    N = cube.shape[0]
    print("cube shape",cube.size)
    dx_v = 1./N
    xyz = np.stack(np.mgrid[0:1:dx_v,0:1:dx_v,0:1:dx_v])
    return xyz, dx_v

def make_phi_theta(xyz_p,center,projax):
    """ presently not very general"""
    center.shape=(3,1,1,1)
    #xyz_p = xyz-center
    #_r = np.sqrt((xyz_p**2).sum(axis=0))
    r = np.sqrt(xyz_p[0]**2 + xyz_p[1]**2 + xyz_p[2]**2)
    theta = np.arccos( xyz_p[2]/r )
    phi   = np.arctan2( xyz_p[1], xyz_p[0])
    #the x,y,z axes in healpix form a spherical axis system with r along z', phi along y', theta along x'
    #this looks the same whether viewing from inside out or outside in 
    #get theta,phi angles from projax 
    r_proj = np.sqrt(projax[0]**2 + projax[1]**2 + projax[2]**2)
    projax /= r_proj
    theta = np.arccos(projax[2]) 
    phi = np.arctan2(projax[1],projax[0])
    #now get the rotated coordinate axes note that r hat is unit magnitude by definition here since projax is normalized  
    z_p = projax
    y_p = [-1*np.sin(phi),np.cos(phi),0] 
    x_p = [np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -1*np.sin(theta)]
    #now for every point, we can find the x,y,z using new axes 
    #xyz_p=np.array([1,0,0])
    #x_new = xyz_p[0]*x_p[0] + xyz_p[1]*y_p[0] + xyz_p[2]*z_p[0]
    #y_new = xyz_p[0]*x_p[1] + xyz_p[1]*y_p[1] + xyz_p[2]*z_p[1]
    #z_new = xyz_p[0]*x_p[2] + xyz_p[1]*y_p[2] + xyz_p[2]*z_p[2]
    x_new = xyz_p[0]*x_p[0] + xyz_p[1]*x_p[1] + xyz_p[2]*x_p[2]
    y_new = xyz_p[0]*y_p[0] + xyz_p[1]*y_p[1] + xyz_p[2]*y_p[2]
    z_new = xyz_p[0]*z_p[0] + xyz_p[1]*z_p[1] + xyz_p[2]*z_p[2]
    r_new = np.sqrt(x_new**2+y_new**2+z_new**2)
    import pdb
    #theta_new = np.arccos(z_new/r_new) right but not what we want
    #phi_new = np.arctan2(y_new,x_new)
    #theta_new = np.arctan2(x_new,np.abs(z_new))
    #dcc change here
    theta_new = np.arctan2(x_new,np.sqrt(z_new**2+y_new**2))
    phi_new = np.arctan2(y_new,z_new)
    xyz_new = [x_new,y_new,z_new]
    if 0:
        pdb.set_trace()
    if 0:
        projax=0
        plt.clf()
        plt.imshow(theta_new.mean(axis=projax).T)
        plt.colorbar()
        plt.savefig("%s/plots/theta_new.png"%os.environ["HOME"])
        plt.clf()
        plt.imshow(phi_new.mean(axis=projax).T)
        plt.colorbar()
        plt.savefig("%s/plots/phi_new.png"%os.environ["HOME"])
        plt.clf()
        plt.imshow(y_new.mean(axis=projax).T)
        plt.colorbar()
        plt.savefig("%s/plots/y_new.png"%os.environ["HOME"])
        plt.clf()
        plt.imshow(x_new.mean(axis=projax).T)
        plt.colorbar()
        plt.savefig("%s/plots/x_new.png"%os.environ["HOME"])
        plt.clf()
        plt.imshow(z_new.mean(axis=projax).T)
        plt.colorbar()
        plt.savefig("%s/plots/z_new.png"%os.environ["HOME"])
        plt.clf()
        rho_theta = np.sqrt(xyz_new[0]**2 + xyz_new[2]**2)
        plt.imshow(rho_theta.mean(axis=projax).T)
        plt.colorbar()
        plt.savefig("%s/plots/rho_theta.png"%os.environ["HOME"])
        plt.clf()
        rho_phi = np.sqrt(xyz_new[1]**2 + xyz_new[2]**2)
        plt.imshow(rho_phi.mean(axis=projax).T)
        plt.colorbar()
        plt.savefig("%s/plots/rho_phi.png"%os.environ["HOME"])
    return xyz_new, phi_new, theta_new
      

    

def rotate_test(projax): 
    new_z_ax = projax /np.sqrt(projax[0]**2 + projax[1]**2 + projax[2]**2)
    rot_z = [-1*projax[1], projax[0], 0] #rotation axis
    #Normalize rotation axis
    u = rot_z / np.sqrt(rot_z[0]**2 + rot_z[1]**2 + rot_z[2]**2)
    #Now get the angle of the rotation from the dot product between the projax and the original z axis
    rotang = np.arccos(projax[2] / np.sqrt(projax[0]**2 + projax[1]**2 + projax[2]**2))
    #Now the result of the rotation for the new x and y axes comes from a rotation matrix
    rot_x_x = np.cos(rotang) + u[0]**2*(1 - np.cos(rotang))
    rot_x_y = u[0]*u[1]*(1 - np.cos(rotang))
    rot_x_z = -1*u[1] * np.sin(rotang)
    rot_x = [rot_x_x, rot_x_y, rot_x_z]
    rot_x /= np.sqrt(rot_x[0]**2 + rot_x[1]**2 + rot_x[2]**2) #I think this is redundant 
    print(rot_x[0]**2 + rot_x[1]**2 + rot_x[2]**2)
    rot_y_x = u[0]*u[1]*(1 - np.cos(rotang))
    rot_y_y = np.cos(rotang) + u[1]**2*(1 - np.cos(rotang))
    rot_y_z = u[0]*np.sin(rotang)
    rot_y = [rot_y_x, rot_y_y, rot_y_z]
    rot_y /= np.sqrt(rot_y[0]**2 + rot_y[1]**2 + rot_y[2]**2)
    print(rot_y[0]**2 + rot_y[1]**2 + rot_y[2]**2)
    return rot_x, rot_y, new_z_ax


class s1p():
    def __init__(self,ad,field, center=nar([0.5,-0.2,0.5]), projax=[0,1,0], verbose=False, Nbins=None):
        self.ad = ad
        self.field = field
        self.proj(center, projax, verbose, Nbins)
        
    def proj(self, center=nar([0.5,-0.2,0.5]),projax=[0,0,1], verbose=False, Nbins=None):
        ad = self.ad
        field = self.field
        if Nbins is None:
            Nbins = ad["x"].shape()
        xyz, dx_v = [ad["x"].in_units("code_length"), ad["y"].in_units("code_length"), ad["z"].in_units("code_length")], [ad["dx"].in_units("code_length"), ad["dy"].in_units("code_length"), ad["dz"].in_units("code_length")]
        xyz[0] = xyz[0] - center[0]*ad["x"].in_units("code_length").uq
        xyz[1] = xyz[1] - center[1]*ad["x"].in_units("code_length").uq
        xyz[2] = xyz[2] - center[2]*ad["x"].in_units("code_length").uq
        import pdb
        print("xyz mins", xyz[0].min(), xyz[1].min(), xyz[2].min())
        print("xyz maxs", xyz[0].max(), xyz[1].max(), xyz[2].max())
        #center works, projax is a dummy variable.  Presently hard coded along y.
        xyz_new,phi,theta = make_phi_theta(xyz,center,projax)
        ok3 = np.sqrt( xyz_new[0]**2+xyz_new[1]**2+xyz_new[2]**2) > .005
        print("xyz shape", np.shape(xyz_new))
        print(ok3.shape)
        print((~ok3).sum())
        got=False
        phi_collector = [] 
        collector = []
        theta_collector = []
        shift_array=[-0.5, 0.5]
        for sX in shift_array:
            for sY in shift_array:
                for sZ in shift_array:
                    #print(sX,sY,sZ)
                    #shift = np.stack([sX,sY,sZ])*dx_v
                    #shift.shape=(3,1,1,1)
                    shift = [sX*dx_v[0], sY*dx_v[1], sZ*dx_v[2]]
                    print(sX, sY, sZ)
                    xyz_shift = [xyz[0] + shift[0], xyz[1] + shift[1], xyz[2] + shift[2]]
                    xyz_new_shift, this_phi, this_theta = make_phi_theta(xyz_shift, center, projax)
                    collector.append(xyz_new_shift)
                    theta_collector.append(this_theta)
                    phi_collector.append(this_phi)
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
        print(np.shape(min_phi))
        #2.) Quantize the phi and theta bin arrays
        #the eps is to make sure the last point isn't in the Nbins+1 bin.
        eps = 1e-15
        deltaphi=max_phi - min_phi
        deltatheta = max_theta - min_theta
        ok1 = (deltaphi > np.pi/2)#*(deltatheta < np.pi/2)
        print("ok1 shape", ok1.shape)
        #dcc no don't do this.
        #ok3 = okph*okth
        #ok3 = np.ones_like(max_phi,dtype='bool')
        import pdb 
        #pdb.set_trace()
        max_phi[ok1] -= 2*np.pi
        #smallest and largest angles in the bin
        print("min_theta shape", min_theta.shape)
        minmin_theta = min_theta[ok3].min()-eps
        minmin_phi  = min_phi[ok3].min() - eps
        maxmax_theta = max_theta[ok3].max() +eps
        maxmax_phi = max_phi[ok3].max() + eps
        print("theta min/max" , minmin_theta + np.pi/2, maxmax_theta+np.pi/2)
        print("phi min/max", minmin_phi + np.pi, maxmax_phi + np.pi)
        #do i need to do this to get the right width of each bin?
        #phi_diff = np.minimum(np.abs(maxmax_phi - minmin_phi),np.abs(-2*np.pi + (maxmax_phi - minmin_phi)))
        dtheta = (maxmax_theta-minmin_theta)/(Nbins)
        dphi   = (maxmax_phi - minmin_phi)/(Nbins)
        #dphi   = (phi_diff)/(Nbins)
        #get the min and max bin for theta and phi
        deltaphi1 = max_phi - min_phi
        deltaphi2 = (max_phi - 2*np.pi) - min_phi
        actdeltaphi = np.minimum(np.abs(deltaphi1),np.abs(deltaphi2))
        import pdb 
        min_theta_bin =((min_theta - minmin_theta)//dtheta).astype('int')
        min_phi_bin =  ((min_phi - minmin_phi)//dphi).astype('int')
        #min_phi_bin =  (np.minimum((min_phi - minmin_phi), 2*np.pi - (min_phi - minmin_phi))//dphi).astype('int')
        max_theta_bin =((max_theta - minmin_theta)//dtheta).astype('int')
        max_phi_bin =  ((max_phi - minmin_phi)//dphi).astype('int')
        #max_phi_bin =  (np.minimum(max_phi - minmin_phi, 2*np.pi - (max_phi - minmin_phi))//dphi).astype('int')
        #number of populated bins in theta and phi
        Ntheta = max_theta_bin-min_theta_bin
        import pdb
        Nphi   = max_phi_bin - min_phi_bin #dcc maybe needs a +1? No, a max.
        ones = np.ones_like(Ntheta)
        Ntheta = np.maximum(Ntheta, ones)
        Nphi = np.maximum(Nphi, ones)
        Nphi[~ok3] = 0
        Ntheta[~ok3] = 0
        Nphi_flat = Nphi.flatten()
        max_phi_bin_flat = max_phi_bin.flatten()
        min_phi_bin_flat = min_phi_bin.flatten()
        max_phi_flat = max_phi.flatten()
        min_phi_flat = min_phi.flatten()
        import pdb
        plt.clf()
        plt.plot(Nphi_flat)
        plt.savefig("Nphiflaty.png")
        plt.close()
        plt.plot(max_phi_flat)
        plt.savefig("maxphiflaty.png")
        plt.close()
        plt.plot(min_phi_flat)
        plt.savefig("minphiflaty.png")
        ind = np.argmax(Nphi_flat)
        if 0: 
            projax = 2
            plt.clf() 
            plt.imshow(Ntheta.max(axis=projax))
            plt.colorbar()
            plt.savefig("%s/plots/Nthetaproj.png"%os.environ["HOME"])
        if 0:
            plt.clf()
            plt.plot(xyz[0][ok3],Ntheta[ok3])
            plt.savefig("%s/plots/xNtheta.png"%os.environ["HOME"])
            plt.clf()
            plt.plot(xyz[1][ok3],Ntheta[ok3])
            plt.savefig("%s/plots/yNtheta.png"%os.environ["HOME"])
            plt.clf()
            plt.plot(xyz[2][ok3],Ntheta[ok3])
            plt.savefig("%s/plots/zNtheta.png"%os.environ["HOME"])
            plt.clf()
            plt.plot(xyz[0][ok3],Nphi[ok3])
            plt.savefig("%s/plots/xNphi.png"%os.environ["HOME"])
            plt.clf()
            plt.plot(xyz[1][ok3],Nphi[ok3])
            plt.savefig("%s/plots/yNphi.png"%os.environ["HOME"])
            plt.clf()
            plt.plot(xyz[2][ok3],Nphi[ok3])
            plt.savefig("%s/plots/zNphi.png"%os.environ["HOME"])
        ind_min = np.argmin(Nphi_flat)
        Ntheta_max = Ntheta.max()
        Nphi_max = Nphi.max() #why is this so big?
        Ntheta_bins=Nbins
        Nphi_bins = Nbins
        #pdb.set_trace()

        #set up target array H and mask array F
        #and theta and phi bins
        H = np.zeros(Ntheta_bins*Nphi_bins)
        import yt
        H = yt.YTArray(H, field.units)
        F = np.zeros(Ntheta_bins*Nphi_bins)
        F[:]=np.nan
        theta_bins = np.linspace( minmin_theta, maxmax_theta, Ntheta_bins)
        phi_bins = np.linspace( minmin_phi, maxmax_phi, Nphi_bins)
        theta_cen = theta_bins + 0.5*(theta_bins[1]-theta_bins[0])
        phi_cen = phi_bins + 0.5*(phi_bins[1]-phi_bins[0])
        coordPhi, coordTheta = np.meshgrid(phi_cen,theta_cen, indexing='ij')
        coordPhi = yt.YTArray(coordPhi, 'dimensionless')
        coordTheta = yt.YTArray(coordTheta, 'dimensionless')
        import healpy as hp 
        #coordPhi, coordTheta = hp.pix2ang(Nbins, np.arange(hp.nside2npix(Nbins)))
        import loop
        import pdb
        theta_middle = (max_theta + min_theta) / 2  #midpoints of angles in each cell for determining healpix indices 
        phi_middle = (max_phi +  min_phi) / 2      #in loop.big_loop
        theta_middle += np.pi/2 
        phi_middle += np.pi
        H,F,m = loop.big_loop(H,F,field,Nphi_max,Ntheta_max,min_phi_bin,min_theta_bin,Ntheta_bins,Nphi,Ntheta,theta_middle, phi_middle,verbose=True)
        print(np.shape(H))
        self.map = m
        #everything from here on isn't used for healpix output and currently doesn't produce anything meaningful
        H.shape = Nphi_bins,Ntheta_bins
        F.shape = Nphi_bins,Ntheta_bins

        self.coordPhi=coordPhi
        self.coordTheta=coordTheta
        self.H = H
        self.mask = F
        #pdb.set_trace()
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
    import pdb 
    p=ax0.pcolormesh(coordPhi[::-1], coordTheta, H, cmap=cmap,norm=nrm )
    print(np.shape(p))
    fig.colorbar(p,ax=ax0)
    
    fig.savefig('%s/test22'%plot_dir)

def E_B_maps(Q_map, U_map): 
    I_map = np.sqrt(Q_map*Q_map + U_map*U_map)
    import healpy as hp 
    alms = hp.sphtfunc.map2alm([I_map, Q_map, U_map])
    nside = hp.get_nside(I_map)
    E_map = hp.sphtfunc.alm2map(alms[1],nside)
    B_map = hp.sphtfunc.alm2map(alms[2],nside)
    plt.clf()
    hp.mollview(E_map)
    hp.graticule()
    plt.savefig("Emap.png")
    plt.clf()
    hp.mollview(B_map)
    hp.graticule()
    plt.savefig("Bmap.png")
    return E_map, B_map

def rotationtest(): 
    #testing rotation for any axis projection
    rot_x, rot_y, rot_z = rotate_test([0, 1, 0])
    x = [1, 0, 0]
    y = [0, 1, 0]
    z = [0, 0, 1] 
    ax = plt.figure().add_subplot(projection='3d')
    ax.quiver([0,0,0],[0,0,0],[0,0,0],x,y,z, normalize=True,length = .1)
    ax.quiver([0,0,0],[0,0,0],[0,0,0],rot_x,rot_y,rot_z, normalize=True, length = .1)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.savefig('rotationtest.png')
    plt.clf()
    ax = plt.figure().add_subplot(projection='3d')
#ax.quiver([0,0,0],[0,0,0],[0,0,0],x,y,z, normalize=True,length = .1)
    ax.quiver([0,0,0],[0,0,0],[0,0,0],rot_x,rot_y,rot_z, normalize=True, length = .1, color = ['red','blue','green'])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.savefig('rotationtestfinal.png')
    plt.clf()
    ax = plt.figure().add_subplot(projection='3d')
    ax.quiver([0,0,0],[0,0,0],[0,0,0],x,y,z, normalize=True,length = .1, color = ['red','blue','green'])
#ax.quiver([0,0,0],[0,0,0],[0,0,0],rot_x,rot_y,rot_z, normalize=True, length = .1, color = 'red')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.savefig('rotationtestinit.png')

def angletest():
    import projector.proj as proj
    N = 128 #set cube edge
    cube = proj.make_cube(128) #call to function that makes the cyclinders for proj test 
    xyz, dx_v = make_xyz(cube) #get array for coordinates and the dx between each cell
    r, theta, phi = make_phi_theta(xyz, nar([0.5,-0.2,0.5]),nar([1,1,0])) #get r, theta, phi from general rotation about given axis
    #now need to make the image but I want the angles to be plotted instead of something like density
    center = nar([0.5,0.2,0.5])
    projax = nar([1,1,0])
    #I think I might be able to do this in maker.py instead

def plot_xyz(): 
    import projector.proj as proj 
    N = 128 
    cube = proj.make_cube(128)
    xyz,dx_v = make_xyz(cube)
    r, theta, phi, xyz = make_phi_theta(xyz, nar([0.5,-0.2,0.5]), nar([0.,1.,0.]))
    plt.close('all')
    im = np.mean(xyz[0],axis=0)
    plt.imshow(im, origin='lower')
    plt.colorbar()
    plt.savefig('x.png')
    plt.close('all')
    im = np.mean(xyz[1],axis=1)
    plt.imshow(im, origin='lower')
    plt.colorbar()
    plt.savefig('y.png')
    plt.close('all')
    im = np.mean(xyz[2],axis=2)
    plt.imshow(im, origin='lower')
    plt.colorbar()
    plt.savefig('z.png')
    plt.close('all')
