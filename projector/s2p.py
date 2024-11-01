
from starter2 import *
from scipy.ndimage import gaussian_filter
import pdb
import projector.proj as proj


def project(cube, xyz, dxyz, proj_center, proj_axis):

    
    Nz = cube.size
    proj_center.shape=3,1
    xyz= xyz - proj_center
    if 1:
        #this shift puts them in an order that makes sense to draw lines
        shifter = np.array([[[-0.5, -0.5, +0.5, +0.5,  0.5, +0.5, -0.5, -0.5]],
                            [[-0.5, -0.5, -0.5, -0.5,  0.5, +0.5, +0.5, +0.5]],
                            [[-0.5,  0.5, +0.5, -0.5, -0.5, +0.5, +0.5, -0.5]]])

    if 0:
        #the result from the programatic build
        shifter = np.array([[[-0.5, -0.5, -0.5, -0.5,  0.5,  0.5,  0.5,  0.5]],
                            [[-0.5, -0.5,  0.5,  0.5, -0.5, -0.5,  0.5,  0.5]],
                            [[-0.5,  0.5, -0.5,  0.5, -0.5,  0.5, -0.5,  0.5]]])
    if 0:
        #build
        shift_array=[-0.5, 0.5]
        shifter = None
        xyz = xyz-proj_center
        for sX in shift_array:
            for sY in shift_array:
                for sZ in shift_array:
                    this_shift = np.array([sX,sY,sZ])
                    this_shift.shape = this_shift.size,1
                    if shifter is None:
                        shifter = this_shift
                    else:
                        shifter = np.hstack([shifter, this_shift])
        shifter.shape = 3,1,8
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
    #ec_theta = np.take_along_axis(theta,asrt,axis=1)[...,2:]
    #ec_phi   = np.take_along_axis(phi,asrt,axis=1)[...,2:]
    plt.clf()
    #plt.scatter(phi, theta, c='r')
    #plt.scatter(corners_new[0].flatten(),corners_new[1].flatten(),c='r')
    #plt.scatter(corners[0].flatten(),corners[1].flatten(),c='g')
    plt.scatter(phi_new,theta_new)
    plt.plot(phi_new[0,:], theta_new[0,:])
    plt.savefig('%s/theta_phi'%plot_dir)
    plt.clf()
    plt.scatter(corners_new[0].flatten(),corners_new[1].flatten(),c='r')
    plt.plot(corners_new[0,0,:], corners_new[1,0,:])
    plt.savefig('%s/corners'%plot_dir)



    #pdb.set_trace()
    return corners_new

