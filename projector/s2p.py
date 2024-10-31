
from starter2 import *
from scipy.ndimage import gaussian_filter
import pdb
import projector.proj as proj


def project(cube, xyz, dxyz, proj_center, proj_axis):

    
    Nz = cube.size
    shift_array=[-0.5, 0.5]
    shifter = None
    proj_center.shape=3,1
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
    pdb.set_trace()

    return corners_new

