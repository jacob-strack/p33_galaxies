import yt
from starter2 import *
import pdb
import projector.proj as proj
from yt.fields.api import ValidateParameter
reload(proj)

#stride one projector
import projector.s1p as s1p
reload(s1p)
from scipy.ndimage import gaussian_filter

def _Q_int(field, data): 
    #new basis vectors in terms of domain coordinates 
    projax_x = data.get_field_parameter("projax_x")
    projax_y = data.get_field_parameter("projax_y")
    projax_z = data.get_field_parameter("projax_z")
    print(data.has_field_parameter("projax"))
    projax = yt.YTArray([projax_x, projax_y, projax_z], ds.unit_system["length"])
    theta = np.arccos(projax[2]) 
    phi = np.arctan2(projax[1],projax[0])
    z_p = projax
    y_p = [-1*np.sin(phi), np.cos(phi), 0]
    x_p = [np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -1*np.sin(theta)]
    B = [data["Bx"], data["By"], data["Bz"]]
    Bx_new = B[0]*x_p[0] + B[1]*x_p[1] + B[2]*x_p[2] 
    By_new = B[0]*y_p[0] + B[1]*y_p[1] + B[2]*y_p[2] 
    Bz_new = B[0]*z_p[0] + B[1]*z_p[1] + B[2]*z_p[2]
    B_sq = B[0]**2 + B[1]**2 + B[2]**2
    return data["density"]*(Bx_new**2 + By_new**2) / B_sq

def _U_int(field, data): 
    #new basis vectors in terms of domain coordinates 
    projax_x = data.get_field_parameter("projax_x")
    projax_y = data.get_field_parameter("projax_y")
    projax_z = data.get_field_parameter("projax_z")
    projax = yt.YTArray([projax_x, projax_y, projax_z], ds.unit_system["length"])
    theta = np.arccos(projax[2]) 
    phi = np.arctan2(projax[1],projax[0])
    z_p = projax
    y_p = [-1*np.sin(phi), np.cos(phi), 0]
    x_p = [np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -1*np.sin(theta)]
    B = [data["Bx"], data["By"], data["Bz"]]
    Bx_new = B[0]*x_p[0] + B[1]*x_p[1] + B[2]*x_p[2] 
    By_new = B[0]*y_p[0] + B[1]*y_p[1] + B[2]*y_p[2] 
    Bz_new = B[0]*z_p[0] + B[1]*z_p[1] + B[2]*z_p[2]
    B_sq = B[0]**2 + B[1]**2 + B[2]**2
    return data["density"]*2*Bx_new*By_new / B_sq

import healpy as hp 
N = 256
Nbins = N*4
ds = yt.load("~/scratch/DD0026/DD0026")
projax = [0.,0.,1.]
ad = ds.all_data()
ad.set_field_parameter("projax_x", projax[0])
ad.set_field_parameter("projax_y", projax[1])
ad.set_field_parameter("projax_z", projax[2])
ds.add_field(name=("gas", "Q_integrand"), function=_Q_int, sampling_type = "cell", units = "g/cm**3", validators = ValidateParameter(["projax_x", "projax_y", "projax_z"]))
ds.add_field(name=("gas","U_integrand"), function=_U_int, sampling_type = "cell", units = "g/cm**3", validators = ValidateParameter(["projax_x","projax_y","projax_z"]))
print("test", ad.get_field_parameter("projax_z"))
field = ad["U_integrand"].in_units("code_density")
#ppp= s1p.s1p(cube, center=nar([0.5,-0.2,0.5]), verbose=True,Nbins=1024)
#testing if projax works...
ppp= s1p.s1p(ad,field, center=nar([0.5,0.5,0.5]), projax=[1,0,0], verbose=True, Nbins=N)
H = gaussian_filter(ppp.H,2)
import healpy as hp 
#H = ppp.H

print("H valid value", hp.pixelfunc.maptype(ppp.map))
print("map shape", np.shape(ppp.map))
s1p.plot_image(ppp.coordPhi, ppp.coordTheta, H, "%s/test5"%plot_dir, mask=ppp.mask)
print(np.shape(H))
#put other field here when making E,B maps 
field2 = ad["Q_integrand"].in_units("code_density")
p2 = s1p.s1p(ad,field2, center=nar([0.5, 0.5, 0.495]), projax=[1,0,0], verbose=True, Nbins = N)
H2 = gaussian_filter(p2.H,2) 

print("H valid value", hp.pixelfunc.maptype(ppp.map))
print("H2 valid value", hp.pixelfunc.maptype(p2.map))
#map stuff 
E,B = s1p.E_B_maps(ppp.map,p2.map)
if 0:
    #tests
    #proj.test3()
    #proj.test5()
    proj.test4()


