
from starter2 import *

import projector.proj as proj
reload(proj)

#stride one projector
import projector.s1p as s1p
reload(s1p)
from scipy.ndimage import gaussian_filter

N = 128
Nbins = N*4
cube = proj.make_cube(N)

#ppp= s1p.s1p(cube, center=nar([0.5,-0.2,0.5]), verbose=True,Nbins=1024)
ppp= s1p.s1p(cube, center=nar([0.5,-0.2,0.5]), verbose=True, Nbins=256)
H = gaussian_filter(ppp.H,2)
#H = ppp.H

s1p.plot_image(ppp.coordPhi, ppp.coordTheta, H, "%s/test5"%plot_dir, mask=ppp.mask)

if 0:
    #tests
    #proj.test3()
    #proj.test5()
    proj.test4()
