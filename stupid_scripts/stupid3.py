import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

# Set HEALPix parameters
nside = 4  # HEALPix resolution parameter
ipix = 42   # Pixel index to plot boundaries for

# Get the vertices of the pixel boundaries in Cartesian coordinates
xyz = hp.boundaries(nside, ipix, step=1)

# Convert Cartesian coordinates to spherical coordinates
theta, phi = hp.vec2ang(xyz.T)
print(theta)
print(phi)

# Plot a blank map using Mollweide projection
hp.mollview(np.arange(hp.nside2npix(nside)), title="HEALPix Pixel Boundary", cmap="Blues")

# Overlay the boundary line of the pixel
#hp.projplot(theta, phi, 'r-', linewidth=2, label=f'Pixel {ipix} Boundary')
hp.projplot(theta, phi,'bo-')
hp.projscatter([0],[0],c='r')
hp.projscatter([0],[np.pi],c='orange')
hp.projscatter([np.pi],[np.pi],c='yellow')
hp.projscatter([0.1*np.pi],[0.1*np.pi],c='green')
hp.projscatter([0.1*np.pi],[-0.1*np.pi],c='cyan')
hp.projscatter([0.5*np.pi],[-0.5*np.pi],c='blue')
hp.projscatter([0.5*np.pi],[np.pi],c='k',marker='*')
hp.projscatter([0.5*np.pi],[0],c='k',marker='s')
hp.projscatter([0.5*np.pi],[np.pi-1e-10],c='k',marker='^')
#hp.projplot(np.array([0.35,0.36]),np.array([6.05,5.89]))
#hp.projplot([0.35,0.36],[6.05,5.89],'r-')
#hp.projplot(np.array([0,np.pi]),np.array([0,2*np.pi]),'r-')
# Add legend
#plt.legend(loc="upper left")

plt.savefig('%s/stupid3'%plot_dir)

