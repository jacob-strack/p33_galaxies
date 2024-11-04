import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

# Set HEALPix parameters
nside = 16  # HEALPix resolution parameter
ipix = 42   # Pixel index to plot boundaries for

# Get the vertices of the pixel boundaries in Cartesian coordinates
xyz = hp.boundaries(nside, ipix, step=1)

# Convert Cartesian coordinates to spherical coordinates
theta, phi = hp.vec2ang(xyz.T)

# Plot a blank map using Mollweide projection
hp.mollview(np.zeros(hp.nside2npix(nside)), title="HEALPix Pixel Boundary", cmap="Blues")

# Overlay the boundary line of the pixel
hp.projplot(theta, phi, 'r-', linewidth=2, label=f'Pixel {ipix} Boundary')

# Add legend
plt.legend(loc="upper left")

plt.savefig('%s/stupid4'%plot_dir)

