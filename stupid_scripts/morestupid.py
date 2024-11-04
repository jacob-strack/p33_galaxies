import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

# Set HEALPix parameters
nside = 32  # HEALPix resolution parameter
npix = hp.nside2npix(nside)  # Total number of pixels

# Create a sample map with random values
m = np.random.rand(npix)

# Select a pixel to highlight
ipix = [200,10]  # Arbitrary pixel index to mark

# Convert the pixel index to spherical coordinates (theta, phi)
theta, phi = hp.pix2ang(nside, ipix)

# Plot the map using Mollweide projection
hp.mollview(m, title="Random HEALPix Map with Highlighted Pixel", cmap="viridis")

# Overlay the pixel using projplot
hp.projplot(theta, phi, 'ro', markersize=10, label=f'Pixel {ipix}',linewidth=3)
plt.legend(loc="upper left")

plt.savefig('%s/morestupid'%plot_dir)
