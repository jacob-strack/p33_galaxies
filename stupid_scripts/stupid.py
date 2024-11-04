import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

# Set the HEALPix parameters
nside = 16  # Resolution parameter
ipix = 42   # Pixel index to plot boundaries for

# Get the vertices of the pixel boundaries in Cartesian coordinates
xyz = hp.boundaries(nside, ipix, step=1)

# Convert Cartesian coordinates to spherical coordinates
theta, phi = hp.vec2ang(xyz.T)

# Convert theta, phi to RA and Dec in degrees for plotting
ra = np.degrees(phi)  # Right ascension
dec = 90 - np.degrees(theta)  # Declination

# Close the loop by appending the first vertex at the end
ra = np.append(ra, ra[0])
dec = np.append(dec, dec[0])

# Plot the boundaries
plt.figure(figsize=(6, 6))
plt.plot(ra, dec, marker='o', linestyle='-', color='blue')
plt.title(f'Pixel Boundaries for pixel {ipix} (nside={nside})')
plt.xlabel("Right Ascension (degrees)")
plt.ylabel("Declination (degrees)")
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig('%s/stupid'%plot_dir)

