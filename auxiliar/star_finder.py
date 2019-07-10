# -*- coding: utf-8 -*-

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.stats import sigma_clipped_stats  
from photutils import DAOStarFinder, CircularAperture, find_peaks
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
import numpy as np

hdulist = fits.open('ngc7412B.fits')
hdulist.info()
image_data = hdulist[0].data

mean, median, std = sigma_clipped_stats(image_data, sigma=3.0)
daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
sources = daofind(image_data - median) 
for col in sources.colnames:    
	sources[col].info.format = '%.8g'
  
positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r=4.)
norm = ImageNormalize(stretch=SqrtStretch())
plt.imshow(image_data, cmap='Greys', origin='lower', norm=norm)
apertures.plot(color='blue', lw=1.5, alpha=0.5)
plt.show()

#For a better visualization for the image:

plt.imshow(image_data, cmap='Greys', origin='lower', norm=LogNorm())
apertures.plot(color='blue', lw=1.5, alpha=0.5)
plt.show()
