# -*- coding: utf-8 -*-
from astropy.io import fits
from astropy.modeling import Fittable2DModel
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
From Func import moments, fitgaussian

#Calling the fits file
hdulist = fits.open('file_name.fits')
hdulist.info()

#Excluding the nan data
image_data = hdulist[0].data
image_data2 = np.nan_to_num(image_data)
   
#Creating the data and plotting the gaussian fit

moments(image_data2)

plt.matshow(image_data2, cmap=plt.cm.gist_earth_r)

params = fitgaussian(image_data2)
fit = gaussian(*params)

plt.contour(fit(*np.indices(image_data2.shape)), cmap=plt.cm.copper)
ax = plt.gca()
(height, x, y, width_x, width_y) = params

plt.text(0.95, 0.05, """
x : %.1f
y : %.1f
width_x : %.1f
width_y : %.1f""" %(x, y, width_x, width_y),
        fontsize=16, horizontalalignment='right',
        verticalalignment='bottom', transform=ax.transAxes)
plt.show()

hdulist.close()
