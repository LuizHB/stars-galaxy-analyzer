# -*- coding: utf-8 -*-
#!python number=enable

from astropy.io import fits
import numpy as np
from sgolay2D import sgolay2d
from matplotlib.pyplot import matshow

hdulist = fits.open('file_name.fits')

hdulist.info()

image_data = hdulist[0].data
image_data2 = np.nan_to_num(image_data)

#smooth
data_smooth = sgolay2d(image_data2, window_size=29, order=4)

#plot smooth
matshow(data_smooth)

hdulist.close()
