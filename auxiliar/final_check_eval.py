from astropy.io import fits
from astropy.modeling import Fittable2DModel
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from matplotlib.pyplot import matshow
import scipy.signal
from Functions import dist, twoD_Gaussian, moments

hdulist = fits.open('estrela_cut.fits')

hdulist.info()

image_data = hdulist[0].data
image_data2 = np.nan_to_num(image_data)

xx = np.linspace(0, 482, 483)
yy = np.linspace(0, 449, 450)
xx, yy = np.meshgrid(xx, yy)

moments(image_data2)

height = moments(image_data2)[0] 
x = moments(image_data2)[1]
y = moments(image_data2)[2]
width_x = moments(image_data2)[3]
width_y = moments(image_data2)[4]

initial_guess = (height,x,y,width_x,width_y,0,140)
popt, pcov = optimize.curve_fit(twoD_Gaussian, (xx, yy), image_data2.ravel(), p0=initial_guess)

Xin, Yin = np.mgrid[0:450, 0:483]
r_r = dist(Xin,Yin,popt[1],popt[2],1,1)
r_i = width_x
r_m = np.sqrt(2)*r_i
r_e = np.sqrt(3)*r_i

for index, val in enumerate(np.arange(0,r_r.max())):
    mask = np.bitwise_and(r_r > r_m, r_r <= r_e)
    mask2= np.bitwise_and(r_r > 0, r_r < r_i)
    sky = np.median(image_data2[mask])
    image_data2[mask2]=sky
	
plt.imshow(image_data2)
plt.show()

hdu = fits.PrimaryHDU(image_data2)
hdu1 = fits.HDUList([hdu])
hdu1.writeto('subtracted_region.fits')
