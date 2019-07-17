# -*- coding: utf-8 -*-
#!python number=enable

from astropy.io import fits
from astropy.modeling import Fittable2DModel
import matplotlib.pyplot as plt
import numpy as np
import copy
from scipy import optimize

#Calling image file
hdulist = fits.open('file_name.fits')
hdulist.info()

#Excluding nan from data
image_data = hdulist[0].data
image_data2 = np.nan_to_num(image_data)

def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x, y) = xdata_tuple                                                        
    xo = float(xo)                                                              
    yo = float(yo)                                                              
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)   
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)    
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)   
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)         
                        + c*((y-yo)**2)))                                   
    return g.ravel()
     
def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y
    
#Creating a grid of zeros
xx = np.linspace(0,image_data2.shape[1] -1, image_data2.shape[1])
yy = np.linspace(0,image_data2.shape[0] -1, image_data2.shape[0])
xx, yy = np.meshgrid(xx, yy)
    
moments(image_data2)

height = moments(image_data2)[0] 
x = moments(image_data2)[1]
y = moments(image_data2)[2]
width_x = moments(image_data2)[3]
width_y = moments(image_data2)[4]    

#sky and theta used have values found in the data header, data image or outside the file
initial_guess = (height,x,y,width_x,width_y,sky,theta)
popt, pcov = optimize.curve_fit(twoD_Gaussian, (xx, yy), image_data2.ravel(), p0=initial_guess)    
    
data_fit = twoD_Gaussian((xx, yy), *popt)

fig, ax = plt.subplots(1, 1)
ax.hold(True)
ax.imshow(image_data2.reshape(image_data2.shape[0],image_data2.shape[1]), cmap=plt.cm.jet, origin='bottom',
    extent=(xx.min(), xx.max(), yy.min(), yy.max()))
ax.contour(xx, yy, data_fit.reshape(image_data2.shape[0], image_data2.shape[1]), 8, colors='w')
plt.show()
