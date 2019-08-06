# -*- coding: utf-8 -*-

import numpy as np 
from scipy.interpolate import splrep,splev

def mat_theta(xg,yg,theta):
    """
    Method to calculate the angles's matrix of the galaxy
    """
    mask = np.bitwise_and(xg >= 0, yg >= 0 )
    theta[mask] = np.abs(np.arctan(yg[mask]/xg[mask]))
    mask = np.bitwise_and(xg <= 0, yg >= 0 )
    theta[mask] = np.pi - np.abs(np.arctan(yg[mask]/xg[mask]))
    mask = np.bitwise_and(xg <= 0, yg <= 0 )
    theta[mask] = np.pi + np.abs(np.arctan(yg[mask]/xg[mask]))
    mask = np.bitwise_and(xg >= 0, yg <= 0 )
    theta[mask] = 2*np.pi - np.abs(np.arctan(yg[mask]/xg[mask]))
    return theta
    
    
def rot_ceu(theta_a,phi_a):
    """
    Method to calculate the rotation in the sky plane starting from a angle phi in radians 
    """
    phi_0 = np.pi*3./2.
    ang = phi_0 - phi_a
    theta_a = theta_a + ang
    mask = np.where(theta_a>2*np.pi)
    theta_a[mask] = theta_a[mask] - 2*np.pi
    return theta_a

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
