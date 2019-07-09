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
