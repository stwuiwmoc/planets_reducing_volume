# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 23:18:45 2021

@author: swimc
"""

import numpy as np
import scipy as sp
import sklearn as sl
import time

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import mpl_toolkits.axes_grid1

def cz1(p, cx, cy, aqx, phi):        
    a = 1/(2*p)
    aqz = a * ( aqx**2 + 0**2 )
    theta = np.arctan(2*a*aqx)
    
    D = -4*a**2*cx**2*np.sin(phi)**2*np.sin(theta)**2 - 8*a**2*cx*cy*np.sin(phi)*np.sin(theta)**2*np.cos(phi) + 4*a**2*cy**2*np.sin(phi)**2*np.sin(theta)**2 - 4*a**2*cy**2*np.sin(theta)**2 + 2*a*aqx*np.sin(2*theta) + 4*a*aqz*np.sin(theta)**2 + 4*a*cx*np.sin(theta)*np.cos(phi) - 4*a*cy*np.sin(phi)*np.sin(theta) - np.sin(theta)**2 + 1

    cz1 = (4*a*aqx*np.sin(theta) - a*cx*(np.sin(phi - 2*theta) - np.sin(phi + 2*theta)) - a*cy*(np.cos(phi - 2*theta) - np.cos(phi + 2*theta)) - 2*np.sqrt(D) + 2*np.cos(theta))/(4*a*np.sin(theta)**2)
    return cz1

def func(t, o, p, q):
    x = t[0]
    y = t[1]
    return cz1(p, x, y, q, o)

if __name__ == '__main__':
    
    # parametar --------------------------------------------------------------
    # p 曲率半径, q 軸外し距離, o 回転角phi（np.pi不要）
    # a = (a0, da, a_num) -> a0 変化幅の中心, da 中心からどのくらい幅をとるか, a_num 何分割か
    """
    o = (0,)
    p = (8666, 1, 11)
    q = (1800, 1, 11)
    """
    o0 = 0
    p0 = 8666
    q0 = 1800
    
    o_t = 0 + np.pi/2
    p_t = 8666 - 2
    q_t = 1800 + 1
    
    #-------------------------------------------------------------------------
    
    px = 1023
    m1_radi = 1850/2
    
    x_arr = np.linspace(-m1_radi, m1_radi, px)
    y_arr = np.linspace(-m1_radi, m1_radi, px)
    
    xx, yy = np.meshgrid(x_arr, y_arr)
    
    t = np.array([xx.flatten(), yy.flatten()])
    
    z_t = cz1(p_t, xx, yy, q_t, o_t)
    
    z_t_1d = np.squeeze(z_t.reshape((px**2, 1)))
    
    param, covariance = sp.optimize.curve_fit(func, t, z_t_1d, p0=[o0, p0, q0])
    
    
    