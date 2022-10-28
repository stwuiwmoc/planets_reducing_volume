# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 23:02:24 2021

@author: swimc

円形にmaskしたときにrms = 1 um になるような proper へのrmsの入力値を最小二乗で探す
"""
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import scipy as sp
import os
import time
import proper as pr

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import mpl_toolkits.axes_grid1

def mkfolder(suffix = ""):
    """    
    Parameters
    ----------
    suffix : str, optional
        The default is "".

    Returns
    -------
    str ( script name + suffix )
    """
    name = os.path.basename(__file__)
    name = name.replace(".py", "") + suffix
    os.makedirs(name, exist_ok=True)
    return name


def zer_raw(zer_order, zer_num, zer_rms):
    wavelength = 500*10**-9
    diam = m1_radi * 1000 * 2 # m 単位の直径に直す

    zer_poly_arr = np.zeros(zer_order)
    zer_poly_arr[zer_num-1] = zer_rms
    
    wavestruct = pr.prop_begin(diam, wavelength, px, 1)
    wfe = pr.prop_zernikes(wavestruct, [zer_num], [zer_poly_arr[zer_num-1]])
    return wfe

def rms_micron(array, digits):
    #input is [mm], output is [um]
    data = array*1000
    sigma = np.nansum(data**2)
    num = np.count_nonzero(np.logical_not(np.isnan(array)))
    rms = np.sqrt(sigma / num)
    rms_str = str(round(rms, digits))
    return rms_str, rms

def zer_func(X):
    zer_rms = X * 10 ** -6
    zer = tf * zer_raw(zer_order, zer_num, zer_rms) * 1000
    rms = rms_micron(zer, 2)
    f = (rms[1] - 2)**2
    return f
    
def image_plot(fig, title, position, c, c_scale, cb_micron=True):
    cmap = cm.jet
    fs = 15

    ax = fig.add_subplot(position)
    ax.imshow(c, interpolation="nearest", cmap=cmap, vmin=np.nanmin(c_scale), vmax=np.nanmax(c_scale), origin="lower", extent=[-925, 925, -925, 925])
    ax.set_title(title, fontsize=fs)
  
    divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes('right', '5%', pad='3%')
    
    if cb_micron==True:
        norm = Normalize(vmin=np.nanmin(c_scale*1000), vmax=np.nanmax(c_scale*1000)) # c_scaleがmm表記だと見にくいのでμ表記に変更
        cbar_title = r"[$\mu$m]"
    else:
        norm = Normalize(vmin=np.nanmin(c_scale), vmax=np.nanmax(c_scale))
        cbar_title = "[mm]"
    mappable = cm.ScalarMappable(norm = norm, cmap = cm.jet)
    cbar = fig.colorbar(mappable, ax=ax, cax=cax)
    cbar.set_label(cbar_title, fontsize=fs)
        
    return ax

if __name__ == '__main__':
    
    # parametar --------------------------------------------------------------
    m1_radi = 1850/2
    
    zer_order, zer_num = 21, 4
    init = 1
    #-------------------------------------------------------------------------    
    px = 1023
    
    xx, yy = np.meshgrid(np.linspace(-m1_radi, m1_radi, px),np.linspace(-m1_radi, m1_radi, px))

    tf = np.where(xx**2+yy**2<m1_radi**2, True, np.nan)
    
    rms_arr = np.empty(zer_order)
    start_time = time.time()
    
    
    for i in range(1, zer_order+1):
        zer_num = i    
        print("\n",i,"\n")
        result = minimize(fun=zer_func, x0=init, method="Powell")
        
        rms_m = result["x"] * 10** -6
        rms_arr[i-1] = result["x"][0]        
    
        
        zer_m = tf * zer_raw(zer_order, zer_num, rms_m) * 1000
                
        fig = plt.figure()
        title = str(zer_num)+"\n"\
            + "Input RMS = " + str(result["x"][0]) + r" [$\mu$m]" + "\n"\
                + "Output RMS = " + rms_micron(zer_m, 2)[0] + r" [$\mu$m]"
        ax = image_plot(fig, title, 111, zer_m, zer_m)
        
        fig.tight_layout()
        picname = mkfolder() + "/z" + str(zer_num) + ".png"
        fig.savefig(picname)
        #fig.clf()
    
    
        end_time = time.time()
        print(end_time - start_time)
    
    print(rms_arr)
    for i in rms_arr:
        print(i, end=", ")