# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 21:33:15 2021

@author: swimc

rms 1 um のzernike各項で与えられる鏡面誤差に対して研磨量を最小化する
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
    diam = m1_radi / 1000 * 2 # m 単位の直径に直す

    zer_poly_arr = np.zeros(zer_order)
    zer_poly_arr[zer_num-1] = zer_rms
    
    wavestruct = pr.prop_begin(diam, wavelength, px, 1)
    wfe = pr.prop_zernikes(wavestruct, [zer_num], [zer_poly_arr[zer_num-1]])
    return wfe
  

def cz1(o, p, aqx, cx, cy):
    """
    順番が o, p, q に代わっているので注意

    Parameters
    ----------
    o : float [mrad]
        inputは [mrad]
    p : float [mm]
        曲率半径
    aqx : floot [mm]
        軸外し距離 (off-axis)
    cx : 2D-mesh-array [mm]
    cy : 2D-mesh-array [mm]

    Returns
    -------
    cz1 : 2D-mesh-array
        input された o, p, q での切り取り放物面

    """
    phi = o / 1000 # [mrad] -> [rad] への換算
    
    a = 1/(2*p)
    aqz = a * ( aqx**2 + 0**2 )
    theta = np.arctan(2*a*aqx)
    
    D = -4*a**2*cx**2*np.sin(phi)**2*np.sin(theta)**2 - 8*a**2*cx*cy*np.sin(phi)*np.sin(theta)**2*np.cos(phi) + 4*a**2*cy**2*np.sin(phi)**2*np.sin(theta)**2 - 4*a**2*cy**2*np.sin(theta)**2 + 2*a*aqx*np.sin(2*theta) + 4*a*aqz*np.sin(theta)**2 + 4*a*cx*np.sin(theta)*np.cos(phi) - 4*a*cy*np.sin(phi)*np.sin(theta) - np.sin(theta)**2 + 1

    cz1 = (4*a*aqx*np.sin(theta) - a*cx*(np.sin(phi - 2*theta) - np.sin(phi + 2*theta)) - a*cy*(np.cos(phi - 2*theta) - np.cos(phi + 2*theta)) - 2*np.sqrt(D) + 2*np.cos(theta))/(4*a*np.sin(theta)**2)
    return cz1


def volume_func(X):
    o, p, q = X
    para_t = np.where(tf==True, cz1(o, p, q, xx, yy), np.nan)
    
    para_0 = np.where(tf==True, cz1(o0, p0, q0, xx, yy), np.nan)
    surf = para_0 + raw_for_minimize
    
    err_t = surf - para_t
    
    volume_t = volume(err_t, 0)
    return volume_t[1]

def pv_micron(array, digits):
    #input is [mm], output is [um]
    array = array*1000 # mm -> um
    peak = np.nanmax(array)
    valley = np.nanmin(array)
   
    pv = peak - valley
    pv_str =str(round(pv, digits))
    
    return pv_str, pv
    
def rms_micron(array, digits):
    #input is [mm], output is [um]
    data = array*1000
    sigma = np.nansum(data**2)
    num = np.count_nonzero(np.logical_not(np.isnan(array)))
    rms = np.sqrt(sigma / num)
    rms_str = str(round(rms, digits))
    return rms_str, rms

def volume(array, digits):
    # input is [mm] , output is [mm^3] and [mm]
    a_drop = array[~np.isnan(array)]
    threshold_idx = int(a_drop.size * ignore_offset)
    a_sort = np.sort(a_drop)
    threshold = a_sort[threshold_idx]
    
    offset = threshold
    
    offset_str = str(round(offset*1000, digits))
    
    array = array - offset
    unit_area = m1_radi**2/(px**2)
    unit_volume = array*unit_area
    volume = np.nansum(unit_volume)
    volume_str = str(round(volume, digits))
    return volume_str, volume, offset_str, offset

def assess(array, digits):
    # for plot axis title
    pv_t = pv_micron(array, digits)
    rms_t = rms_micron(array, digits)
    v_t = volume(array, digits)
    return pv_t[0], rms_t[0], v_t[0], v_t[2]

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
    # o 回転角[mrad], p 曲率半径[mm], q 軸外し距離[mm] 
    o0, do = 0, 0
    p0, dp = 8666, 0
    q0, dq = 1800, 0
    
    init = [o0 + do, p0 + dp, q0 + dq]
    
    ignore_offset = 0.02
    m1_radi = 1850/2
    
    zer_order = 12

    # 単位 m で入力
    zer_rms_arr = np.array([1.0 , 1.0010387658048345 , 1.0010387658048348 , 1.002072113267335 , 1.0020165068395825 , 1.0021407244794258 , 1.0031021910183147 , 1.0031021910183142 , 1.003119533322238 , 1.003119533322238 , 1.004118182121102 , 1.0042309263157316])
    zer_rms_arr = zer_rms_arr * 10 ** -6
    #-------------------------------------------------------------------------    
    
    px = 1023
    
    xx, yy = np.meshgrid(np.linspace(-m1_radi, m1_radi, px),np.linspace(-m1_radi, m1_radi, px))

    tf = np.where(xx**2+yy**2<m1_radi**2, True, False)
    
    for i in range(2, zer_order+1):
        
        zer_num = i
        zer_rms = zer_rms_arr[i-1]
        raw_for_minimize = zer_raw(zer_order, zer_num, zer_rms) * 1000 # 単位 mm に変更
        
        # minimize----------------------------------------------------------------
        start_time = time.time()  
        
        result = minimize(fun=volume_func, x0=init, method="Powell")
        
        end_time = time.time()
        
        print(end_time - start_time)
        
        # for plot ------------------------------------------------------------
        o_m, p_m, q_m = result["x"]
        
        para_0 = np.where(tf==True, cz1(o0, p0, q0, xx, yy), np.nan)
        para_m = np.where(tf==True, cz1(o_m, p_m, q_m, xx, yy), np.nan)
        
        diff_0 = para_0 - para_0
        diff_m = para_m - para_0
        
        surf = para_0 + raw_for_minimize
        
        err_0 = surf - para_0
        err_m = surf - para_m
        
        assess_0 = assess(err_0, 2)
        assess_m = assess(err_m, 2)
    
        # figure plot ------------------------------------------------------------
        fig = plt.figure(figsize=(10, 10))
        
        title_para_0 = "Standard Parabola\n"\
            + "p=" + str(p0) + " [mm]\n"\
                + "q=" + str(q0) + " [mm]\n"\
                    + r"$\phi$=" + str(o0) + " [mrad]"
        title_para_m = "Minimized Parabola\n"\
            + "p=" + str(p_m) + " [mm]\n"\
                + "q=" + str(q_m) + " [mm]\n"\
                    + r"$\phi$=" + str(o_m) + " [mrad]"
        
        title_err_0 = "Zernike order " + str(zer_num) + "\n"\
            + "P-V = " + assess_0[0] + r" [$\mu$m]  RMS = " + assess_0[1] + r" [$\mu$m]" + "\n"\
                + "Removed = " + assess_0[2] + r" [mm$^3$]" + "\n"\
                    + " with offset = " + assess_0[3] + r" [$\mu$m]\n"\
                        + "( ignore_offset : " + str(ignore_offset*100) + "%)"
        
        title_err_m = "Zernike order " + str(zer_num) + "\n"\
            + "P-V = " + assess_m[0] + r" [$\mu$m]  RMS = " + assess_m[1] + r" [$\mu$m]" + "\n"\
                + "Removed = " + assess_m[2] + r" [mm$^3$]" + "\n"\
                    + " with offset = " + assess_m[3] + r" [$\mu$m]\n"\
                        + "( ignore_offset : " + str(ignore_offset*100) + "%)\n"
        
        ax_para_0 = image_plot(fig, title_para_0, 221, diff_0, diff_m)
        ax_para_m = image_plot(fig, title_para_m, 222, diff_m, diff_m)
        ax_err_0 = image_plot(fig, title_err_0, 223, err_0, err_0)
        ax_err_m = image_plot(fig, title_err_m, 224, err_m, err_m)
        
        fig.tight_layout()
        
        picname = mkfolder() + "/zer" + str(zer_num) + "init_do" + str(do) + "dp" + str(dp) + "dq" + str(dq) + ".png"
        fig.savefig(picname)
        fig.clf()