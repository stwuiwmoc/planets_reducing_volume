# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 14:13:46 2021

@author: swimc

powell法で体積を最小化する
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
import scipy as sp
import os
import time


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

def dat_read(fname):
    # 鏡面異常を表す dat file を読み出して成形
    data = np.genfromtxt("digitFig01.csv", delimiter = ",", encoding="utf-8_sig")
    
    pixels = data.reshape((1024,1026))
    
    #縦横の5.52182だけが並ぶ行と列を削除して1024*1023に成形
    pixels = np.delete(pixels, 0, 0)
    pixels = np.delete(pixels, [0, 1024, 1025], 1)
    
    #5.52182は円の外なので nanに置き換え
    pixels = np.where(pixels==5.52182, np.nan, pixels)
    pixels = pixels / 1000 # micron -> mm
    
    return pixels

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
    unit_area = m1_radi**2/(array.size)
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
    ignore_offset = 0.03
    o0, do = 0, 0
    p0, dp = 8666, 0
    q0, dq = 1800, 0
    
    init = [o0 + do, p0 + dp, q0 + dq]
    
    #-------------------------------------------------------------------------    
    m1_radi = 1850/2
    px = 1023
    
    xx, yy = np.meshgrid(np.linspace(-m1_radi, m1_radi, px),np.linspace(-m1_radi, m1_radi, px))

    raw = dat_read("digitFig01.csv") # dat file 読み出し
    #raw = raw[:,::-1].T #   反時計回りに 0.5*pi 回転
    #raw = raw[::-1, ::-1] # 反時計回りに 1.0*pi 回転
    #raw = raw[::-1, :].T #  反時計回りに 1.5*pi 回転

    tf = ~np.isnan(raw)
    
    # raw data edit-----------------------------------------------------------
        
    
    mask = np.where(tf==True, 1, np.nan)
    raw_0f = np.where(tf==True, raw, 0)

    filter_num = int(input("Non filter : 1, gaussian : 2, median : 3\nInput filter_num = "))
    
    if filter_num == 1:
        raw_for_minimize = raw
        g_sigma = 0
        m_size = 0
    
    if filter_num == 2:
        g_sigma = int(input("Input gaussian sigma = "))
        gaussian_0f = mask * sp.ndimage.filters.gaussian_filter(raw_0f, g_sigma)
        raw_for_minimize = gaussian_0f    
        m_size = 0
    
    if filter_num == 3:
        m_size = int(input("Input filter size = "))
        median_0f = mask * sp.ndimage.filters.median_filter(raw_0f, m_size)
        raw_for_minimize = median_0f
        g_sigma = 0
    
    else:
        pass
    
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
    
    filter_str = ("Non filter\n",\
                 r"Gaussian filter : $\sigma$ =" + str(g_sigma) + "\n",\
                      "Median filter : size = " + str(m_size) + r"$\times$" + str(m_size) + "\n")
    
    title_err_0 =filter_str[filter_num-1]\
        + "P-V = " + assess_0[0] + r" [$\mu$m]  RMS = " + assess_0[1] + r" [$\mu$m]" + "\n"\
            + "Removed = " + assess_0[2] + r" [mm$^3$]" + "\n"\
                + " with offset = " + assess_0[3] + r" [$\mu$m] ( Ignore Lower " + str(ignore_offset*100) + "% )" 
    
    title_err_m = filter_str[filter_num-1]\
        + "P-V = " + assess_m[0] + r" [$\mu$m]  RMS = " + assess_m[1] + r" [$\mu$m]" + "\n"\
            + "Removed = " + assess_m[2] + r" [mm$^3$]" + "\n"\
                + " with offset = " + assess_m[3] + r" [$\mu$m] ( Ignore Lower " + str(ignore_offset*100) + "% )" 
    
    ax_para_0 = image_plot(fig, title_para_0, 221, diff_0, diff_m)
    ax_para_m = image_plot(fig, title_para_m, 222, diff_m, diff_m)
    ax_err_0 = image_plot(fig, title_err_0, 223, err_0, err_0)
    ax_err_m = image_plot(fig, title_err_m, 224, err_m, err_m)
    
    fig.tight_layout()
    
    picname1 = mkfolder() + "/" + "init_do" + str(do) + "dp" + str(dp) + "dq" + str(dq) + ".png"
    picname2 = mkfolder() + "/g" + str(g_sigma) + "init_do" + str(do) + "dp" + str(dp) + "dq" + str(dq) + ".png"
    picname3 = mkfolder() + "/m" + str(m_size) + "init_do" + str(do) + "dp" + str(dp) + "dq" + str(dq) + ".png"
    
    picname = (picname1, picname2, picname3)
    
    fig.savefig(picname[filter_num-1])