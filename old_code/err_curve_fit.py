# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 19:27:09 2021

@author: swimc
"""

import numpy as np
import pandas as pd
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

def cz1(p, cx, cy, aqx, phi):    
    a = 1/(2*p)
    aqz = a * ( aqx**2 + 0**2 )
    theta = np.arctan(2*a*aqx)
    
    D = -4*a**2*cx**2*np.sin(phi)**2*np.sin(theta)**2 - 8*a**2*cx*cy*np.sin(phi)*np.sin(theta)**2*np.cos(phi) + 4*a**2*cy**2*np.sin(phi)**2*np.sin(theta)**2 - 4*a**2*cy**2*np.sin(theta)**2 + 2*a*aqx*np.sin(2*theta) + 4*a*aqz*np.sin(theta)**2 + 4*a*cx*np.sin(theta)*np.cos(phi) - 4*a*cy*np.sin(phi)*np.sin(theta) - np.sin(theta)**2 + 1

    cz1 = (4*a*aqx*np.sin(theta) - a*cx*(np.sin(phi - 2*theta) - np.sin(phi + 2*theta)) - a*cy*(np.cos(phi - 2*theta) - np.cos(phi + 2*theta)) - 2*np.sqrt(D) + 2*np.cos(theta))/(4*a*np.sin(theta)**2)
    return cz1

def func(independent, o, p, q):
    x = independent[0]
    y = independent[1]
    return cz1(p, x, y, q, o/1000)

def image_plot(fig, title, position, c, c_scale, cb_micron=True):
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

def pv(array, digits):
    peak = np.nanmax(array)
    valley = np.nanmin(array)
   
    pv = peak - valley
    pv_str =str(round(pv, digits))
    
    return pv_str, pv
    
def rms(array, offset, digits):
    data = array + offset
    sigma = np.nansum(data**2)
    num = np.count_nonzero(np.logical_not(np.isnan(array)))
    rms = np.sqrt(sigma / num)
    rms_str = str(round(rms, digits))
    return rms_str, rms

def volume(array, digits):
    array = array - np.nanmin(array)
    unit_area = m1_radi**2/(px**2)
    unit_volume = array*unit_area
    volume = np.nansum(unit_volume)
    volume_str = str(round(volume, digits))
    return volume_str, volume



if __name__ == "__main__":
 
    
    # parametar --------------------------------------------------------------
    # p 曲率半径, q 軸外し距離, o 回転角phi[mrad]
    
    o0, do = 0, 0
    p0, dp = 8666, 0
    q0, dq = 1800, 1
    
    init = [o0 + do, p0 + dp, q0 +dq]
    #-------------------------------------------------------------------------
    
    px = 1023
    m1_radi = 1850/2
    
    raw = dat_read("digitFig01.csv") # dat file 読み出し
    #raw = raw[:,::-1].T #   反時計回りに 0.5*pi 回転
    #raw = raw[::-1, ::-1] # 反時計回りに 1.0*pi 回転
    #raw = raw[::-1, :].T #  反時計回りに 1.5*pi 回転
    
    raw_1d = raw.flatten()
    
    tf = ~np.isnan(raw) # np.nanが False
    tf_1d = tf.flatten()

    x_arr = np.linspace(-m1_radi, m1_radi, px)
    y_arr = np.linspace(-m1_radi, m1_radi, px)
    xx, yy = np.meshgrid(x_arr, y_arr)
    
    ind = np.array([xx.flatten()[tf_1d], yy.flatten()[tf_1d]]) # independent variable 独立変数
    para_0_drop = cz1(p0, ind[0], ind[1], q0, o0)
    
    dep = para_0_drop + raw_1d[tf_1d] # dependent variable 従属変数
    
    param, covariance = sp.optimize.curve_fit(func, ind, dep, p0=init)
    
    o_m, p_m, q_m = param
    
    # for plot----------------------------------------------------------------
    
    para_0 = np.where(tf==True, cz1(p0, xx, yy, q0, o0/1000), np.nan)
    para_m = np.where(tf==True, cz1(p_m, xx, yy, q_m, o_m/1000), np.nan)
    
    diff_0 = para_0 - para_0
    diff_m = para_m - para_0
    
    surf = para_0 + raw
    
    err_0 = surf - para_0
    err_m = surf - para_m
    
    pv_0 = pv(err_0*1000, 2)
    pv_m = pv(err_m*1000, 2)
    
    rms_0 = rms(err_0*1000, 0, 2)
    rms_m = rms(err_m*1000, 0, 2)
    
    volume_0 = volume(err_0, 0)
    volume_m = volume(err_m, 0)
    
    # figure plot ------------------------------------------------------------
    fig = plt.figure(figsize=(10, 10))
    cmap = cm.jet
    cmap.set_bad("white") # np.nanに対する色を指定
    fs = 15
    
    title_para_0 = "Standard Parabola" + "\n"\
        + "RoC : p = " + str(p0) + " [mm]\n"\
            + "  Off-Axis : q = " + str(q0) + " [mm]\n"\
                + r"Clocking : $\phi$ = " + str(o0) + " [mrad]"

    title_para_m = "Minimized Parabola" + "\n"\
        + "RoC : p = " + str(p_m) + " [mm]\n"\
            + "  Off-Axis : q = " + str(q_m) + " [mm]\n"\
                + r"Clocking : $\phi$ = " + str(o_m) + " [mrad]"
    
    title_err_0 = "P-V = " + pv_0[0] + r" [$\mu$m]  RMS = " + rms_0[0] + r" [$\mu$m]" + "\n"\
        + "Removed = " + volume_0[0] + r" [mm$^3$]" + "\n"\
            + " with offset = " + str(round(np.nanmin(err_0*1000), 2)) + r" [$\mu$m]"
            
    title_err_m = "P-V = " + pv_m[0] + r" [$\mu$m]  RMS = " + rms_m[0] + r" [$\mu$m]" + "\n"\
        + "Removed = " + volume_m[0] + r" [mm$^3$]" + "\n"\
            + " with offset = " + str(round(np.nanmin(err_m*1000), 2)) + r" [$\mu$m]"
    
    ax_para_0 = image_plot(fig, title_para_0, 221, diff_0, diff_m)
    ax_err_0 = image_plot(fig, title_err_0, 222, err_0, err_0)
    ax_para_m = image_plot(fig, title_para_m, 223, diff_m, diff_m)
    ax_err_m = image_plot(fig, title_err_m, 224, err_m, err_m)
    fig.tight_layout()
    
    picname = mkfolder() + "/curve_fit_do" + str(do) + "dp" + str(dp) + "dq" + str(dq) + ".png"
    fig.savefig(picname)