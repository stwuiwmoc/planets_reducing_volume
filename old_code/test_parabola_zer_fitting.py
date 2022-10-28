# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 17:36:23 2021

@author: swimc

o, p, q の単位変化に対する鏡面形状の変化をplotし、その変化量をZernikeに近似して各項のrmsをplot
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
    filename = os.path.basename(__file__)
    filename = filename.replace(".py", "") + suffix
    folder = "mkfolder/" + filename + "/"
    os.makedirs(folder, exist_ok=True)
    return folder

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
    offset = np.nanmin(array)
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
    p0, dp = 8667, 0
    q0, dq = 1800, 0
    
    init = [o0 + do, p0 + dp, q0 + dq]
    
    m1_radi= 1850/2
    
    zer_order = 11
    #-------------------------------------------------------------------------    
    
    o_t, p_t, q_t = init
    
    px = 1023
    
    xx, yy = np.meshgrid(np.linspace(-m1_radi, m1_radi, px),np.linspace(-m1_radi, m1_radi, px))

    raw = dat_read("digitFig01.csv") # dat file 読み出し
    #raw = raw[:,::-1].T #   反時計回りに 0.5*pi 回転
    #raw = raw[::-1, ::-1] # 反時計回りに 1.0*pi 回転
    #raw = raw[::-1, :].T #  反時計回りに 1.5*pi 回転

    tf = ~np.isnan(raw)

    para_0 = np.where(tf==True, cz1(o0, p0, q0, xx, yy), np.nan)
    
    para_t = np.where(tf==True, cz1(o_t, p_t, q_t, xx, yy), np.nan)
    
    para_diff = para_t - para_0
    
    assess_t = assess(para_diff, 2)
    
    zer, fit = pr.prop_fit_zernikes(para_diff/1000, tf, px/2, zer_order, xc=px/2, yc=px/2, FIT=True)
    
    # figure plot-------------------------------------------------------------
    fig = plt.figure(figsize=(5,7))
    
    title_diff = "Test Parabola\n"\
        + "p=" + str(p_t) + " [mm] "\
            + "q=" + str(q_t) + " [mm] "\
                + r"$\phi$=" + str(o_t) + " [mrad]"
    
    ax_diff = image_plot(fig, title_diff, 211, para_diff, para_diff)
    ax_x = np.arange(1, zer_order+1)
    ax_zer = fig.add_subplot(212)
    ax_zer.plot(ax_x, zer*10**6, marker="s", linewidth=1)
    ax_zer.grid()
    ax_zer.set_xlabel("zernike order")
    ax_zer.set_xticks(np.arange(1, zer_order+1))
    ax_zer.set_ylabel(r"zernike RMS [$\mu$m]")
    
    fig.tight_layout()
    picname = mkfolder() + "/do" + str(do) + "dp" + str(dp) + "dq" + str(dq) + ".png"
    fig.savefig(picname)
    
    
    