# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 22:54:40 2021

@author: swimc
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
import scipy as sp
import os
import time
import proper as pr
import PIL

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
    pixels = pixels / 10**3 # micron -> mm
    
    return pixels

def zer_2_image(zer_order, zer_arr):
    #imput zer_arr [mm] -> output 2d-array [mm]
    zer_arr = zer_arr/1000
    wavelength = 500*10**-9
    diam = 2 * m1_radi / 1000 # m 単位の直径に直す
    zer_num = np.arange(1, zer_order+1)
    
    wavestruct = pr.prop_begin(diam, wavelength, px, 1)
    wfe = pr.prop_zernikes(wavestruct, zer_num, zer_arr)
    wfe = wfe * 10**3
    return wfe

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
    m1_radi = 1850/2
    px = 1023
    xx, yy = np.meshgrid(np.linspace(-m1_radi, m1_radi, px),np.linspace(-m1_radi, m1_radi, px))
    zer_order = 10
    
    raw = dat_read("digitFig01.csv")
    tf = ~np.isnan(raw)
    mask = np.where(tf==True, 1, np.nan)
    
    err_m = raw
    
    err_zer = pr.prop_fit_zernikes(err_m/1000, tf, px/2, zer_order, xc=px/2, yc=px/2)
    err_zer_drop = np.delete(err_zer, [0], 0)
    om = np.genfromtxt("WT03_zer10_opration_matrix[m].csv", delimiter=",").T
    om_drop = np.delete(om, [0], 0)
    
    om_inv = sp.linalg.pinv2(om_drop * 10**9) * 10**9
    
    force = np.dot(om_inv, err_zer_drop)
    
    reprod_zer_drop = np.dot(om_drop, force)
    reprod_zer = np.insert(reprod_zer_drop, 0, [0])
    reprod = mask * zer_2_image(zer_order, reprod_zer * 10**3)
    
    residual = raw - reprod
    
    fig = plt.figure(figsize=(5,10))
    ax_raw = image_plot(fig, "raw", 211, raw, raw)
    ax_res = image_plot(fig, "residual", 212, residual, raw)
    fig.tight_layout()