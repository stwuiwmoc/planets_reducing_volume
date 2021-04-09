# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 13:54:26 2021

@author: swimc

mesured parameter の谷を減らすための平滑化、どのフィルターでどんなぼかしになるかの検証
"""

import numpy as np
import pandas as pd
import os
import time
import proper as pr
import scipy as sp
import pykrige as krg
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
    pixels = pixels / 1000 # micron -> mm
    
    return pixels

def image_plot(fig, title, position, c, c_scale, cb_micron=True):
    cmap = cm.get_cmap("jet")
    fs = 15
    
    ax = fig.add_subplot(position)
    ax.imshow(c, interpolation="nearest", cmap=cmap, vmin=np.nanmin(c_scale), vmax=np.nanmax(c_scale),\
              origin="lower", extent=[-925, 925, -925, 925])
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
    pv_t = pv_micron(array, digits)
    rms_t = rms_micron(array, digits)
    v_t = volume(array, digits)
    return pv_t[0], rms_t[0], v_t[0], v_t[2]

def kriging_newpx(array, new_px):
    array_1d = array.flatten()
    x = xx.flatten()
    y = yy.flatten()
    
    x_grid = np.linspace(-m1_radi, m1_radi, new_px)
    y_grid = np.linspace(-m1_radi, m1_radi, new_px)
    
    ok_module = krg.ok.OrdinaryKriging(x, y, array_1d)
    z, sigma = ok_module.execute("grid", x_grid, y_grid)
    z = z.reshape((new_px,new_px))
    return z    

def image_resize(array, new_px):
    img = PIL.Image.fromarray(array)
    resized = img.resize((new_px,new_px))
    resized = np.asarray(resized)
    return resized

if __name__ == '__main__':
    
    px = 1023
    m1_radi, dr = 1850/2, 0
    r = m1_radi - dr
    
    
    x_arr = np.linspace(-m1_radi, m1_radi, px)
    y_arr = np.linspace(-m1_radi, m1_radi, px)
    xx, yy = np.meshgrid(x_arr, y_arr)
  
    raw = dat_read("digitFig01.csv")
    #tf = np.where(xx**2+yy**2<m1_radi**2, True, False)
    tf = ~np.isnan(raw)
    mask = np.where(tf==True, 1, np.nan)
    raw_0f = np.where(tf==True, raw, 0)
    
    g_sigma = 10
    gaussian = mask * sp.ndimage.filters.gaussian_filter(raw_0f, g_sigma)
    
    m_size = 10
    mean = mask * sp.ndimage.filters.uniform_filter(raw_0f, m_size)
    
    new_px = 256
    downscale = image_resize(raw_0f, new_px)
    upscale = image_resize(downscale, px)
    
    error = upscale - raw
    
    
    #kriging = kriging_newpx(resize, 100)
    
    # plot -------------------------------------------------------------------
    fig = plt.figure(figsize=(15,10))

    assess_raw = assess(raw, 2)
    assess_g = assess(gaussian, 2)
    assess_m = assess(mean, 2)
    assess_d = assess(downscale, 2)
    assess_u = assess(upscale, 2)
    assess_e = assess(error, 2)
    
    title_raw = "raw : radius = " + str(m1_radi) + "\n"\
        + "P-V = " + assess_raw[0] + r" [$\mu$m]  RMS = " + assess_raw[1] + r" [$\mu$m]" + "\n"\
            + "Removed = " + assess_raw[2] + r" [mm$^3$]" + "\n"\
                + " with offset = " + assess_raw[3] + r" [$\mu$m]"

    title_g = "gaussian " + str(g_sigma) + "\n" + assess_g[0] + " / " + assess_g[1] + "/ " + assess_g[2] + "/ " + assess_g[3]
    title_m = "mean " + str(m_size) + "\n" + assess_m[0] + " / " + assess_m[1] + " / " + assess_m[2] + " / " + assess_m[3]
    title_d = "resize " + str(new_px) + " px\n" + assess_d[0] + " / " + assess_d[1] + " / " + assess_d[2] + " / " + assess_d[3]
    title_u = "resize " + str(px) + " px\n" + assess_u[0] + " / " + assess_u[1] + " / " + assess_u[2] + " / " + assess_u[3]
    title_e = "error " + str(px) + " px\n" + assess_e[0] + " / " + assess_e[1] + " / " + assess_e[2] + " / " + assess_e[3]
    
    ax_raw = image_plot(fig, title_raw, 231, raw, raw)
    ax_g = image_plot(fig, title_g, 232, gaussian, raw_0f)
    ax_m = image_plot(fig, title_m, 233, mean, raw_0f)
    ax_d = image_plot(fig, title_d, 234, downscale, downscale)
    ax_u = image_plot(fig, title_u, 235, upscale, raw)
    ax_e = image_plot(fig, title_e, 236, error, error)
    
    fig.tight_layout()
    picname = mkfolder() + "/g" + str(g_sigma) + "m" + str(m_size) + ".png"
    fig.savefig(picname)