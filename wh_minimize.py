# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 20:11:08 2021

@author: swimc

逆行列を計算し、wh駆動により最小化
"""

import numpy as np
import os
import PIL
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
    a_drop = array[~np.isnan(array)]
    threshold_idx = int(a_drop.size * ignore_offset)
    a_sort = np.sort(a_drop)
    threshold = a_sort[threshold_idx]
    
    offset = threshold
    offset_str = str(round(offset*1000, digits))
    
    array = abs(array - offset)
    unit_area = m1_radi**2/(array.size)
    unit_volume = array*unit_area
    volume = np.nansum(unit_volume)
    volume_str = str(round(volume, digits))
    return volume_str, volume, offset_str, offset

def assess(array, digits):
    pv_t = pv_micron(array, digits)
    rms_t = rms_micron(array, digits)
    v_t = volume(array, digits)
    return pv_t[0], rms_t[0], v_t[0], v_t[2]

def image_resize(array, new_px):
    img = PIL.Image.fromarray(array)
    resized = img.resize((new_px,new_px))
    resized = np.asarray(resized)
    return resized

def zer_raw(zer_order, zer_num):
    # input int, output 2d_array [m]
    # 2um rmsを作るinput
    zer_rms_arr = (10**-6) *  np.array([2.0,2.0020775316096695,2.002077531609669,2.0041442265346694,2.0040330136791655,2.0042814489588507,2.006204382036629,2.006204382036629,2.0062390666444756,2.006239066644476,2.008236364242204,2.0084618526314637,2.0080541567495085,2.0085794355167943,2.008066811409963,2.010253171410069,2.010253171410068,2.0103049772121215,2.0103049772121215,2.0104092451902407,2.0104092451902407])
    
    wavelength = 500*10**-9
    diam = m1_radi / 1000 * 2 # m 単位の直径に直す

    wavestruct = pr.prop_begin(diam, wavelength, px, 1)
    wfe = pr.prop_zernikes(wavestruct, [zer_num], [zer_rms_arr[zer_num-1]])
    return wfe


def zer_2_image(zer_order, zer_arr):
    #imput zer_arr [m] -> output 2d-array [mm]
    wavelength = 500*10**-9
    diam = 2 * m1_radi / 1000 # m 単位の直径に直す
    zer_num = np.arange(1, zer_order+1)
    
    wavestruct = pr.prop_begin(diam, wavelength, px, 1)
    wfe = pr.prop_zernikes(wavestruct, zer_num, zer_arr)
    wfe = wfe * 10**3
    return wfe

def torque_plot(fig, title, position, torque):
    fs = 15
    x = np.arange(1, 12)

    ax = fig.add_subplot(position)    
    ax.plot(x, torque[0:11], color="black", marker="s", linewidth=1)
    ax.plot(x, torque[11:22], color="red", marker="s", linewidth=1)
    ax.plot(x, torque[22:33], color="blue", marker="s", linewidth=1)
    
    ax.set_title(title, fontsize=fs)
    ax.grid()
    ax.set_xlabel("Warping Harness Number", fontsize=fs)
    ax.set_xticks(x)
    ax.set_ylabel("Support Force", fontsize=fs)
    
    return ax

if __name__ == '__main__':
    raw_px = 1023
    px = 1023
    m1_radi = 1850/2
    zer_order = 21
    xx, yy = np.meshgrid(np.linspace(-m1_radi, m1_radi, px),np.linspace(-m1_radi, m1_radi, px))
    ignore_offset = 0.02
    
    raw = dat_read("digitFig01.csv")
    
    zer_num = 7
    
    tf = np.where(xx**2+yy**2<m1_radi**2, True, False)
    mask = np.where(tf==True, 1, np.nan)
    mesured = mask * zer_raw(zer_order, zer_num) * 1000
    """
    #mesured = image_resize(raw, px)
    mesured = raw
    """
    tf = ~np.isnan(mesured)
    mask = np.where(tf==True, 1, np.nan)
    mesured_0f = np.where(tf==True, mesured, 0)
    mesured_z, fit = pr.prop_fit_zernikes(mesured/1000, tf, px/2, zer_order, xc=px/2, yc=px/2, FIT=True)
    fit = mask * fit * 1000
    
    om = np.genfromtxt("zer_opration_matrix[m].csv", delimiter=",").T
    om_0f = np.where(np.isnan(om)==True, 0, om)
    
    om_inv = np.linalg.pinv(om_0f)
    
    torque = np.dot(om_inv, mesured_z)
    
    reprod_z = np.dot(om_0f, torque)
    reprod = mask * zer_2_image(zer_order, reprod_z)
    
    diff_re = mesured-reprod
    diff_f = fit-reprod
    
    #figure plot---------------------------------------------------------------
    fig = plt.figure(figsize=(10,15))
    
    assess_raw = assess(mesured, 2)
    assess_f = assess(fit, 2)
    assess_d_re = assess(diff_re, 2)
    assess_d_f = assess(diff_f, 4)
    
    title_raw = "raw ( ignore_offset : " + str(ignore_offset*100) + "%)\n"\
        + "P-V = " + assess_raw[0] + r" [$\mu$m]  RMS = " + assess_raw[1] + r" [$\mu$m]" + "\n"\
            + "Removed = " + assess_raw[2] + r" [mm$^3$]" + "\n"\
                + " with offset = " + assess_raw[3] + r" [$\mu$m]"
                
    title_f = "zernike fit\n" + assess_f[0] + " / " + assess_f[1] + " / " + assess_f[2] + " / " + assess_f[3]
    
    title_d_re = "raw - Reproduction\n" + assess_d_re[0] + " / " + assess_d_re[1] + " / " + assess_d_re[2] + " / " + assess_d_re[3]
    title_d_f = "fit - Reproduction\n" + assess_d_f[0] + " / " + assess_d_f[1] + " / " + assess_d_f[2] + " / " + assess_d_f[3]
    
    ax_raw = image_plot(fig, title_raw, 321, mesured, mesured)
    ax_fit = image_plot(fig, title_f, 322, fit, mesured)
    
    ax_torque = torque_plot(fig, "", 323, torque)
    
    ax_re = image_plot(fig, "Reproduction", 324, reprod, mesured)
    
    ax_raw_diff = image_plot(fig, title_d_re, 325, diff_re, diff_re)
    ax_fit_diff = image_plot(fig, title_d_f, 326, diff_f, diff_f)
   
    fig.tight_layout()
    