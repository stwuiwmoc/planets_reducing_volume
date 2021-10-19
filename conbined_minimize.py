# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 21:54:36 2021

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

def cz1(o, p, aqx, cx, cy):
    """
    順番が o, p, q に代わっているので注意

    Parameters
    ----------
    o : float [mrad]
        回転角
    p : float [mm]
        曲率半径
    aqx : floot [mm]
        軸外し距離 (off-axis)
    cx : 2D-mesh-array [mm]
    cy : 2D-mesh-array [mm]

    Returns
    -------
    cz1 : 2D-mesh-array [mm]
        input された o, p, q での切り取り放物面
    """
    
    phi = o / 1000 # [mrad] -> [rad] への換算
    
    a = 1/(2*p)
    aqz = a * ( aqx**2 + 0**2 )
    theta = np.arctan(2*a*aqx)
    
    D = -4*a**2*cx**2*np.sin(phi)**2*np.sin(theta)**2 - 8*a**2*cx*cy*np.sin(phi)*np.sin(theta)**2*np.cos(phi) + 4*a**2*cy**2*np.sin(phi)**2*np.sin(theta)**2 - 4*a**2*cy**2*np.sin(theta)**2 + 2*a*aqx*np.sin(2*theta) + 4*a*aqz*np.sin(theta)**2 + 4*a*cx*np.sin(theta)*np.cos(phi) - 4*a*cy*np.sin(phi)*np.sin(theta) - np.sin(theta)**2 + 1

    cz1 = (4*a*aqx*np.sin(theta) - a*cx*(np.sin(phi - 2*theta) - np.sin(phi + 2*theta)) - a*cy*(np.cos(phi - 2*theta) - np.cos(phi + 2*theta)) - 2*np.sqrt(D) + 2*np.cos(theta))/(4*a*np.sin(theta)**2)
    return cz1

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

def volume_func(X):
    o, p, q = X
    para_t = np.where(tf==True, cz1(o, p, q, xx, yy), np.nan)
    
    para_0 = np.where(tf==True, cz1(o0, p0, q0, xx, yy), np.nan)
    surf = para_0 + filtered
    
    err_t = surf - para_t
    
    volume_t = volume(err_t, 0)
    return volume_t[1]

def image_resize(array, new_px):
    img = PIL.Image.fromarray(array)
    resized = img.resize((new_px,new_px))
    resized = np.asarray(resized)
    return resized

def zer_raw(zer_order, zer_num):
    # input int, output 2d_array [mm]
    zer_rms_arr = 0.5 * (10**-6) * np.array([2.0,2.0020775316096695,2.002077531609669,2.0041442265346694,2.0040330136791655,2.0042814489588507,2.006204382036629,2.006204382036629,2.0062390666444756,2.006239066644476,2.008236364242204,2.0084618526314637,2.0080541567495085,2.0085794355167943,2.008066811409963,2.010253171410069,2.010253171410068,2.0103049772121215,2.0103049772121215,2.0104092451902407,2.0104092451902407])
    
    wavelength = 500*10**-9
    diam = m1_radi / 1000 * 2 # m 単位の直径に直す

    wavestruct = pr.prop_begin(diam, wavelength, px, 1)
    wfe = pr.prop_zernikes(wavestruct, [zer_num], [zer_rms_arr[zer_num-1]])
    wfe = wfe * 10**3
    return wfe

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

def nan_drop(masked_array_2d, mask_2d):
    # input 2d-array, mask -> output np.nanを除いた1d-array, 復元用のloc-1d
    array = masked_array_2d
    
    array_1d = array.flatten()
    mask_1d = mask_2d.flatten()
    loc = np.where(np.isnan(mask_1d)==True, np.nan, np.arange(mask_1d.size))
    loc_1d = loc.flatten()
    array_1d_drop = array_1d[~np.isnan(loc_1d)]
    loc_drop = loc_1d[~np.isnan(loc_1d)]
    return array_1d_drop, loc_drop

def drop_recovery(array_drop, loc_drop, pixel):
    # input nanを除いた1d-array, 復元用のloc-drop -> output 復元された 2d-array
    recovery_1d = np.empty(pixel**2) * np.nan
    for i in range(loc_drop.size):
        idx = int(loc_drop[i])
        recovery_1d[idx] = array_drop[i]
    recovery = recovery_1d.reshape((pixel,pixel))
    return recovery

# function for plot ----------------------------------------------------------
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

def torque_plot(fig, title, position, force):
    # input force [mm] 
    fs = 15
    x = np.arange(1, 13)
    x_str = []
    for i in range(12):
        text = ""
        for j in range(3):
            num = str(12*j + i + 1).zfill(2)
            text = text + "\n" + num
        x_str.append(text)

    torque = force
    
    ax = fig.add_subplot(position)
    ax.plot(x, torque[0:12], color="black", marker="s", linewidth=1)
    ax.plot(x, torque[12:24], color="green", marker="o", linewidth=1)
    ax.plot(x, torque[24:36], color="darkviolet", marker="^", linewidth=1)
    
    ax.set_title(title, fontsize=fs)
    ax.grid()
    ax.set_xlabel("Motor Number", fontsize=fs)
    ax.set_xticks(x)
    ax.set_xticklabels(x_str)
    ax.set_ylabel("Motor drive amount [mm]", fontsize=fs)
    
    ax.hlines([5, -5], xmin=1, xmax=12, color ="red", linestyle="dashed")
    return ax

if __name__ == '__main__':
    ## parametar --------------------------------------------------------------
   
    # o 回転角[mrad], p 曲率半径[mm], q 軸外し距離[mm]     
    o0, do = 0, 0
    p0, dp = 8667, 0
    q0, dq = 1800, 0
    
    # フチを無視する長さ [mm] 半径に対して計算
    ignore_radi = 50
    
    # 研磨量計算時に下位 %を無視して offset を設定 入力は％ではなく小数
    ignore_offset = 0.02
    
    # constant ---------------------------------------------------------------
    m1_radi = 1850/2
    px = 1023
    xx, yy = np.meshgrid(np.linspace(-m1_radi, m1_radi, px),np.linspace(-m1_radi, m1_radi, px))
    zer_order = 10
    
    ## option -----------------------------------------------------------------
    option_loop = int(input("Input : 0, Loop : 1\nInput Code type -> "))

    if option_loop == 0:
        
        option_raw = int(input("Mesured data : 0, Zernike : 1 ~ 21\nInput data_type -> "))
        option_filter = int(input("filter select\nNon filter : 0, Mean : 1, Gaussian : 2, Median : 3\nInput filter_num -> "))
        filter_param = int(input("Input filter size or sigma  -> "))
        option_ideal = int(input("Ideal Surface Minimize\nSkip : 0, Minimize : 1\nInput -> "))
        option_wh = int(input("Warping Harness Minimize\nSkip : 0, Minimize(zer) : 4\nInput -> "))
        loop_range = [option_raw, option_raw + 1]
        
    if option_loop == 1:    
        option_filter = int(input("filter select\nNon filter : 0, Mean : 1, Gaussian : 2, Median : 3\nInput filter_num -> "))
        filter_param = int(input("Input filter size or sigma  -> "))
        option_ideal = int(input("Ideal Surface Minimize\nSkip : 0, Minimize : 1\nInput -> "))
        option_wh = int(input("Warping Harness Minimize\nSkip : 0, Minimize(zer) : 4\nInput -> "))
        
        loop_range = [2, zer_order + 1]
        
    else:
        pass
    
    for i in range(loop_range[0], loop_range[1]):
        print(i)
        option_raw = i


        ## 対象データ入力 -----------------------------------------------------------
        if option_raw == -1:    
            raw = dat_read("raw_data/diff.csv")
            raw_0f = raw

        if option_raw == 0:    
            #raw = dat_read("digitFig01.csv") # dat file 読み出し
            raw = np.loadtxt("mkfolder/stitch2mesh/diff.csv") * 1e-6
            print(raw.max())
            #raw = raw[:,::-1].T #   反時計回りに 0.5*pi 回転
            #raw = raw[::-1, ::-1] # 反時計回りに 1.0*pi 回転
            #raw = raw[::-1, :].T #  反時計回りに 1.5*pi 回転
            tf_raw = ~np.isnan(raw)
            raw_0f = np.where(tf_raw==True, raw, 0)
            """
            img = PIL.Image.fromarray(raw_0f)
            img_45 = img.rotate(-45)
            raw_0f = np.asarray(img_45)
            """

        else:
            
            raw = zer_raw(zer_order, option_raw)
            raw = np.where(xx**2+yy**2<m1_radi**2, raw, np.nan)
            raw_0f = np.where(~np.isnan(raw)==True, raw, 0)
            """
            tf = np.where(xx**2+yy**2<m1_radi**2, True, False)

            """
        
        ## フチの処理 ------------------------------------------------------------
        valid_radi = m1_radi - ignore_radi # 有効半径
        tf = np.where(xx**2+yy**2<valid_radi**2, True, False)
        mask = np.where(tf==True, 1, np.nan)
        
        
        ## フィルタ処理 -----------------------------------------------------------  
        raw_0f = np.where(tf==True, raw_0f, 0)
        
        if option_filter == 0:
            filter_param = 0
            filtered = mask * raw
            
        if option_filter == 1:
            filtered = mask * sp.ndimage.filters.uniform_filter(raw_0f, size=filter_param)
        
        if option_filter == 2:
            filtered = mask * sp.ndimage.filters.gaussian_filter(raw_0f, filter_param)
        
        if option_filter == 3:
            filtered = mask * sp.ndimage.filters.median_filter(raw_0f, filter_param)
            
        else:
            pass
        
        ## 理想鏡面での最小化--------------------------------------------------
        
        init = [o0 + do, p0 + dp, q0 + dq] # 理想鏡面最小化での初期値
        
        if option_ideal == 1:
            print("\nStart Ideal surface minimize")
            
            start_time = time.time()    
        
            result = minimize(fun=volume_func, x0=init, method="Powell")
            
            end_time = time.time()
            print(end_time - start_time)
                
            o_m, p_m, q_m = result["x"]
            
        if option_ideal == 0:
            o_m, p_m, q_m = o0, p0, q0
            
        else:
            pass
    
        para_0 = np.where(tf==True, cz1(o0, p0, q0, xx, yy), np.nan)
        para_m = np.where(tf==True, cz1(o_m, p_m, q_m, xx, yy), np.nan)
        
        para_diff = para_m - para_0
        surf = para_0 + filtered
        
        err_0 = surf - para_0
        err_m = surf - para_m
        
        offset = volume(err_m, 2)[3]
        
        #err_m = err_m - offset
        
        ## WH minimize-------------------------------------------------------------
        px_s = 256
        xx_256, yy_256 = np.meshgrid(np.linspace(-m1_radi, m1_radi, px_s),np.linspace(-m1_radi, m1_radi, px_s))
        tf_256 = np.where(xx_256**2+yy_256**2<=valid_radi**2, True, False)  
        mask_256 = np.where(tf_256==True, 1, np.nan)
        err_256 = mask_256 * image_resize(np.nan_to_num(err_m, nan=0), px_s)
                    
        if option_wh == 2:
            print("\nStart Warping Harness minimize")    
            
            err_drop, loc_drop = nan_drop(err_256, mask_256)
            
            om = np.genfromtxt("WT03_256_opration_matrix[m].csv", delimiter=",").T
            om_inv = sp.linalg.pinv(om * 10**9) * 10**9
            force = np.dot(om_inv, err_drop/1000)
            
            reprod_drop = np.dot(om, force)
            reprod_256 = drop_recovery(reprod_drop, loc_drop, px_s) * 10**3
            reprod = mask * image_resize(reprod_256, px)
                 
        if option_wh == 4:
            print("\nStart Warping Harness minimize")
            
            err_meaned = err_256 - np.nanmean(err_256)
            err_zer = pr.prop_fit_zernikes(err_meaned/1000, tf_256, px_s/2, zer_order, xc=px_s/2, yc=px_s/2)
            err_zer_drop = np.delete(err_zer, [0], 0) # piston, tilt成分を除去 
            
            om = np.genfromtxt("WT06_zer10_opration_matrix[m].csv", delimiter=",").T
            om_drop = np.delete(om, [0], 0) # piston, tilt成分を除去 
            
            om_inv = sp.linalg.pinv2(om_drop * 10**9) * 10**9
            
            force = np.dot(om_inv, err_zer_drop)
            
            reprod_zer_drop = np.dot(om_drop, force)
            #reprod_zer = np.insert(reprod_zer_drop, 0, err_zer[0:3])
            reprod_zer = np.insert(reprod_zer_drop, 0, [0])
            
            reprod = mask * zer_2_image(zer_order, reprod_zer * 10**3)
        
        
        if option_wh == 0:
            force = np.zeros(36)
            reprod = mask * np.zeros((px,px))
        
        else:
            pass
        
        residual = err_m - reprod
        ## figure plot ------------------------------------------------------------
        assess_raw = assess(raw, 2)
        assess_f = assess(filtered, 2)
        assess_i = assess(err_m, 2)
        assess_res = assess(residual, 4)
        
        title_raw_str = [ "Zernike " + str(i) for i in range(0,22)]
        title_raw_str[0] = "Raw"
        
        title_raw = title_raw_str[option_raw] + "\n"\
            + "Aperture = " + str(2*m1_radi) + " [mm]\n"\
                + "P-V = " + assess_raw[0] + r" [$\mu$m]  RMS = " + assess_raw[1] + r" [$\mu$m]" + "\n"\
                    + "Removed = " + assess_raw[2] + r" [mm$^3$]" + "\n"\
                        + " with offset = " + assess_raw[3] + r" [$\mu$m]" + "\n"\
                            + "( Ignore lower " + str(ignore_offset*100) + " %)"
        
        title_filter_str = ["Non filter", 
                      "Mean filter (size " + str(filter_param) + r"$\times$" + str(filter_param) + "px)",
                      r"Gaussian filter ($\sigma$ = " + str(filter_param),
                      "Median filter (size " + str(filter_param) + r"$\times$" + str(filter_param) + "px)"]
        
        title_f = title_filter_str[option_filter] +"\n"\
            + "Effective Aperture = " + str(2*valid_radi) + " [mm]\n"\
                + "P-V = " + assess_f[0] + r" [$\mu$m]  RMS = " + assess_f[1] + r" [$\mu$m]" + "\n"\
                    + "Removed = " + assess_f[2] + r" [mm$^3$]" + "\n"\
                        + " with offset = " + assess_f[3] + r" [$\mu$m]" + "\n"\
                            + "( Ignore lower " + str(ignore_offset*100) + " %)"
        
        title_d = "\n" + "Minimized Parabola\n"\
            + "p=" + str(round(p_m, 3)) + " [mm]\n"\
                + "q=" + str(round(q_m, 3)) + " [mm]\n"\
                    + r"$\phi$=" + str(round(o_m, 3)) + " [mrad]"
        
        title_i_str = ["Skipped (Same as above result)",
                       "Ideal Surface Minimize"]
        
        title_i = "\n" + title_i_str[option_ideal] + "\n"\
            + "P-V = " + assess_i[0] + r" [$\mu$m]  RMS = " + assess_i[1] + r" [$\mu$m]" + "\n"\
                + "Removed = " + assess_i[2] + r" [mm$^3$]" + "\n"\
                    + " with offset = " + assess_i[3] + r" [$\mu$m]" + "\n"\
                        + "( Ignore lower " + str(ignore_offset*100) + " %)"
    
        title_rep = "Active Support Reproduction"
        
        title_res_str = ["Skipped (Same as above result)",
                         "Warping Harness Minimize (zernike)",
                         "Warping Harness Minimize (256^2)",
                         "Warping Harness Minimize (zer_inv + 256^2)",
                         "Active Support Minimize"]
        
        
        title_res = "\n" + title_res_str[option_wh] + "\n"\
            + "P-V = " + assess_res[0] + r" [$\mu$m]  RMS = " + assess_res[1] + r" [$\mu$m]" + "\n"\
                + "Removed = " + assess_res[2] + r" [mm$^3$]" + "\n"\
                    + " with offset = " + assess_res[3] + r" [$\mu$m]" + "\n"\
                        + "( Ignore lower " + str(ignore_offset*100) + " %)"
        
        title_t = "\n" + title_res_str[option_wh] + "\n"\
            + "Black = WH 1-12 / Green = WH 13-24 / Purple = WH 25-36" + "\n"\
                + "Step RMS : " + str(round(force.std(), 2))
        
        fig = plt.figure(figsize=(10,22))
        ax_raw = image_plot(fig, title_raw, 421, raw, raw)
        ax_f = image_plot(fig, title_f, 422, filtered, raw)
        ax_d = image_plot(fig, title_d, 423, para_diff, raw)
        ax_i = image_plot(fig, title_i, 424, err_m, err_m)
        ax_rep = image_plot(fig, title_rep, 425, reprod, err_m)
        ax_res = image_plot(fig, title_res, 426, residual, residual)
        ax_torque = torque_plot(fig, title_t, 414, force)
        
        fig.tight_layout()
        #picname = mkfolder() + "r" + str(option_raw).zfill(2) + "_f" + str(option_filter) + str(filter_param) + "_id" + str(option_ideal) + "_wh" + str(option_wh) + ".png"
        picname = mkfolder() + "wh" + str(option_wh) + "_r" + str(option_raw).zfill(2) + "_id" + str(option_ideal) + "_f" + str(option_filter) + str(filter_param) + ".png"
        #picname = mkfolder() + "wh" + str(option_wh) + "_r" + str(option_raw).zfill(2) + "-rotate45_id" + str(option_ideal) + "_f" + str(option_filter) + str(filter_param) + ".png"
        
        fig.savefig(picname)
        
        print("\nReduced rate = ", 1 - float(assess_res[2])/float(assess_f[2]), "\n")
        if option_ideal == 0:
            fig.show()
        
        if option_ideal == 1:
            fig.clf()
        
        
        fig2 = plt.figure(figsize=(25, 5))
        gs = fig2.add_gridspec(1, 17)
        
        ax2_i = image_plot(fig2, "Target_shape", gs[0, 0:4], err_m, err_m)
        ax2_rep = image_plot(fig2, title_rep, gs[0, 5:9], reprod, err_m)
        ax2_torque = torque_plot(fig2, title_t, gs[0, 10:], force)
        fig.tight_layout()
    