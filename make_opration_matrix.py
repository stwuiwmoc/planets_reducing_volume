# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 19:24:38 2021

@author: swimc

wh による33個の単位駆動を それぞれzernike fittingして、zernike 作用行列を作成
256^2 * 33 の作用行列を作成
それぞれの作用行列をcsvに保存

f01 - f33 について、もとの691点、krigingした元データ、zernike fittingで再現した形状、各項のrmsをplot
"""

import numpy as np
import pandas as pd
import os
import pykrige as krg
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

def read(filename):
    #skip行数の設定
    skip = 0
    with open(filename) as f:
        while True:
            line = f.readline()
            if line[0] == '%':
                skip += 1
            else:
                break
            #エラーの処理
            if len(line) == 0:
                break

    #データの読み出し
    df = pd.read_csv(filename, sep='\s+', skiprows=skip, header=None) #\s+...スペース数に関わらず区切る
    df.columns = ["x", "y", "z", "color", "dx", "dy", "dz" ]
    df = df * 10**3 # [m] -> [mm]
    return df

def kriging(df_0, dfxx):
    # input [mm] -> output [mm]    
    #mesh型のデータを格子点gridに補完
    df_0 = df_0
    dfxx = dfxx
    
    x = dfxx["x"]
    y = dfxx["y"]
    dw = dfxx["dz"] - df_0["dz"]
    
    ok_module = krg.ok.OrdinaryKriging(x, y, dw)
    z, sigmasq = ok_module.execute("grid", x_arr, y_arr) #格子点にdwをfitting
    return z

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
        norm = Normalize(vmin=np.nanmin(c_scale*10**6), vmax=np.nanmax(c_scale*10**6))
        cbar_title = "[nm]"
    mappable = cm.ScalarMappable(norm = norm, cmap = cm.jet)
    cbar = fig.colorbar(mappable, ax=ax, cax=cax)
    cbar.set_label(cbar_title, fontsize=fs)
        
    return ax

def df_plot(fig, title, position, df0, dfxx):
    cmap = cm.jet
    fs = 15
    
    x = df0["x"]
    y = df0["y"]
    c = (dfxx["dz"] - df0["dz"]) * 10**6 # [mm] -> [nm]
    
    ax = fig.add_subplot(position)
    ax.scatter(x, y, c=c, cmap=cm.jet)
    ax.set_title(title, fontsize=fs)
    
    divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes('right', '5%', pad='3%')
    norm = Normalize(vmin=np.nanmin(c), vmax=np.nanmax(c))
    cbar_title = "[nm]"
    
    mappable = cm.ScalarMappable(norm = norm, cmap = cm.jet)
    cbar = fig.colorbar(mappable, ax=ax, cax=cax)
    cbar.set_label(cbar_title, fontsize=fs)
    return ax

def zer_term_plot(fig, title, position, zer):
    # input zer [m], plot zer[nm]
    ax = fig.add_subplot(position)
    fs = 15
    
    x = np.arange(1, zer_order+1)
    ax.plot(x, zer*10**9, marker="s", linewidth=1)
    ax.grid()
    ax.set_xlabel("zernike order", fontsize=fs)
    ax.set_xticks(x)
    ax.set_ylabel("zernike RMS [nm]", fontsize=fs)
    ax.set_title(title, fontsize=fs)
    return ax

if __name__ == '__main__':
    file_num = 33
    px = 256
    m1_radi = 1850/2    
    zer_order = 21
    
    x_arr = y_arr = np.linspace(-m1_radi, m1_radi, px)
    xx, yy = np.meshgrid(x_arr, y_arr)
    
    mask = np.where(xx**2+yy**2<=m1_radi**2, 1, np.nan)
    tf = ~np.isnan(mask)
    
    zer_opration_matrix = np.empty((file_num, zer_order))
    opration_matrix = np.empty((file_num, tf.sum()))
    
    df0 = read("_Fxx/PM3.5_36ptAxWT03_F00.smesh.txt")

    for i in range(0, file_num):
        
        num = str(i+1).zfill(2)
        
        data_fname = "_Fxx/PM3.5_36ptAxWT03_F" + num + ".smesh.txt"
        dfxx = read(data_fname)
        
        diff = tf * kriging(df0, dfxx) / 1000 # [mm] -> [m]
        diff_drop, loc_drop = nan_drop(diff, mask)
        opration_matrix[i] = diff_drop
        
        zer, fit = pr.prop_fit_zernikes(diff, tf, px/2, zer_order, xc=px/2, yc=px/2, FIT=True) 
        zer_opration_matrix[i] = zer
        
        diff = mask * diff * 1000 # [m] -> [mm]
        fit = mask * fit * 1000 # [m] -> [mm]

        fig = plt.figure(figsize=(11,10))
        ax_df = df_plot(fig, "raw : " + num, 221, df0, dfxx)        
        ax_dz = image_plot(fig, "Kriging", 222, diff, diff, False)
        ax_zer = zer_term_plot(fig, "zernike terms", 223, zer)
        ax_fit = image_plot(fig, "zernike fitting", 224, fit, diff, False)
        
        """
        fig = plt.figure(figsize=(8,7))
        ax_for_ppt = image_plot(fig, "F"+ num, 111, mask*diff, mask*diff, False)
        fig.tight_layout()
        """
        
        picname = mkfolder() + "WT03_F" + num + ".png"
        fig.savefig(picname)
        fig.clf()
        
        print(num)
    
    save_fname = "WT03_256_opration_matrix[m].csv"
    np.savetxt(save_fname, opration_matrix, delimiter=",")
    
    zer_save_fname = "WT03_zer_opration_matrix[m].csv"
    np.savetxt(zer_save_fname, zer_opration_matrix, delimiter=",")
    """
    """