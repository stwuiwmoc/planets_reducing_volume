# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 19:52:48 2021

@author: swimc
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
import PIL

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


if __name__ == '__main__':
    file_num = 33
    px = 256
    m1_radi = 1850/2    
    zer_order = 21
    angle = 240
    
    x_arr = y_arr = np.linspace(-m1_radi, m1_radi, px)
    xx, yy = np.meshgrid(x_arr, y_arr)
    
    mask = np.where(xx**2+yy**2<=m1_radi**2, 1, np.nan)
    tf = ~np.isnan(mask)
    
    i = 1
    
    num_i = str(i).zfill(2)
    data_fname_raw = "_Fxx/PM3.5_36ptAxWT02_F" + num_i + ".smesh.txt"
    dfxx_raw = read(data_fname_raw)
    df0 = read("_Fxx/PM3.5_36ptAxWT02_F00.smesh.txt")
    
    surf_raw = mask * kriging(df0, dfxx_raw)
    img = PIL.Image.fromarray(surf_raw)
    img_120 = img.rotate(angle)
    surf_raw_120 = mask * np.asarray(img_120)
    
    diff_arr = np.empty(file_num)
    
    for j in range(0, file_num):
        num_j = str(j+1).zfill(2)
        
        data_fname = "_Fxx/PM3.5_36ptAxWT02_F" + num_j + ".smesh.txt"
        dfxx = read(data_fname)
        
        surf = mask * kriging(df0, dfxx)
        
        diff = abs(surf - surf_raw_120)
        
        diff_sum = np.nansum(diff)
        
        diff_arr[j] = diff_sum
        
        print(num_j, diff_sum)
        
    j_min = np.argmin(diff_arr)
    num_min = str(j_min+1).zfill(2)
    df_min = read("_Fxx/PM3.5_36ptAxWT02_F" + num_min + ".smesh.txt")
    surf_min = mask * kriging(df0, df_min)
    
    
    fig = plt.figure(figsize=(10,10))
    ax_raw = image_plot(fig, "raw " + num_i, 221, surf_raw, surf_raw, False)
    ax_120 = image_plot(fig, str(angle), 222, surf_raw_120, surf_raw, False)
    ax_min = image_plot(fig, num_min, 223, surf_min, surf_raw, False)
    ax_diff = image_plot(fig, "difference", 224, surf_min-surf_raw, surf_raw, False)
    
    fig.tight_layout()
    
    picname = mkfolder() + "/" + num_i + "_" + str(angle) + ".png"
    fig.savefig(picname)