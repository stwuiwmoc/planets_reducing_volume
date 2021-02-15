# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 17:24:17 2021

@author: swimc
"""
import time

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

def drop_recovery(array_drop, mask_2d, pixel):
    mask_1d = mask_2d.flatten()
    loc = np.where(np.isnan(mask_1d)==True, np.nan, np.arange(mask_1d.size))
    loc_1d = loc.flatten()
    loc_drop = loc_1d[~np.isnan(loc_1d)]
    
    recovery_1d = np.empty(px**2) * np.nan
    for i in range(loc_drop.size):
        idx = int(loc_drop[i])
        recovery_1d[idx] = array_drop[i]
    recovery = recovery_1d.reshape((pixel,pixel))
    return recovery

if __name__ == '__main__':

    """
    px = 256
    m1_radi = 1850/2
    raw = np.linspace(100, 500, px**2).reshape((px, px))
    
    x_arr = y_arr = np.linspace(-m1_radi, m1_radi, px)
    xx, yy = np.meshgrid(x_arr, y_arr)
    
    raw = np.where(xx**2+yy**2<=m1_radi**2, raw, np.nan)
    raw_1d = raw.flatten()
    
    loc = np.where(np.isnan(raw_1d)==True, np.nan, np.arange(raw_1d.size))
    
    raw_drop = raw_1d[~np.isnan(raw_1d)]
    loc_drop = loc[~np.isnan(raw_1d)]
    
    raw_calc = raw_drop + np.arange(raw_drop.size)
    
    raw_re_1d = np.empty(px**2) * np.nan
    start = time.time()
    
    for i in range(loc_drop.size):
        idx = int(loc_drop[i])
        raw_re_1d[idx] = raw_drop[i]
    
    raw_re = raw_re_1d.reshape(px,px)
    end = time.time()
    print(end-start)
    """
    
    # 作用行列用のテスト
    file_num = 33
    px = 256
    m1_radi = 1850/2    
    zer_order = 21
    
    x_arr = y_arr = np.linspace(-m1_radi, m1_radi, px)
    xx, yy = np.meshgrid(x_arr, y_arr)
    
    mask = np.where(xx**2+yy**2<=m1_radi**2, 1, np.nan)
    tf = ~np.isnan(mask)
    
    df0 = read("_Fxx/PM3.5_36ptAxWT03_F00.smesh.txt")
    
    dfxx = read("_Fxx/PM3.5_36ptAxWT03_F01.smesh.txt")
    
    diff = mask * kriging(df0, dfxx)
    
    diff_drop, loc_drop = nan_drop(diff, mask)
    
    """
    diff_re_1d = np.empty(px**2) * np.nan
    for i in range(loc_drop.size):
        idx = int(loc_drop[i])
        diff_re_1d[idx] = diff_drop[i]
    
    diff_re = diff_re_1d.reshape((px,px))
    """
    diff_re = drop_recovery(diff_drop, mask)