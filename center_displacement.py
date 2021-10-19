# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 11:18:36 2021

@author: swimc
cz1の計算だとoffaxis距離が実物の鏡（真円）の中心になっており、鏡の0時方向（＋ｘ）と6時方向（ーｘ）の端点（ｒ=925）での
ｚ座標が異なる

しかし実物の鏡では最もへこみが大きい部分は鏡の円の中心から少しずれており、また、0時6時方向の端点のz方向高さはほぼ同じ
接点がどの程度鏡の円の中心からずれれば0時6時方向のｚ方向高さが一致するのかを計算

"""

import numpy as np
import matplotlib.pyplot as plt
import conbined_minimize as cb
import pandas as pd

import conbined_minimize as cb

def mkfolder(suffix = ""):
    import os
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
    theta = np.arctan(2*a*aqx)
    aqz = a * ( aqx**2 + 0**2 )

    D = -4*a**2*cx**2*np.sin(phi)**2*np.sin(theta)**2 - 8*a**2*cx*cy*np.sin(phi)*np.sin(theta)**2*np.cos(phi) + 4*a**2*cy**2*np.sin(phi)**2*np.sin(theta)**2 - 4*a**2*cy**2*np.sin(theta)**2 + 2*a*aqx*np.sin(2*theta) + 4*a*aqz*np.sin(theta)**2 + 4*a*cx*np.sin(theta)*np.cos(phi) - 4*a*cy*np.sin(phi)*np.sin(theta) - np.sin(theta)**2 + 1

    cz1 = (4*a*aqx*np.sin(theta) - a*cx*(np.sin(phi - 2*theta) - np.sin(phi + 2*theta)) - a*cy*(np.cos(phi - 2*theta) - np.cos(phi + 2*theta)) - 2*np.sqrt(D) + 2*np.cos(theta))/(4*a*np.sin(theta)**2)
    return cz1

def calc_tilt(p, q):
    """
    cz1の傾き計算部分のみ

    Parameters
    ----------
    p : float [mm]
        曲率半径
    q : floot [mm]
        軸外し距離 (off-axis)
    
    Returns
    -------
    tilt [rad]
    None.

    """
    a = 1/(2*p)
    theta = np.arctan(2*a*q)
    return theta
    

if __name__ == '__main__':
    m1_radi = 925
    
    p0 = 8667
    q0 = 1800
    o0 = 0
    
    dq = np.arange(-20, 20, 0.1)
    #dq = 10
    
    theta_new = calc_tilt(p0, q0+dq)
    s = dq / np.cos(theta_new)
    
    dir_0 = cz1(o0, p0, q0+dq, m1_radi-s, 0)
    dir_6 = cz1(o0, p0, q0+dq, -m1_radi-s, 0)
    dir_diff = dir_0 - dir_6
    diff_min_idx = np.argmin(abs(dir_diff))
    diff_min_dq = dq[diff_min_idx]
    diff_min_s = s[diff_min_idx]
    print("dq ", diff_min_dq)
    print("s", diff_min_s)
    
    fig_dq = plt.figure(figsize=(7,7))
    gsdq = fig_dq.add_gridspec(2, 1)
    axdq_dir = fig_dq.add_subplot(gsdq[0,0])
    axdq_dir.plot(dq, dir_0, label="dir_0")
    axdq_dir.plot(dq, dir_6, label="dir_6")
    axdq_dir.grid()
    axdq_dir.legend()
    
    axdq_diff = fig_dq.add_subplot(gsdq[1,0])
    axdq_diff.plot(dq, dir_diff, label="dir_0 - dir_6")
    axdq_diff.grid()
    axdq_diff.vlines(diff_min_dq, dq.min(), dq.max())
    axdq_diff.set_xlabel("dq [mm]")
    
    fig_s = plt.figure(figsize=(7,7))
    gss = fig_s.add_gridspec(2, 1)
    axs_dir = fig_s.add_subplot(gss[0,0])
    axs_dir.plot(s, dir_0)
    axs_dir.plot(s, dir_6)
    axs_dir.grid()

    axs_diff = fig_s.add_subplot(gss[1,0])
    axs_diff.plot(s, dir_diff)
    axs_diff.grid()    
    axs_diff.vlines(diff_min_s, s.min(), s.max())
    axs_diff.set_xlabel("s [mm]")
    
    
    ## figure check
    px = 1023
    xx, yy = np.meshgrid(np.linspace(-m1_radi, m1_radi, px),np.linspace(-m1_radi, m1_radi, px))
    outer_mask = np.where(xx**2 + yy**2 >m1_radi**2, np.nan, 1)
    inner_mask = np.where(xx**2 + yy**2 <(m1_radi-200)**2, np.nan, 1)
    
    ideal = outer_mask * inner_mask * cb.cz1(o0, p0, q0, xx, yy)
    displace = outer_mask * inner_mask * cb.cz1(o0, p0, q0+diff_min_dq, xx-diff_min_s, yy)
    
    fig = plt.figure(figsize=(7,14))
    ax1 = cb.image_plot(fig, "ideal", 211, ideal, ideal)
    ax2 = cb.image_plot(fig, "displace", 212, displace, displace)
    
    x_arr = np.linspace(-m1_radi, m1_radi, px)
    outer_mask_1d = np.where(x_arr**2 > m1_radi**2, np.nan, 1)
    ideal_1d = outer_mask_1d * cb.cz1(o0, p0, q0, x_arr, 0)
    displace_1d = outer_mask_1d * cb.cz1(o0, p0, q0+diff_min_dq, x_arr-diff_min_s, 0)
    fig1d = plt.figure(figsize = (7,7))
    ax1d1 = fig1d.add_subplot(211)
    ax1d1.plot(x_arr, ideal_1d, label="ideal")
    ax1d1.plot(x_arr, displace_1d, label="displace")
    ax1d1.grid()
    ax1d1.legend()
    ax1d1.set_title("x-z plot (y = 0)")
    
    ax1d2 = fig1d.add_subplot(212)
    ax1d2.plot(x_arr, ideal_1d, label="ideal")
    ax1d2.plot(x_arr, displace_1d, label="displace")
    ax1d2.set_ylim(40, 48)
    ax1d2.grid()
    ax1d2.legend()