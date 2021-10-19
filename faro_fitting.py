# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 15:34:36 2021

@author: swimc
"""

import numpy as np
import math
import pandas as pd
import scipy as sp
import time

import conbined_minimize as cb

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import mpl_toolkits.axes_grid1

def read_faro(fname, center):
    raw = pd.read_table(fname, header=None)
    ndarray = raw[[2,3,4]].values
    ndarray = ndarray.T
    ndarray[0] = ndarray[0] - center[0]
    ndarray[1] = ndarray[1] - center[1]
    return ndarray

def faro_interpolate(ndarray):
    x_old = ndarray[0]
    y_old = ndarray[1]
    c_old = ndarray[2]
    
    x_arr = y_arr = np.linspace(-m1_radi, m1_radi, px)
    xx_new, yy_new = np.meshgrid(x_arr, y_arr)
    
    xy_old = np.stack([x_old, y_old], axis=1)
    c_new = sp.interpolate.griddata(xy_old, c_old, (xx_new, yy_new), method="cubic", fill_value=np.nan)
    return c_new
   
def CircleFitting(x,y, sigma=False):
    """最小二乗法による円フィッティングをする関数
        input: x,y 円フィッティングする点群

        output  cxe 中心x座標
                cye 中心y座標
                re  半径

        参考
        元の関数
        https://myenigma.hatenablog.com/entry/2015/09/07/214600#Python%E3%82%B5%E3%83%B3%E3%83%97%E3%83%AB%E3%83%97%E3%83%AD%E3%82%B0%E3%83%A9%E3%83%A0
        一般式による最小二乗法（円の最小二乗法）　画像処理ソリューション
        http://imagingsolution.blog107.fc2.com/blog-entry-16.html
    """

    sumx  = sum(x)
    sumy  = sum(y)
    sumx2 = sum([ix ** 2 for ix in x])
    sumy2 = sum([iy ** 2 for iy in y])
    sumxy = sum([ix * iy for (ix,iy) in zip(x,y)])

    F = np.array([[sumx2,sumxy,sumx],
                  [sumxy,sumy2,sumy],
                  [sumx,sumy,len(x)]])

    G = np.array([[-sum([ix ** 3 + ix*iy **2 for (ix,iy) in zip(x,y)])],
                  [-sum([ix ** 2 *iy + iy **3 for (ix,iy) in zip(x,y)])],
                  [-sum([ix ** 2 + iy **2 for (ix,iy) in zip(x,y)])]])

    T=np.linalg.inv(F).dot(G)

    cxe=float(T[0]/-2)
    cye=float(T[1]/-2)
    re=math.sqrt(cxe**2+cye**2-T[2])
    #print (cxe,cye,re)
    if sigma==False:    
        return (cxe,cye,re)
    if sigma==True:
        sigma = sum(( (x - cxe)**2 + (y - cye)**2 - re**2 ) **2 )
        return (cxe,cye,re, sigma)
    
def rotation_func(X, Param):
    rot1, rot2, rot3 = X
    edge = Param
    rot = sp.spatial.transform.Rotation.from_rotvec([rot1, rot2, rot3])
    surf_rot = rot.apply(edge.T).T
    x_rot, y_rot, z_rot = surf_rot
    cxe, cye, re, sigma = CircleFitting(x_rot, y_rot, True)
    return sigma
    
def RotOptimize(edge):
    result_rot = sp.optimize.minimize(rotation_func, x0=(0,0,0), 
                                      args=(edge,), method="Powell")

    rot = sp.spatial.transform.Rotation.from_rotvec(result_rot.x)
    edge_rot = rot.apply(edge.T).T
    x_rot, y_rot, z_rot = edge_rot
    cxe, cye, re = CircleFitting(x_rot, y_rot)
    theta = np.linspace(0, 2*np.pi, 1000)
    opt_circle = [re * np.cos(theta) + cxe, re * np.sin(theta) + cye]    
    result_dict = {
        }
    return result_rot, edge_rot, (cxe, cye, re), opt_circle
    
    


def minimize_func(X, Param):
    o, p, q, z= X
    surf, x, y = Param
    ideal = cb.cz1(o, p, q, x, y) + z
    return np.nansum(abs(surf - ideal))

def cz1_result(OptimizeResult, x, y):
    o, p, q, z = OptimizeResult["x"]
    return cb.cz1(o, p, q, x, y) + z

def Optimize(surf):
    x = surf[0]
    y = surf[1]
    param = (surf, x, y)
    result = sp.optimize.minimize(minimize_func, x0=(0.1,8667,1800, 0), 
                                   args=(param,), method="Powell")
    
    surf_diff = [x, y, surf[2] - cz1_result(result, x, y)]
    return result, surf_diff
 
    
    
    
    

def scatter2d(fig, title, position, ndarray):
    cmap = cm.jet
    fs = 15
    
    x = ndarray[0]
    y = ndarray[1]
    c = ndarray[2]
    
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

if __name__ == '__main__':
    px = 512
    m1_radi = 1850/2
    ignore_radi = 200
    valid_radi = m1_radi - ignore_radi # 有効半径
        
    xx, yy = np.meshgrid(np.linspace(-m1_radi, m1_radi, px),np.linspace(-m1_radi, m1_radi, px))
    mask = np.where(xx**2+yy**2<valid_radi**2, 1, np.nan)
        
    
    side_raw1 = read_faro("_PLANETS_faro/鏡側面（円筒）.txt",[0,0])
    
    result_rot1, side1_rot, circleparam1, euler1, opt_circle1 = RotOptimize(side_raw1)
    
    surf1_raw = read_faro("_PLANETS_faro/鏡面上1回目.txt", center_fit)
    surf2_raw = read_faro("_PLANETS_faro/鏡面上（アキシャルパッド上）.txt", center_fit)
    
    
    surf1_pre = surf1_raw
    surf1_pre[2] = surf1_raw[2] - surf1_raw[2].min()
    surf2_pre = surf2_raw
    surf2_pre[2] = surf2_raw[2] - surf2_raw[2].min()
    surf3_pre = np.concatenate([surf1_pre, surf2_pre], axis=1)
    
    start_time = time.time()
    result1, surf1_diff = Optimize(surf1_pre)
    result2, surf2_diff = Optimize(surf2_pre)
    result3, surf3_diff = Optimize(surf3_pre)
    
    surf1_min = mask * cz1_result(result1, xx, yy)
    surf2_min = mask * cz1_result(result2, xx, yy)
    surf3_min = mask * cz1_result(result3, xx, yy)
    
    print(str(time.time() - start_time))
    ## for plot
    fig3 = plt.figure(figsize=(5,10))
    gs3 = fig3.add_gridspec(2,1)
    ax_side1 = scatter2d(fig3, "side", gs3[0,0], side_raw1)
    ax_side1.scatter([-1000,-1000, 1000, 1000],[1000,-1000, 1000, -1000])
    theta = np.linspace(0, 2*np.pi, 1000)
    x_sideplot = center_fit[2] * np.cos(theta)
    y_sideplot = center_fit[2] * np.sin(theta)
    ax_side2 = scatter2d(fig3, "side", gs3[1,0], side_raw2)
    ax_side2.scatter(x_sideplot, y_sideplot, s=1)
    
    fig7 = plt.figure(figsize = (5,10))
    gs7 = fig7.add_gridspec(2,1)
    ax7_side = scatter2d(fig7, "side 1", gs7[0,0], side1_rot)
    ax7_side.scatter(opt_circle1[0], opt_circle1[1], s=1)

    fig1 = plt.figure(figsize = (5,10))
    gs1 = fig1.add_gridspec(2,1)
    ax_surf1 = scatter2d(fig1, "surf 1", gs1[0,0], surf1_raw)
    ax_surf2 = scatter2d(fig1, "surf axialpad", gs1[1,0], surf2_raw)
    
    
    """
    fig2 = plt.figure(figsize = (5,10))
    gs2 = fig2.add_gridspec(2,1)
    ax_fit1 = cb.image_plot(fig2, "surf 1 spline", gs2[0,0], surf1_fit, surf1_fit, cb_micron=False)
    ax_fit2 = cb.image_plot(fig2, "surf axialpad spline", gs2[1,0], surf2_fit, surf2_fit, cb_micron=False)
    """
    fig4 = plt.figure(figsize=(5, 15))
    gs4 = fig4.add_gridspec(3,1)
    ax_minimize1 = cb.image_plot(fig4, "fitted", gs4[0,0], surf1_min, surf1_min, cb_micron=False)
    ax_minimize2 = cb.image_plot(fig4, "axialpad fitted", gs4[1,0], surf2_min, surf2_min, cb_micron=False)
    ax_minimize3 = cb.image_plot(fig4, "fitted 1 and 2", gs4[2,0], surf3_min, surf3_min, cb_micron=False)
    
    """
    fig5 = plt.figure(figsize=(5, 10))
    gs5 = fig5.add_gridspec(2,1)
    ax_diff1 = cb.image_plot(fig5, "diff", gs5[0,0], surf1_min - ideal, surf1_min - ideal, cb_micron=False)
    ax_diff2 = cb.image_plot(fig5, "diff", gs5[1,0], surf2_min - ideal, surf2_min - ideal, cb_micron=False)
    """
    
    fig6 = plt.figure(figsize=(5,15))
    gs6 = fig6.add_gridspec(3,1)
    ax6_1 = scatter2d(fig6, "", gs6[0,0], surf1_diff)
    ax6_2 = scatter2d(fig6, "", gs6[1,0], surf2_diff)
    ax6_3 = scatter2d(fig6, "", gs6[2,0], surf3_diff)
    