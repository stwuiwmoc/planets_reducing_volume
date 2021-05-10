# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 16:55:29 2021

@author: swimc
"""

import numpy as np

import conbined_minimize as mymin

def tilt_calc(r_mm, theta_deg):
    theta_rad = np.deg2rad(theta_deg)
    
    x1 = r_mm * np.cos(theta_rad)
    y1 = r_mm * np.sin(theta_rad)
    
    x2 = (r_mm + dr) * np.cos(theta_rad)
    y2 = (r_mm + dr) * np.sin(theta_rad)
    
    z1 = mymin.cz1(0, p0, q0, x1, y1)
    z2 = mymin.cz1(0, p0, q0, x2, y2)
    
    a = (z2 - z1) / dr
    
    tilt_deg = np.rad2deg(np.arctan(a))
    return tilt_deg

if __name__ == '__main__':
    
    theta0 = 0 # 鏡の現実の角度を補正する [deg]
    dr = 0.01 # 傾きを計算するための動径方向の微小距離 [mm]
    
    p0 = 8667
    q0 = 1800
    
    #r = float(input("Input radius [mm] -> "))
    #theta = float(input("Input theta [deg] -> ")) + theta0
    
    tilt_center = tilt_calc(0, 0)
    tilt_limb_xplus = tilt_calc(900, 0)
    tilt_limb_xminus = tilt_calc(900, 180)