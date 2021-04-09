# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 22:45:51 2021

@author: swimc
"""

import numpy as np

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

def cz1_theta(o, p, aqx, cx, cy, theta_deg):
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
    theta_deg : [degree]

    Returns
    -------
    cz1 : 2D-mesh-array [mm]
        input された o, p, q での切り取り放物面
    """

    phi = o / 1000 # [mrad] -> [rad] への換算

    a = 1/(2*p)
    #theta = np.arctan(2*a*aqx)
    theta = np.deg2rad(theta_deg)
    aqz = a * ( aqx**2 + 0**2 )

    D = -4*a**2*cx**2*np.sin(phi)**2*np.sin(theta)**2 - 8*a**2*cx*cy*np.sin(phi)*np.sin(theta)**2*np.cos(phi) + 4*a**2*cy**2*np.sin(phi)**2*np.sin(theta)**2 - 4*a**2*cy**2*np.sin(theta)**2 + 2*a*aqx*np.sin(2*theta) + 4*a*aqz*np.sin(theta)**2 + 4*a*cx*np.sin(theta)*np.cos(phi) - 4*a*cy*np.sin(phi)*np.sin(theta) - np.sin(theta)**2 + 1

    cz1 = (4*a*aqx*np.sin(theta) - a*cx*(np.sin(phi - 2*theta) - np.sin(phi + 2*theta)) - a*cy*(np.cos(phi - 2*theta) - np.cos(phi + 2*theta)) - 2*np.sqrt(D) + 2*np.cos(theta))/(4*a*np.sin(theta)**2)
    return cz1

def solve(a, b, c):
    D = b**2 - 4 * a * c
    x = (-b + np.sqrt(D)) / (2 * a)
    return x

def hanamura_z(theta_deg, p, x, y):
    theta = np.deg2rad(theta_deg)

    x0 = p * np.tan(theta)
    y0 = (p/2) * np.tan(theta)**2

    X = x * np.cos(theta) + x0
    Z = x * np.sin(theta) + y0

    a = 1 / (2 * p)
    b = np.tan( np.pi / 2 - theta)
    c = a * y**2 - Z - b * X
    nodex = solve(a, b, c)
    z = (X - nodex) / np.sin(theta)
    return z

if __name__ == '__main__':
    # parametar --------------------------------------------------------------

    # o 回転角[mrad], p 曲率半径[mm], q 軸外し距離[mm]
    o0, do = 0, 0
    p0, dp = 8667, 0
    q0, dq = 1800, 1

    # constant ---------------------------------------------------------------
    m1_radi = 1850/2
    """
    print("deg is 11.89")
    print("r = ", m1_radi, "\n")
    print("0 kn : ", cz1(0, p0, 1800, m1_radi, 0))
    print("0 yh : ", hanamura_z(np.rad2deg(np.arctan(q0/p0)), p0, m1_radi, 0), "\n")

    print("6 kn : ", cz1(0, p0, 1800, -m1_radi, 0))
    print("6 yh : ", hanamura_z(np.rad2deg(np.arctan(q0/p0)), p0, -m1_radi, 0), "\n")

    print("3 kn : ", cz1(0, p0, 1800, 0, m1_radi))
    print("3 yh : ", hanamura_z(np.rad2deg(np.arctan(q0/p0)), p0, 0, m1_radi), "\n")

    print("9 kn : ", cz1(0, p0, 1800, 0, -m1_radi))
    print("9 yh : ", hanamura_z(np.rad2deg(np.arctan(q0/p0)), p0, 0, -m1_radi), "\n")

    print("deg is 11.73")
    print("r = ", m1_radi, "\n")
    print("0 kn : ", cz1_theta(0, p0, 1800, m1_radi, 0, 11.73))
    print("0 yh : ", hanamura_z(11.73, p0, m1_radi, 0), "\n")

    print("6 kn : ", cz1_theta(0, p0, 1800, -m1_radi, 0, 11.73))
    print("6 yh : ", hanamura_z(11.73, p0, -m1_radi, 0), "\n")

    print("3 kn : ", cz1_theta(0, p0, 1800, 0, m1_radi, 11.73))
    print("3 yh : ", hanamura_z(11.73, p0, 0, m1_radi), "\n")

    print("9 kn : ", cz1_theta(0, p0, 1800, 0, -m1_radi, 11.73))
    print("9 yh : ", hanamura_z(11.73, p0, 0, -m1_radi), "\n")
    """
    print("q=",q0)
    print("0 kn : ", cz1(0, p0, q0, m1_radi, 0))
    print("0 yh : ", hanamura_z(np.rad2deg(np.arctan(q0/p0)), p0, m1_radi, 0), "\n")

    print("q=",q0+dq)
    print("0 kn : ", cz1(0, p0, q0+dq, m1_radi, 0))
    print("0 yh : ", hanamura_z(np.rad2deg(np.arctan((q0+dq)/p0)), p0, m1_radi, 0), "\n")
