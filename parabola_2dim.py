# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 16:22:02 2021

@author: swimc
"""

import numpy as np
import matplotlib.pyplot as plt

def line(a, x0):
    dydx = 2 * a * x0
    y0 = a * x0**2
    angle = np.arctan(dydx)
    x_arr = np.arange(x0-925*np.cos(angle), x0+925*np.cos(angle), 0.1)
    
    y_line = dydx * (x_arr - x0) + y0
    return x_arr, y_line


if __name__ == '__main__':
    a0 = 1/(2*8667)
    x = np.arange(0, 7000, 0.1)
    y = a0 * x**2
    
    q0 = 1800
    q1 = 6000
    
    fig1 = plt.figure(figsize=(25,25))
    ax_q = fig1.add_subplot()
    ax_q.plot(x, y)
    ax_q.plot(line(a0, q0)[0], line(a0, q0)[1], linewidth=3)
    ax_q.plot(line(a0, q1)[0], line(a0, q1)[1], linewidth=3)
    
    ax_q.vlines(q0, ymin=y.min(), ymax=y.max(), color="darkgray", linestyles="dashdot")
    ax_q.vlines(q1, ymin=y.min(), ymax=y.max(), color="darkgray", linestyles="dashdot")
    ax_q.set_aspect("equal")
    
    ax_q.set_xlabel("y (Off-Axis Distance) [mm]")
    
    x1 = np.arange(0, 3000, 0.1)
    a1 = 1/(2 * 20000)
    fig2 = plt.figure(figsize=(25,25))
    ax_p = fig2.add_subplot()
    
    ax_p.plot(x1, a0*x1**2)
    ax_p.plot(x1, a1*x1**2)
    ax_p.plot(line(a0, q0)[0], line(a0, q0)[1], linewidth=3)
    ax_p.plot(line(a1, q0)[0], line(a1, q0)[1], linewidth=3)
    
    ax_p.vlines(q0, ymin=(a0*x1**2).min(), ymax=(a0*x1**2).max(), color="darkgray", linestyles="dashdot")
    
    ax_p.set_aspect("equal")
    
    
    