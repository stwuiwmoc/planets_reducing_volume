# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 10:24:50 2021

@author: swimc
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

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
    import os
    filename = os.path.basename(__file__)
    filename = filename.replace(".py", "") + suffix
    folder = "mkfolder/" + filename + "/"
    os.makedirs(folder, exist_ok=True)
    return folder

def interpolate(raw, x_mesh, y_mesh):
    x_old = raw[:,0]
    y_old = raw[:,1]
    dz_old = raw[:,2]
    
    xy_old = np.stack([x_old, y_old], axis=1)
    dz_new = sp.interpolate.griddata(xy_old, dz_old, (x_mesh, y_mesh), method="linear", fill_value = 0)
    return dz_new

if __name__ == "__main__":
    m1_radi = 1850/2
    ignore_radi = 25
    px = 1023
    varid_radi = m1_radi - ignore_radi
    
    fname = "raw_data/diff.txt"
    raw = np.loadtxt(fname)
    
    xx, yy = np.meshgrid(np.linspace(-m1_radi, m1_radi, px),np.linspace(-m1_radi, m1_radi, px))
    z_mesh = interpolate(raw, xx, yy)
    
    plt.imshow(z_mesh)
    
    save_fname = mkfolder() + "diff.csv"
    np.savetxt(save_fname, z_mesh)
    