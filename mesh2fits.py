# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 17:06:43 2021

@author: swimc
"""

import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import PIL
import datetime

import stitch2mesh as s2m

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

def image_rotate(mesh, degree):
    im = PIL.Image.fromarray(mesh)
    im_rotate = im.rotate(degree)
    return im_rotate

if __name__ == '__main__':
    m1_radi = 1850/2
    varid_radi = 900
    px_new = 1001
    fits_radi = 1000 # 1pixel = 2mm で、1001*1001px の fits作成
    rotate_angle = 10.814 #+ 0.979 # [deg] 今回の計測はTT座標系基準のデータなので、TT座標系と光軸方向のずれは研磨用fitsに入れる必要は無い
    
    fname_rawdata = "mkfolder/stitch2mesh/zer10_1116xm130ym830WTp03.v3.hei_dense.csv"
    rawdata = np.loadtxt(fname_rawdata)
    
    ## 回転 ------------------------------------------------------------------
    tf_old = ~np.isnan(rawdata)
    mask_old = np.where(tf_old==True, 1, np.nan)
    zz_rotate = mask_old * image_rotate(rawdata*tf_old, -1 * rotate_angle)
    zz_old = np.where(~np.isnan(zz_rotate)==True, zz_rotate, 0)
    
    ## px変更 ----------------------------------------------------------------
    px_old = len(rawdata)
    xx_old, yy_old = np.meshgrid(np.linspace(-m1_radi, m1_radi, px_old),np.linspace(-m1_radi, m1_radi, px_old))
    
    old = np.stack([xx_old.flatten(), yy_old.flatten(), zz_old.flatten()], axis=1)
    
    xx_new, yy_new = np.meshgrid(np.linspace(-fits_radi, fits_radi, px_new),np.linspace(-fits_radi, fits_radi, px_new))
    tf_new = np.where(xx_new**2+yy_new**2 < varid_radi**2, True, False)
    mask_new = np.where(tf_new==True, 1, np.nan)
    zz_new = mask_new * s2m.interpolate(old, xx_new, yy_new)
    
    
    ## fits書き込み -----------------------------------------------------------
    fname_header = "raw_data/0923xm130_1007ym830.hei.v2_dense.z3.fits/0923xm130_1007ym830.hei.v2_dense.z3.fits"
    fname_txtfile = fname_rawdata[27:-4]
    
    hdu = fits.open(fname_header)[0]
    
    hdu.header["file_den"] = fname_txtfile + ".txt"
    hdu.header["vertex_a"] = "{:.04f}".format(rotate_angle)
    hdu.header["tsmpl1"] = datetime.datetime.now().strftime("%Y%m%d %H%M%S")
    hdu.header["tmskpl1"] = datetime.datetime.now().strftime("%Y%m%d %H%M%S")
    
    hdu.data = (zz_new * 1e6).astype("float32") # mm -> nm
    hdu.writeto(mkfolder() + fname_txtfile + ".fits", overwrite=True)
    
    