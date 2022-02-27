# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 10:24:50 2021

@author: swimc

stitch後のデータをmeshに変換し、zernike低次項を除く
"""
# %%
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import proper as pr

import conbined_minimize as cb


def mkfolder(suffix=""):
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
    x_old = raw[:, 0]
    y_old = raw[:, 1]
    dz_old = raw[:, 2]

    xy_old = np.stack([x_old, y_old], axis=1)
    dz_new = sp.interpolate.griddata(
        xy_old, dz_old, (x_mesh, y_mesh), method="cubic", fill_value=np.nan)
    return dz_new


def zer2mesh(zernike_arr, zernike_order, radi, pixel):
    # input zer_arr [m] -> output 2d-mesh [mm]
    wavelength = 500e-9
    diam = 2 * radi * 1e3  # m 単位の直径に直す
    zernike_num = np.arange(1, zernike_order + 1)

    wavestruct = pr.prop_begin(diam, wavelength, pixel, 1)
    wfe = pr.prop_zernikes(wavestruct, zernike_num, zernike_arr)
    wfe = wfe * 1e3  # output mm
    return wfe


if __name__ == "__main__":
    m1_radi = 1850 / 2
    ignore_radi = 25
    px = 1023

    varid_radi = m1_radi - ignore_radi
    zer_fit_order = 10

    fname = "raw_data/0207xm130All2_0208cir_0208ykagomeAllCc.v4.8.hei_dense.txt"
    shaped = np.loadtxt(fname)[:, 1:4]

    xx, yy = np.meshgrid(
        np.linspace(-m1_radi, m1_radi, px),
        np.linspace(-m1_radi, m1_radi, px))

    z_mesh = interpolate(shaped, xx, yy) * 1e-6  # mm単位に合わせる

    tf = ~np.isnan(z_mesh)
    mask = np.where(tf, 1, np.nan)

    zer, fit = pr.prop_fit_zernikes(
        tf * z_mesh * 1e-3,  # properでの指定は単位 [m]
        tf, px / 2, zer_fit_order, xc=px / 2, yc=px / 2, FIT=True)

    # mask_varid = np.where(xx**2+yy**2<varid_radi**2, True, np.nan)

    zer10mesh = zer2mesh(zer, zer_fit_order, m1_radi, px)
    diff_zer10 = (z_mesh - zer10mesh) * mask

    zer03mesh = zer2mesh(zer[:3], 3, m1_radi, px)
    diff_zer03 = (z_mesh - zer03mesh) * mask

    fig = plt.figure(figsize=(10, 15))
    ax_raw = cb.image_plot(fig, fname[9:-4], 321, z_mesh, z_mesh)
    ax_zerfit = cb.image_plot(
        fig,
        "zernike(1~" + str(zer_fit_order) + ") fitting",
        322,
        fit * mask * 1e3,
        fit * mask * 1e3)

    ax_zer03 = cb.image_plot(
        fig, "zernike(<=" + str(3) + ")", 323, zer03mesh * mask, zer03mesh * mask)
    ax_zer10 = cb.image_plot(
        fig, "zernike(<=" + str(zer_fit_order) + ")", 324, zer10mesh * mask, zer10mesh * mask)

    ax_zer03diff = cb.image_plot(
        fig, "zernike(<=" + str(3) + ") removed", 325, diff_zer03, diff_zer03)
    ax_zer10diff = cb.image_plot(
        fig, "zernike(<=" + str(zer_fit_order) + ") removed", 326, diff_zer10, diff_zer10)
    fig.tight_layout()

    fig.savefig(mkfolder() + fname[9:-4] + ".png")

    np.savetxt(
        mkfolder() + "zer03_" + fname[9:-4] + ".csv",
        diff_zer03,
        delimiter=",")
    np.savetxt(
        mkfolder() + "zer10_" + fname[9:-4] + ".csv",
        diff_zer10,
        delimiter=",")

    np.savetxt(
        mkfolder() + fname[9:-4] + ".csv",
        z_mesh * mask,
        delimiter=",")
