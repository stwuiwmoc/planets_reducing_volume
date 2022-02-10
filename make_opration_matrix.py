# -*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import os
import proper as pr
import scipy as sp

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import mpl_toolkits.axes_grid1


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
    filename = os.path.basename(__file__)
    filename = filename.replace(".py", "") + suffix
    folder = "mkfolder/" + filename + "/"
    os.makedirs(folder, exist_ok=True)
    return folder


def read(filename):
    # skip行数の設定
    skip = 0
    with open(filename) as f:
        while True:
            line = f.readline()
            if line[0] == '%':
                skip += 1
            else:
                break
            # エラーの処理
            if len(line) == 0:
                break

    # データの読み出し
    df = pd.read_csv(
        filename,
        sep='\\s+',
        skiprows=skip,
        header=None)  # \s+...スペース数に関わらず区切る
    df.columns = ["x", "y", "z", "color", "dx", "dy", "dz"]
    df = df * 10**3  # [m] -> [mm]
    return df


def fem_interpolate(df_0, dfxx):
    # input [mm] -> output [mm]
    # mesh型のデータを格子点gridに補完
    x_old = dfxx["x"]
    y_old = dfxx["y"]
    dw_old = dfxx["dz"] - df_0["dz"]

    xy_old = np.stack([x_old, y_old], axis=1)
    dw_new = sp.interpolate.griddata(
        xy_old, dw_old, (xx, yy), method="linear", fill_value=0)
    return dw_new


def image_plot(fig, title, position, c, c_scale, cb_micron=True):
    cmap = cm.jet
    fs = 15

    ax = fig.add_subplot(position)
    ax.imshow(c, interpolation="nearest", cmap=cmap, vmin=np.nanmin(
        c_scale), vmax=np.nanmax(c_scale), origin="lower", extent=[-925, 925, -925, 925])
    ax.set_title(title, fontsize=fs)

    divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes('right', '5%', pad='3%')

    if cb_micron:
        norm = Normalize(
            vmin=np.nanmin(c_scale * 1000),
            vmax=np.nanmax(c_scale * 1000))  # c_scaleがmm表記だと見にくいのでμ表記に変更
        cbar_title = r"[$\mu$m]"
    else:
        norm = Normalize(
            vmin=np.nanmin(c_scale * 10**6),
            vmax=np.nanmax(c_scale * 10**6))
        cbar_title = "[nm]"
    mappable = cm.ScalarMappable(norm=norm, cmap=cm.jet)
    cbar = fig.colorbar(mappable, ax=ax, cax=cax)
    cbar.set_label(cbar_title, fontsize=fs)

    return ax


def df_plot(fig, title, position, df0, dfxx):
    fs = 15

    x = df0["x"]
    y = df0["y"]
    c = (dfxx["dz"] - df0["dz"]) * 10**6  # [mm] -> [nm]

    ax = fig.add_subplot(position)
    ax.scatter(x, y, c=c, cmap=cm.jet)
    ax.set_title(title, fontsize=fs)

    divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes('right', '5%', pad='3%')
    norm = Normalize(vmin=np.nanmin(c), vmax=np.nanmax(c))
    cbar_title = "[nm]"

    mappable = cm.ScalarMappable(norm=norm, cmap=cm.jet)
    cbar = fig.colorbar(mappable, ax=ax, cax=cax)
    cbar.set_label(cbar_title, fontsize=fs)
    return ax


def zer_term_plot(fig, title, position, zer):
    # input zer [m], plot zer[nm]
    ax = fig.add_subplot(position)
    fs = 15

    x = np.arange(1, zer_order + 1)
    ax.plot(x, zer * 10**9, marker="s", linewidth=1)
    ax.grid()
    ax.set_xlabel("zernike order", fontsize=fs)
    ax.set_xticks(x)
    ax.set_ylabel("zernike RMS [nm]", fontsize=fs)
    ax.set_title(title, fontsize=fs)
    return ax


if __name__ == '__main__':
    file_num = 36
    px = 512
    m1_radi = 1850 / 2
    zer_order = 11

    x_arr = y_arr = np.linspace(-m1_radi, m1_radi, px)
    xx, yy = np.meshgrid(x_arr, y_arr)

    mask = np.where(xx**2 + yy**2 <= m1_radi**2, 1, np.nan)
    tf = ~np.isnan(mask)

    zer_opration_matrix = np.empty((file_num, zer_order))
    opration_matrix = np.empty((file_num, tf.sum()))

    df0 = read("raw_data/Fxx/PM3.5_36ptAxWT06_F00.smesh.txt")
    act_tuning = np.array([
        37., 37., 20., 20., -38., 38.,
        37., 37., 20., 20., -38., 38.,
        37., 37., 20., 20., -38., 38.,
        37., 37., 20., 20., -38., 38.,
        37., 37., 20., 20., -38., 38.,
        37., 37., 20., 20., -38., 38.])

    """
    act_tuning = np.array([
        1, 1, 1, 1, -1, 1,
        1, 1, 1, 1, -1, 1,
        1, 1, 1, 1, -1, 1,
        1, 1, 1, 1, -1, 1,
        1, 1, 1, 1, -1, 1,
        1, 1, 1, 1, -1, 1])
    """

    for i in range(0, file_num):

        num = str(i + 1).zfill(2)

        data_fname = "raw_data/Fxx/PM3.5_36ptAxWT06_F" + num + ".smesh.txt"
        dfxx = read(data_fname)

        diff = act_tuning[i] * tf * fem_interpolate(df0, dfxx) / 1000  # [mm] -> [m]
        zer, fit = pr.prop_fit_zernikes(
            diff, tf, px / 2, zer_order, xc=px / 2, yc=px / 2, FIT=True)
        zer_opration_matrix[i] = zer

        diff = mask * diff * 1000  # [m] -> [mm]
        fit = mask * fit * 1000  # [m] -> [mm]

        fig1 = plt.figure(figsize=(11, 10))
        gs1 = fig1.add_gridspec(2, 2)
        ax_df = df_plot(fig1, "act" + num + " : 0.05 [Nm]", gs1[0, 0], df0, dfxx)
        ax_dz = image_plot(
            fig1,
            "act" + num + " : 1 [mm] ( x" + str(act_tuning[i]) + " )",
            gs1[0, 1], diff, diff, False)

        ax_zer = zer_term_plot(fig1, "zernike terms", gs1[1, 0], zer)
        ax_fit = image_plot(fig1, "zernike fitting", gs1[1, 1], fit, diff, False)

        fig1.tight_layout()

        picname = mkfolder() + "WT06_F" + num + ".png"
        fig1.savefig(picname)
        fig1.clf()

        print(num)

    zer_save_fname = mkfolder() + "WT06_zer" + str(zer_order) + "_opration_matrix[m].csv"
    np.savetxt(zer_save_fname, zer_opration_matrix, delimiter=",")
