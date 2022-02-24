# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 17:38:50 2021

@author: swimc
"""
# %%

import numpy as np
import matplotlib.pyplot as plt

import planets_optimize_myclass as pom
import importlib


def mkfolder(suffix=""):
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


if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=50e-3,
        pixel_number=256,
        zernike_max_degree=11,
        offset_height_percent=2)

    torque_value_limit = 5

    target_surface = pom.StitchedCsvToSurface(
        constants=CONSTS,
        original_stitched_csv_fpath="mkfolder/stitch2mesh/zer03_0207xm130All2_0208cir_0208ykagomeAllCc.v4.8.hei_dense.csv",
        deformed_stitched_csv_fpath="")

    zernike_removed_surface = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=target_surface.surface,
        removing_zernike_number_list=[1, 2, 3, 4, 5, 6])

    reproducted_torque = pom.ZernikeToTorque(
        constants=CONSTS,
        target_zernike_value_array=zernike_removed_surface.zernike_value_array,
        ignore_zernike_number_list=[1, 2, 3, 4, 5, 6],
        restructed_torque_value=torque_value_limit)

    reproducted_zernike = pom.TorqueToZernike(
        constants=CONSTS,
        torque_value_array=reproducted_torque.torque_value_array)

    reproducted_surface = pom.ZernikeToSurface(
        constants=CONSTS,
        zernike_value_array=reproducted_zernike.zernike_value_array)

    reproducted_zernike_removed_surface = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=reproducted_surface.surface,
        removing_zernike_number_list=zernike_removed_surface.removing_zernike_number_list)

    result_surface = pom.Surface(
        constants=CONSTS,
        surface=zernike_removed_surface.surface - reproducted_zernike_removed_surface.surface)

    volume_reduciton_rate = 1 - result_surface.volume / zernike_removed_surface.volume

    cbar_min_percent_ = 55
    cbar_max_percent_ = 95

    fig1 = plt.figure(figsize=(12, 16))
    gs1 = fig1.add_gridspec(10, 2)

    ax11 = zernike_removed_surface.make_image_plot(
        fig1, gs1[0:3, 0],
        None, cbar_min_percent_, cbar_max_percent_, 3, 3)
    ax11.set_title("measured (zer ≦ 6 is removed)\n" + ax11.get_title())

    ax17 = reproducted_surface.make_image_plot(fig1, gs1[0:3, 1])
    ax17.set_title("WH reproducted surface\n" + ax17.get_title())

    ax16 = reproducted_zernike_removed_surface.make_image_plot(
        fig1, gs1[3:6, 0],
        zernike_removed_surface.surface, cbar_min_percent_, cbar_max_percent_, 3, 3)
    ax16.set_title("WH reproducted surface (Z ≦ 6 is removed)\n" + ax16.get_title())

    ax12 = result_surface.make_image_plot(
        fig1, gs1[3:6, 1],
        zernike_removed_surface.surface, cbar_min_percent_, cbar_max_percent_, 3, 3)
    ax12.set_title(ax12.get_title() + "\nvolume_reduction = -" + str(round(volume_reduciton_rate * 1e2, 1)) + "%")

    ax14 = fig1.add_subplot(gs1[6:8, :])
    ax14_xaxis = np.arange(CONSTS.zernike_max_degree) + 1
    ax14.plot(
        ax14_xaxis, target_surface.zernike_value_array,
        marker="s", label="measured_surface_zernike")
    ax14.plot(
        ax14_xaxis, zernike_removed_surface.zernike_value_array,
        marker="s", label="target_zernike (Z ≦ 6 is removed)")
    ax14.plot(
        ax14_xaxis, reproducted_zernike.zernike_value_array,
        marker="s", label="WH_reproducted_zernike")

    ax14.legend()
    ax14.set_xticks(ax14_xaxis)
    ax14.set_xticklabels(ax14_xaxis)
    ax14.grid()
    ax14.set_ylim(
        target_surface.zernike_value_array.min() * 1.3,
        target_surface.zernike_value_array.max() * 1.3)

    ax13 = reproducted_zernike.make_torque_plot(fig1, gs1[8:, :])
    ax13.set_title(ax13.get_title() + " Limit : " + str(torque_value_limit) + " [mm]")

    fig1.tight_layout()
    fig1.savefig(mkfolder() + "fig1.png")
