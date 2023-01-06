# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 17:38:50 2021

@author: swimc
"""
# %%

import importlib

import matplotlib.pyplot as plt
import numpy as np

import planets_optimize_myclass as pom


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
        offset_height_percent=2,
        alpha_array=np.array([
            -74., +74., -40., -40., +76., -76.,
            -74., +74., -40., -40., +76., -76.,
            -74., +74., -40., -40., +76., -76.,
            -74., +74., -40., -40., +76., -76.,
            -74., +74., -40., -40., +76., -76.,
            -74., +74., -40., -40., +76., -76.,
        ])
    )

    torque_value_limit = 5

    cbar_min_percent_ = 60
    cbar_max_percent_ = 90
    zaxis_bottom_percent_ = 50

    target_surface = pom.KagiStitchToSurface(
        constants=CONSTS,
        txt_fpath="raw_data/0207xm130All2_0208cir_0208ykagomeAllCc.v4.8.hei_dense.txt")

    zernike_removed_surface = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=target_surface.surface,
        removing_zernike_number_list=[1, 2, 3, 4, 5, 6])

    reproducted_torque = pom.ZernikeToTorque(
        constants=CONSTS,
        target_zernike_value_array=zernike_removed_surface.zernike_value_array,
        ignore_zernike_number_list=[1, 2, 3, 4, 5, 6],
        restructed_torque_value=torque_value_limit)

    reproducted_surface = pom.TorqueToSurface(
        constants=CONSTS,
        torque_value_array=reproducted_torque.torque_value_array
    )

    reproducted_zernike_removed_surface = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=reproducted_surface.surface,
        removing_zernike_number_list=zernike_removed_surface.removing_zernike_number_list)

    result_surface = pom.Surface(
        constants=CONSTS,
        surface=zernike_removed_surface.surface - reproducted_zernike_removed_surface.surface)

    volume_reduciton_rate = 1 - result_surface.volume / zernike_removed_surface.volume

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

    ax14 = target_surface.make_zernike_value_plot(
        fig1, gs1[6:8, :],
        label_="measured_surface_zernike")

    ax14_xaxis = np.arange(CONSTS.zernike_max_degree) + 1
    ax14.plot(
        ax14_xaxis, zernike_removed_surface.zernike_value_array,
        marker="s", label="target_zernike (Z ≦ 6 is removed)")
    ax14.plot(
        ax14_xaxis, np.dot(CONSTS.operation_matrix_A, reproducted_torque.torque_value_array),
        marker="s", label="WH_reproducted_zernike")

    ax14.legend()
    ax14.set_ylim(
        target_surface.zernike_value_array.min() * 1.3,
        target_surface.zernike_value_array.max() * 1.3)

    ax13 = reproducted_torque.make_torque_plot(fig1, gs1[8:, :])
    ax13.set_title(ax13.get_title() + " Limit : " + str(torque_value_limit) + " [mm]")

    fig1.tight_layout()
    fig1.savefig(mkfolder() + "fig1.png")

    flat_surface = pom.Surface(
        constants=CONSTS,
        surface=np.ones([CONSTS.pixel_number, CONSTS.pixel_number]) * (zernike_removed_surface.pv * 0.3 + np.nanmin(zernike_removed_surface.surface)))

    fig2 = plt.figure(figsize=(16, 12))
    gs2 = fig2.add_gridspec(2, 2)

    ax21 = zernike_removed_surface.make_3d_plot(
        fig2, gs2[0, 1],
        None, cbar_min_percent_, cbar_max_percent_, zaxis_bottom_percent_)

    ax22 = flat_surface.make_3d_plot(
        fig2, gs2[0, 0],
        zernike_removed_surface.surface, cbar_min_percent_, cbar_max_percent_, zaxis_bottom_percent_)

    ax23 = result_surface.make_3d_plot(
        fig2, gs2[1, 1],
        zernike_removed_surface.surface, cbar_min_percent_, cbar_max_percent_, zaxis_bottom_percent_)

    ax24 = reproducted_zernike_removed_surface.make_3d_plot(
        fig2, gs2[1, 0],
        zernike_removed_surface.surface, cbar_min_percent_, cbar_max_percent_, zaxis_bottom_percent_)

    fig2.tight_layout()
