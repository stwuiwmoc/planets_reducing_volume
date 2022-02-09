# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 17:38:50 2021

@author: swimc
"""
# %%

import numpy as np
import matplotlib.pyplot as plt

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
    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=50e-3,
        pixel_number=256,
        zernike_max_degree=10,
        offset_height_percent=2)

    target_surface = pom.StitchedCsvToSurface(
        constants=CONSTS,
        original_stitched_csv_fpath="mkfolder/stitch2mesh/zer03_0131xm130allcc_0201cirAll.v4.8.hei_dense.csv",
        deformed_stitched_csv_fpath="")

    zernike_removed_surface = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=target_surface.surface,
        removing_zernike_number_list=[1, 2, 3, 4, 5, 6])

    reproducted_torque = pom.ZernikeToTorque(
        constants=CONSTS,
        target_zernike_number_list=np.arange(1, 11),
        target_zernike_value_array=zernike_removed_surface.zernike_value_array,
        ignore_zernike_number_list=[1])

    reproducted_zernike = pom.TorqueToZernike(
        constants=CONSTS,
        torque_value_array=reproducted_torque.torque_value_array,
        ignore_zernike_number_list=reproducted_torque.ignore_zernike_number_list)

    reproducted_surface = pom.ZernikeToSurface(
        constants=CONSTS,
        zernike_number_list=reproducted_zernike.remaining_zernike_number_list,
        zernike_value_array=reproducted_zernike.remaining_reproducted_zernike_value_array)

    result_surface = pom.Surface(
        constants=CONSTS,
        surface=target_surface.surface - reproducted_surface.surface)

    cbar_min_percent_ = 15
    cbar_max_percent_ = 90

    fig1 = plt.figure(figsize=(12, 12))
    gs1 = fig1.add_gridspec(3, 2)

    ax11 = zernike_removed_surface.make_image_plot(
        fig1, gs1[0:2, 0],
        None, cbar_min_percent_, cbar_max_percent_)

    ax12 = result_surface.make_image_plot(
        fig1, gs1[0:2, 1],
        zernike_removed_surface.surface, cbar_min_percent_, cbar_max_percent_)

    ax13 = reproducted_zernike.make_torque_plot(fig1, gs1[2, :])
    fig1.tight_layout()
    fig1.savefig(mkfolder() + "fig1.png")
