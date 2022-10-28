# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 17:29:18 2021

@author: swimc
"""
import importlib

import matplotlib.pyplot as plt

# %%
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
        pixel_number=1024,
        zernike_max_degree=11,
        offset_height_percent=0)

    fname_natural = "raw_data/0228xm130_3degintAll_0301ym870-510CirAll_0228ykagomeACc.v4.6.hei_dense.txt"
    fname_active = "raw_data/0301xm130_3degintC_0302ym910-510cirAll.v4.6.hei_dense.txt"

    wt_natural = pom.KagiStitchToSurface(
        constants=CONSTS,
        txt_fpath=fname_natural)

    filter_natural = pom.FilteredSurface(
        constants=CONSTS,
        inputed_surface=wt_natural.surface,
        filter_parameter=0)

    zer03_natural = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=filter_natural.surface,
        removing_zernike_number_list=[1, 2, 3])

    zer06_natural = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=filter_natural.surface,
        removing_zernike_number_list=[1, 2, 3, 4, 5, 6])

    wt_active = pom.KagiStitchToSurface(
        constants=CONSTS,
        txt_fpath=fname_active)

    filter_active = pom.FilteredSurface(
        constants=CONSTS,
        inputed_surface=wt_active.surface,
        filter_parameter=0)

    zer03_active = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=filter_active.surface,
        removing_zernike_number_list=[1, 2, 3])

    zer06_active = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=filter_active.surface,
        removing_zernike_number_list=[1, 2, 3, 4, 5, 6])

    zer03_diff = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=filter_active.surface - filter_natural.surface,
        removing_zernike_number_list=[1, 2, 3])

    zer06_diff = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=filter_active.surface - filter_natural.surface,
        removing_zernike_number_list=[1, 2, 3, 4, 5, 6])

    fig1 = plt.figure(figsize=(12, 12))
    gs1 = fig1.add_gridspec(3, 3)
    ax11 = wt_natural.make_zernike_value_plot(fig1, gs1[0, 0])
    ax12 = zer06_natural.make_image_plot(
        fig1, gs1[0, 1], cbar_surface=zer03_active.surface)
    ax12.set_title("z≦6 removed\nA\n" + ax12.get_title())
    ax13 = zer03_natural.make_image_plot(
        fig1, gs1[0, 2], cbar_surface=zer03_active.surface)
    ax13.set_title("z≦3 removed\nA\n" + ax13.get_title())

    ax14 = wt_active.make_zernike_value_plot(fig1, gs1[1, 0])
    ax15 = zer06_active.make_image_plot(
        fig1, gs1[1, 1], cbar_surface=zer03_active.surface)
    ax15.set_title("B\n" + ax15.get_title())
    ax16 = zer03_active.make_image_plot(
        fig1, gs1[1, 2], cbar_surface=zer03_active.surface)
    ax16.set_title("B\n" + ax16.get_title())

    ax17 = zer06_diff.make_image_plot(
        fig1, gs1[2, 1], cbar_surface=zer03_active.surface)
    ax17.set_title("B - A\n" + ax17.get_title())
    ax18 = zer03_diff.make_image_plot(
        fig1, gs1[2, 2], cbar_surface=zer03_active.surface)
    ax18.set_title("B - A\n" + ax18.get_title())

    fig1.suptitle("A : " + fname_natural[9:] + "\nB : " + fname_active[9:])

    fig1.tight_layout()
    fig1.savefig(mkfolder() + "fig1.png")

    fig2 = plt.figure(figsize=(7, 10))
    gs2 = fig2.add_gridspec(3, 1)
    zer06_diff.make_image_plot(
        figure=fig2, position=gs2[0:2, 0], cbar_min_percent=20, cbar_max_percent=80)
    zer06_diff.make_circle_path_plot(figure=fig2, position=gs2[2, 0])
    fig2.tight_layout()
