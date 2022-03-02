# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 17:07:48 2022

@author: swimc
"""
# %%

import matplotlib.pyplot as plt
import importlib

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
        ignore_radius=25e-3,
        pixel_number=1024,
        zernike_max_degree=10,
        offset_height_percent=2)

    OAP = pom.OapConstants(
        ideal_radius_of_curvature=8667e-3,
        ideal_off_axis_distance=1800e-3,
        ideal_clocking_angle_rad=0,
        delta_radius_of_curvature=0,
        delta_off_axis_distance=0,
        delta_clocking_angle_rad=0)

    exelis = pom.ExelisCsvToSurface(
        constants=CONSTS)

    # smoothing
    filtered_exelis = pom.FilteredSurface(
        constants=CONSTS,
        inputed_surface=exelis.surface,
        filter_parameter=100)

    ideal_oap = pom.OapSurface(
        constants=CONSTS,
        radius_of_curvature=OAP.ideal_radius_of_curvature,
        off_axis_distance=OAP.ideal_off_axis_distance,
        clocking_angle_rad=OAP.ideal_clocking_angle_rad)

    # oap minimize
    oap_minimize = pom.OapMinimize(
        constants=CONSTS,
        oap_constants=OAP,
        inputed_surface=filtered_exelis.surface)

    optimized_oap = pom.OapSurface(
        constants=CONSTS,
        radius_of_curvature=oap_minimize.result_parameters[0],
        off_axis_distance=oap_minimize.result_parameters[1],
        clocking_angle_rad=oap_minimize.result_parameters[2])

    oap_optimized_exelis = pom.Surface(
        constants=CONSTS,
        surface=filtered_exelis.surface +
        ideal_oap.surface -
        optimized_oap.surface)

    # wh minimize
    zernike_removed = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=filtered_exelis.surface,
        removing_zernike_number_list=[1, 2, 3, 4, 5, 6])

    fitting_torque = pom.ZernikeToTorque(
        constants=CONSTS,
        target_zernike_value_array=zernike_removed.zernike_value_array,
        ignore_zernike_number_list=[1, 2, 3, 4, 5, 6],
        restructed_torque_value=5)

    wh_optimized_zernike = pom.TorqueToZernike(
        constants=CONSTS,
        torque_value_array=fitting_torque.torque_value_array)

    wh_optimized_surface = pom.ZernikeToSurface(
        constants=CONSTS,
        zernike_value_array=wh_optimized_zernike.zernike_value_array)

    wh_optimized_exelis = pom.Surface(
        constants=CONSTS,
        surface=filtered_exelis.surface - wh_optimized_surface.surface)

    fig = plt.figure(figsize=(10, 22))
    gs = fig.add_gridspec(4, 2)

    exelis.make_image_plot(figure=fig, position=gs[0, 0])
    filtered_exelis.make_image_plot(figure=fig, position=gs[0, 1])
    oap_optimized_exelis.make_image_plot(figure=fig, position=gs[1, 1])
    wh_optimized_surface.make_image_plot(figure=fig, position=gs[2, 0])
    wh_optimized_exelis.make_image_plot(figure=fig, position=gs[2, 1])
    wh_optimized_zernike.make_torque_plot(figure=fig, position=gs[3, :])

    fig.savefig(mkfolder() + "test.png")
