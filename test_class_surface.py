# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 20:12:37 2021

@author: swimc
"""
# %%

import importlib

import numpy as np

import planets_optimize_myclass as pom

if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=25e-3,
        pixel_number=1024,
        zernike_max_degree=11,
        offset_height_percent=2)

    exelis = pom.ExelisCsvToSurface(
        constants=CONSTS)

    removed = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=exelis.surface,
        removing_zernike_number_list=[4, 5])

    oap = pom.OapSurface(
        constants=CONSTS,
        radius_of_curvature=8667e-3,
        off_axis_distance=1800e-3,
        clocking_angle_rad=0)

    target_surface = pom.ZernikeToSurface(
        constants=CONSTS,
        zernike_value_array=np.array([0, 0, 0, 0, 0, 0, 2e-7, 0, 0, 0, 0]))

    filtered = pom.FilteredSurface(
        constants=CONSTS,
        inputed_surface=exelis.surface,
        filter_parameter=100)

    fem = pom.FemTxtToSurface(
        constants=CONSTS,
        wt_version=6,
        wt_number=1)

    parent = pom.Surface(
        constants=CONSTS,
        surface=exelis.surface)
