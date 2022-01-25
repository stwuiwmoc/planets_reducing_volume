# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 17:38:50 2021

@author: swimc
"""
# %%

import numpy as np

import planets_optimize_myclass as pom

if __name__ == "__main__":
    CONSTS = pom.Constants(physical_radius=925e-3,
                           ignore_radius=25e-3,
                           pixel_number=256,
                           zernike_max_degree=10,
                           offset_height_percent=2)

    target_surface = pom.ZernikeToSurface(constants=CONSTS,
                                          zernike_number_list=[7],
                                          zernike_value_array=np.array([2e-7]))

    reproducted_torque = pom.ZernikeToTorque(
        constants=CONSTS,
        target_zernike_number_list=target_surface.zernike_number_list,
        target_zernike_value_array=target_surface.zernike_value_array,
        ignore_zernike_number_list=[1])

    reproducted_zernike = pom.TorqueToZernike(
        constants=CONSTS,
        torque_value_array=reproducted_torque.torque_value_array,
        restructed_torque_value=5,
        ignore_zernike_number_list=reproducted_torque.ignore_zernike_number_list)

    reproducted_surface = pom.ZernikeToSurface(
        constants=CONSTS,
        zernike_number_list=reproducted_zernike.remaining_zernike_number_list,
        zernike_value_array=reproducted_zernike.remaining_reproducted_zernike_value_array)

    reproducted_restructed_surface = pom.ZernikeToSurface(
        constants=CONSTS,
        zernike_number_list=reproducted_zernike.remaining_zernike_number_list,
        zernike_value_array=reproducted_zernike.remaining_reproducted_restructed_zernike_value_array)
