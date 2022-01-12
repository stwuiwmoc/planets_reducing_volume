# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 17:01:34 2021

@author: swimc
"""
# %%

import numpy as np
import matplotlib.pyplot as plt

import planets_optimize_myclass as pom

def make_full_torque_value_array(torque_number_list, torque_value_aray):
    full_value_array = np.zeros(36)
    idx_array = np.array(torque_number_list) -1
    for i in range(len(idx_array)):
        idx = idx_array[i]
        full_value_array[idx] = torque_value_aray[i]
    return full_value_array


if __name__ == "__main__":
    CONSTS = pom.Constants(physical_radius=925e-3, 
                           ignore_radius=25e-3,
                           pixel_number=256,
                           zernike_max_degree = 10,
                           offset_height_percent=0)
    """
    # for parameter study
    original_torque_value_array = make_full_torque_value_array([1,7,13,19,25,31],
                                                               [5,-5,5,-5,5,-5])
    original_torque_value_array = make_full_torque_value_array([2,8,14,20,26,32],
                                                               [5,-5,5,-5,5,-5])
    original_torque_value_array = make_full_torque_value_array([4,10,16,22,28,34],
                                                               [5,-5,5,-5,5,-5])
    original_torque_value_array = make_full_torque_value_array([5,11,17,23,29,35],
                                                               [5,-5,5,-5,5,-5])    
    original_torque_value_array = make_full_torque_value_array([6,12,18,24,30,36],
                                                               [5,-5,5,-5,5,-5])
    
    original_torque_value_array = make_full_torque_value_array([1,7,13,19,25,31, 3,9,15,21,27,33, 4,10,16,22,28,34, 5,11,17,23,29,35, 6,12,18,24,30,36],
                                                               [-5,5,-5,5,-5,5, -5,5,-5,5,-5,5, -5,5,-5,5,-5,5, -5,5,-5,5,-5,5, 5,-5,5,-5,5,-5])
    """


    original_torque_value_array = make_full_torque_value_array([3,9,15,21,27,33],
                                                               [5,-5,5,-5,5,-5])
    wh_deformed_zernike = pom.TorqueToZernike(constants=CONSTS,
                                              torque_value_array=original_torque_value_array,
                                              restructed_torque_value=5,
                                              ignore_zernike_number_list=[1,2,3,4,5,6])
    
    wh_deformed_surface = pom.ZernikeToSurface(constants=CONSTS,
                                               zernike_number_list=wh_deformed_zernike.remaining_zernike_number_list,
                                               zernike_value_array=wh_deformed_zernike.remaining_reproducted_zernike_value_array)
        
    fig = plt.figure(figsize=(15,7))
    gs = fig.add_gridspec(2,3)
    wh_deformed_zernike.make_torque_plot(figure=fig, position=gs[0,2])
    wh_deformed_surface.make_image_plot(figure=fig, position=gs[0:2,0:2])
    wh_deformed_surface.make_circle_path_plot(figure=fig, position=gs[1,2])
    fig.tight_layout()

    