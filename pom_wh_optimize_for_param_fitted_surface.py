# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 17:52:07 2021

@author: swimc
"""

import numpy as np
import matplotlib.pyplot as plt

import planets_optimize_myclass as pom

if __name__ == "__main__":
    CONSTS = pom.Constants(physical_radius=925e-3, 
                           ignore_radius=25e-3,
                           pixel_number=256,
                           zernike_max_degree = 10)
    
    
    ## combined_minimize.pyでopqのフィッティングが済んだ段階でのzernikeの値を代入
    param_minimized_zernike_value = np.array([ 8.11304455e-08, -1.01586202e-07, -3.38411547e-07,  3.02566783e-07, 2.10233957e-07, -2.01693302e-07, -6.40135092e-08,  1.15529214e-08, 3.01199936e-07, -1.78044987e-08])
    
    param_minimized_zernike_in_cb = pom.ZernikeToSurface(constants=CONSTS, 
                                                         zernike_number_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                                                         zernike_value_array = param_minimized_zernike_value)
    
    reproducted_torque = pom.ZernikeToTorque(constants=CONSTS,
                                             target_zernike_number_list = param_minimized_zernike_in_cb.zernike_number_list,
                                             target_zernike_value_array = param_minimized_zernike_in_cb.zernike_value_array,
                                             ignore_zernike_number_list=[1])
    
    reproducted_zernike = pom.TorqueToZernike(constants=CONSTS, 
                                              torque_value_array=reproducted_torque.torque_value_array, 
                                              restructed_torque_value=5, 
                                              ignore_zernike_number_list=reproducted_torque.ignore_zernike_number_list)
    
    
    reproducted_surface = pom.ZernikeToSurface(constants=CONSTS,
                                               zernike_number_list=reproducted_zernike.remaining_zernike_number_list,
                                               zernike_value_array=reproducted_zernike.remaining_reproducted_zernike_value_array)
    
    reproducted_restructed_surface = pom.ZernikeToSurface(constants=CONSTS,
                                                          zernike_number_list=reproducted_zernike.remaining_zernike_number_list,
                                                          zernike_value_array=reproducted_zernike.remaining_reproducted_restructed_zernike_value_array)
    
    fig = plt.figure(figsize=(15,10))
    gs = fig.add_gridspec(2,3)

    param_minimized_zernike_in_cb.make_image_plot(figure=fig,
                                                  position=gs[0,0])
    reproducted_zernike.make_torque_plot(figure=fig,
                                         position=gs[1,0:3])
    reproducted_surface.make_image_plot(figure=fig,
                                        position=gs[0,1])
    reproducted_restructed_surface.make_image_plot(figure=fig,
                                                   position=gs[0,2])
    fig.tight_layout()