# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 17:07:48 2022

@author: swimc
"""

import numpy as np
import matplotlib.pyplot as plt

import planets_optimize_myclass as pom

if __name__ == "__main__":
    CONSTS = pom.Constants(physical_radius=925e-3, 
                           ignore_radius=25e-3,
                           pixel_number=1024,
                           zernike_max_degree=10,
                           offset_height_percent=2)
    
    OAP = pom.OapConstants(ideal_radius_of_curvature=8667e-3, 
                           ideal_off_axis_distance=1800e-3, 
                           ideal_clocking_angle_rad=0, 
                           delta_radius_of_curvature=0, 
                           delta_off_axis_distance=0, 
                           delta_clocking_angle_rad=0)
  
    exelis = pom.StitchedCsvToSurface(constants=CONSTS, 
                                      original_stitched_csv_fpath="mkfolder/exelis_rawdata_edit/exelis_reshaped.csv",
                                      deformed_stitched_csv_fpath="")
        
    ## smoothing
    filtered_exelis = pom.FilteredSurface(constants=CONSTS, 
                                          inputed_surface=exelis.surface, 
                                          filter_parameter=100)
    
    ideal_oap = pom.OapSurface(constants=CONSTS, 
                               radius_of_curvature=OAP.ideal_radius_of_curvature, 
                               off_axis_distance=OAP.ideal_off_axis_distance, 
                               clocking_angle_rad=OAP.ideal_clocking_angle_rad)
    
    ## oap minimize
    oap_minimize = pom.OapMinimize(constants=CONSTS, 
                                   oap_constants=OAP, 
                                   inputed_surface=filtered_exelis.surface)
    
    optimized_oap = pom.OapSurface(constants=CONSTS, 
                                   radius_of_curvature=oap_minimize.result_parameters[0], 
                                   off_axis_distance=oap_minimize.result_parameters[1], 
                                   clocking_angle_rad=oap_minimize.result_parameters[2])
    
    oap_optimized_exelis = pom.Surface(constants=CONSTS, 
                                       surface=filtered_exelis.surface+ideal_oap.surface-optimized_oap.surface)
    
    ## wh minimize
    zernike_removed = pom.ZernikeRemovedSurface(constants=CONSTS, 
                                                inputed_surface=oap_optimized_exelis.surface, 
                                                removing_zernike_number_list=[1,2,3,4,5,6])
    
    fitting_torque = pom.ZernikeToTorque(constants=CONSTS, 
                                         target_zernike_number_list=[i+1 for i in range(CONSTS.zernike_max_degree)], 
                                         target_zernike_value_array=zernike_removed.zernike_value_array, 
                                         ignore_zernike_number_list=zernike_removed.removing_zernike_number_list)
    
    wh_optimized_zernike = pom.TorqueToZernike(constants=CONSTS, 
                                               torque_value_array=fitting_torque.torque_value_array, 
                                               restructed_torque_value=5, 
                                               ignore_zernike_number_list=zernike_removed.removing_zernike_number_list)
    
    wh_optimized_surface = pom.ZernikeToSurface(constants=CONSTS, 
                                                zernike_number_list=wh_optimized_zernike.remaining_zernike_number_list, 
                                                zernike_value_array=wh_optimized_zernike.remaining_reproducted_zernike_value_array)
    
    wh_optimized_exelis = pom.Surface(constants=CONSTS, 
                                      surface=zernike_removed.inputed_surface-wh_optimized_surface.surface)