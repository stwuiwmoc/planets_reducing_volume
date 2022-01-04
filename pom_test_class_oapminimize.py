# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 21:32:06 2021

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
                                      
    
    filtered_exelis = pom.FilteredSurface(constants=CONSTS, 
                                          inputed_surface=exelis.surface, 
                                          filter_parameter=100)
                                          
    
    ideal_oap = pom.OapSurface(constants=CONSTS, 
                               radius_of_curvature=OAP.ideal_radius_of_curvature, 
                               off_axis_distance=OAP.ideal_off_axis_distance, 
                               clocking_angle_rad=OAP.ideal_clocking_angle_rad)
                               
    
    oap_minimize = pom.OapMinimize(constants=CONSTS, 
                                   oap_constants=OAP, 
                                   inputed_surface=filtered_exelis.surface)
                                   
    
    optimized_oap = pom.OapSurface(constants=CONSTS, 
                                   radius_of_curvature=oap_minimize.result_parameters[0], 
                                   off_axis_distance=oap_minimize.result_parameters[1], 
                                   clocking_angle_rad=oap_minimize.result_parameters[2])
                                   
    
    oap_optimized_exelis = pom.Surface(constants=CONSTS, 
                                       surface=filtered_exelis.surface+ideal_oap.surface-optimized_oap.surface)
                                       
    