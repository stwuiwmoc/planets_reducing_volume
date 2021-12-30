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
                           zernike_max_degree=10)
    
    OAP = pom.OapConstants(ideal_radius_of_curvature=8667, 
                           ideal_off_axis_distance=1800, 
                           ideal_clocking_angle_rad=0, 
                           delta_radius_of_curvature=0, 
                           delta_off_axis_distance=0, 
                           delta_clocking_angle_rad=0)
  
    exelis = pom.StitchedCsvToSurface(constants=CONSTS, 
                                      original_stitched_csv_fpath="mkfolder/exelis_rawdata_edit/exelis_reshaped.csv",
                                      deformed_stitched_csv_fpath="")
    
