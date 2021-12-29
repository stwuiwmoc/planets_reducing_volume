# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 20:12:37 2021

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
    
    exelis = pom.StitchedCsvToSurface(constants=CONSTS, 
                                      original_stitched_csv_fpath="mkfolder/exelis_rawdata_edit/exelis_reshaped.csv",
                                      deformed_stitched_csv_fpath="")
    
    removed = pom.ZernikeRemovedSurface(constants=CONSTS, 
                                        inputed_surface=exelis.surface, 
                                        removing_zernike_number_list=[4,5])

    removed.make_image_plot()
    
    target_surface = pom.ZernikeToSurface(constants = CONSTS, 
                                          zernike_number_list = [7],
                                          zernike_value_array = np.array([2e-7])) 
    
    filtered = pom.FilteredSurface(constants=CONSTS, 
                                   inputed_surface=exelis.surface, 
                                   filter_parameter=100)
    
    parent = pom.Surface(constants=CONSTS, 
                         surface=exelis.surface)