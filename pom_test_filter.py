# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 19:15:52 2021

@author: swimc
"""

import numpy as np
import matplotlib.pyplot as plt

import planets_optimize_myclass as pom

if __name__ == "__main__":
    CONSTS = pom.Constants(physical_radius=925e-3, 
                           ignore_radius=75e-3,
                           pixel_number=1024,
                           zernike_max_degree=10)
    
    exelis = pom.StitchedCsvToSurface(constants=CONSTS, 
                                      original_stitched_csv_fpath="mkfolder/exelis_rawdata_edit/exelis_reshaped.csv",
                                      None_or_deformed_stitched_csv_fpath=None)
    
    filtered = pom.FilteredSurface(constants=CONSTS, 
                                   inputed_surface=exelis.surface, 
                                   filter_parameter=100)
    
    filtered.make_image_plot()