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
    
    exelis = pom.StitchedCsvToSurface(constants=CONSTS, 
                                      original_stitched_csv_fpath="mkfolder/exelis_rawdata_edit/exelis_reshaped.csv",
                                      deformed_stitched_csv_fpath="")
    
