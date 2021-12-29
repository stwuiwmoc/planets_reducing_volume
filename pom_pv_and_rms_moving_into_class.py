# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 11:24:33 2021

@author: swimc
"""

import planets_optimize_myclass as pom

if __name__ == "__main__":
    CONSTS = pom.Constants(physical_radius=0.925, 
                           ignore_radius=0.25, 
                           pixel_number=256, 
                           zernike_max_degree=10)
    
    parent = pom.ZernikeToSurface(constants=CONSTS, 
                                  zernike_number_list=[7], 
                                  zernike_value_array=[1e-6])
    print(parent.pv)
    
    child = pom.StitchedCsvToSurface(constants=CONSTS, 
                                     original_stitched_csv_fpath="mkfolder/exelis_rawdata_edit/exelis_reshaped.csv", 
                                     None_or_deformed_stitched_csv_fpath=None)
    
    print(child.pv)