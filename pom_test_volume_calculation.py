# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 12:19:46 2021

@author: swimc
"""
import numpy as np

import planets_optimize_myclass as pom

def volume_check(radius, pv, ignore_lower_height_percent):
    a = -radius + 2*radius * ignore_lower_height_percent/100
    theta = np.arcsin(a/radius)
    
    circle = pv*(0.5*radius**2*(-theta - np.sin(2*theta)/2 + np.pi))
    square = pv/(3*radius) * (radius**2 - a**2)**(3/2)
    bottom = pv*ignore_lower_height_percent/100*2*(0.5*radius**2*(-theta - np.sin(2*theta)/2 + np.pi))
    return circle + square - bottom

if __name__ == "__main__":
    CONSTS = pom.Constants(physical_radius=0.925, 
                           ignore_radius=0.025, 
                           pixel_number=1024, 
                           zernike_max_degree=10)
    
    zer2per0 = pom.ZernikeToSurface(constants=CONSTS, 
                                    zernike_number_list=[2], 
                                    zernike_value_array=[1],
                                    ignore_lower_height_percent=0)
    
    zer2per50 = pom.ZernikeToSurface(constants=CONSTS, 
                                     zernike_number_list=[2], 
                                     zernike_value_array=[1],
                                     ignore_lower_height_percent=0)
    
    
    check = volume_check(CONSTS.physical_radius,
                         zer2per50.pv,
                         zer2per50.ignore_lower_height_percent)
    
    print("class ",zer2per50.volume)
    print("check ",check)

    
    exelis = pom.StitchedCsvToSurface(constants=CONSTS, 
                                      original_stitched_csv_fpath="mkfolder/exelis_rawdata_edit/exelis_reshaped.csv", 
                                      None_or_deformed_stitched_csv_fpath=None)
    