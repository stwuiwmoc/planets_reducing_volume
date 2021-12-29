# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 12:19:46 2021

@author: swimc
"""
import numpy as np
import matplotlib.pyplot as plt

import planets_optimize_myclass as pom

def volume_check(radius, pv, ignore_lower_height_percent):
    a = -radius + 2*radius * ignore_lower_height_percent/100
    theta = np.arcsin(a/radius)
    
    square = pv/(3*radius) * (radius**2 - a**2)**(3/2)
    circle = pv/2 * a * radius * (np.pi/2 - theta - np.sin(2*theta)/2) 
    return square - circle

if __name__ == "__main__":
    CONSTS = pom.Constants(physical_radius=0.925, 
                           ignore_radius=0, 
                           pixel_number=1024, 
                           zernike_max_degree=10)
    
    zer2 = pom.ZernikeToSurface(constants=CONSTS, 
                                zernike_number_list=[2], 
                                zernike_value_array=[1.0020721132673347e-6],
                                offset_height_percent=0)
    fig1=plt.figure()
    zer2.make_image_plot(figure=fig1)
    print(zer2.volume)
    
    exelis = pom.StitchedCsvToSurface(constants=CONSTS, 
                                      original_stitched_csv_fpath="mkfolder/exelis_rawdata_edit/exelis_reshaped.csv", 
                                      None_or_deformed_stitched_csv_fpath=None,
                                      offset_height_percent=2)
    fig2=plt.figure()
    exelis.make_image_plot(fig2)