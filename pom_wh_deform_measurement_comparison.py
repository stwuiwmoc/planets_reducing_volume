# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 17:29:18 2021

@author: swimc
"""

import planets_optimize_myclass as pom
import matplotlib.pyplot as plt

if __name__ == "__main__":
    CONSTS = pom.Constants(physical_radius=925e-3, 
                           ignore_radius=25e-3,
                           pixel_number=256,
                           zernike_max_degree=10)
    
    diff = pom.StitchedCsvToSurface(constants=CONSTS,
                                    original_stitched_csv_fpath="mkfolder/stitch2mesh/zer10_1215xm1301214ym870-510cir.v4.22.hei_dense.csv", 
                                    deformed_stitched_csv_fpath="mkfolder/stitch2mesh/zer10_1215xm1301216ym870-510cir.v4.21.hei_dense.csv")
    
    diff.make_image_plot()
    diff.make_circle_path_plot()