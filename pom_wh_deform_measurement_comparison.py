# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 17:29:18 2021

@author: swimc
"""
# %%
import planets_optimize_myclass as pom
import matplotlib.pyplot as plt

if __name__ == "__main__":
    CONSTS = pom.Constants(physical_radius=925e-3, 
                           ignore_radius=25e-3,
                           pixel_number=1024,
                           zernike_max_degree=10,
                           offset_height_percent=2)
    
    diff = pom.StitchedCsvToSurface(constants=CONSTS,
                                    original_stitched_csv_fpath="mkfolder/stitch2mesh/zer03_1215xm1301214ym870-510cir.v4.22.hei_dense.csv", 
                                    deformed_stitched_csv_fpath="mkfolder/stitch2mesh/zer03_1215xm1301216ym870-510cir.v4.21.hei_dense.csv")
    
    zernike_removed = pom.ZernikeRemovedSurface(constants=CONSTS, 
                                                inputed_surface=diff.surface, 
                                                removing_zernike_number_list=[1,2,3,4,5,7,8])
    
    fig1=plt.figure(figsize=(5,10))
    gs1=fig1.add_gridspec(2,1)
    diff.make_image_plot(figure=fig1, position=gs1[0,0])
    
    zernike_removed.make_image_plot(figure=fig1, position=gs1[1,0])
    fig1.tight_layout()
    
    fig2=plt.figure(figsize=(7,10))
    gs2=fig2.add_gridspec(3,1)
    zernike_removed.make_image_plot(figure=fig2, position=gs2[0:2, 0])
    zernike_removed.make_circle_path_plot(figure=fig2, position=gs2[2,0])
    fig2.tight_layout()
    
    plt.plot()
    