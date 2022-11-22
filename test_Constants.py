# %%

import importlib

import matplotlib.pyplot as plt
import numpy as np

import planets_optimize_myclass as pom

if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=25e-3,
        pixel_number=256,
        zernike_max_degree=11,
        offset_height_percent=0)

    surface = pom.FemTxtToSurface(
        constants=CONSTS,
        wt_version=6,
        wt_number=1
    )

    fig1 = plt.figure()
    surface.make_image_plot(fig1, 111)
