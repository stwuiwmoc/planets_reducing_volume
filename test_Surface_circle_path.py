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
        pixel_number=1024,
        zernike_max_degree=11,
        offset_height_percent=2)

    input_zernike_value_array = np.array(
        [0, 1e-6, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    )

    surface = pom.ZernikeToSurface(
        constants=CONSTS,
        zernike_value_array=input_zernike_value_array)

    # plot
    fig1 = plt.figure(figsize=(10, 15))
    gs1 = fig1.add_gridspec(3, 2)

    ax11 = surface.make_image_plot(
        figure=fig1, position=gs1[0, 0]
    )

    fig1.tight_layout()
