# %%
import importlib

import matplotlib.pyplot as plt
import numpy as np

import planets_optimize_myclass as pom

if __name__ == "__main__":
    importlib.reload(pom)

    circle_path_radius = 870e-3

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
    fig1 = plt.figure(figsize=(10, 10))
    gs1 = fig1.add_gridspec(4, 2)

    # とりあえず形状のブロット
    ax11 = surface.make_image_plot(
        figure=fig1, position=gs1[0:2, 0]
    )

    # 円環パスの高さ導出までやったものをプロット
    ax12 = fig1.add_subplot(gs1[2, :])
    ax12.plot(
        surface.calc_circle_path_height(radius=circle_path_radius, angle_division_number=360)[0],
        surface.calc_circle_path_height(radius=circle_path_radius, angle_division_number=360)[1],
    )
    ax12.grid()

    fig1.tight_layout()
