# %%
import importlib

import matplotlib.pyplot as plt
import numpy as np

import planets_optimize_myclass as pom

if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=0,
        pixel_number=1024,
        zernike_max_degree=11,
        offset_height_percent=2,
        alpha_array=np.array([
            37., -37., 20., 20., -38., 38.,
            37., -37., 20., 20., -38., 38.,
            37., -37., 20., 20., -38., 38.,
            37., -37., 20., 20., -38., 38.,
            37., -37., 20., 20., -38., 38.,
            37., -37., 20., 20., -38., 38.
        ])
    )

    # %%
    importlib.reload(pom)
    for i in range(36):
        torque_value_array = np.zeros(36)
        torque_value_array[i] = 1

        surface = pom.TorqueToSurface(
            constants=CONSTS,
            torque_value_array=torque_value_array
        )

        fig1 = plt.figure()
        gs1 = fig1.add_gridspec(1, 1)

        ax11 = surface.make_image_plot(
            figure=fig1,
            position=gs1[0, 0],
            cbar_min_percent=5,
            cbar_max_percent=95
        )
