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
