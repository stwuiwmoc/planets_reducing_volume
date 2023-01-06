# %%
import importlib

import matplotlib.pyplot as plt
import numpy as np

import planets_optimize_myclass as pom


def mkfolder(suffix=""):
    import os
    """
    Parameters
    ----------
    suffix : str, optional
        The default is "".

    Returns
    -------
    str ( script name + suffix )
    """
    filename = os.path.basename(__file__)
    filename = filename.replace(".py", "") + suffix
    folder = "mkfolder/" + filename + "/"
    os.makedirs(folder, exist_ok=True)
    return folder


if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=25e-3,
        pixel_number=1024,
        zernike_max_degree=11,
        offset_height_percent=2,
        alpha_array=np.array([
            -74., +74., -40., -40., +76., -76.,
            -74., +74., -40., -40., +76., -76.,
            -74., +74., -40., -40., +76., -76.,
            -74., +74., -40., -40., +76., -76.,
            -74., +74., -40., -40., +76., -76.,
            -74., +74., -40., -40., +76., -76.,
        ])
    )

    # %%
    importlib.reload(pom)
    for i in range(36):
        torque_value_array = np.zeros(36)
        torque_value_array[i] = +1

        surface = pom.TorqueToSurface(
            constants=CONSTS,
            torque_value_array=torque_value_array
        )

        zernike_removed_surface = pom.ZernikeRemovedSurface(
            constants=CONSTS,
            inputed_surface=surface.surface,
            removing_zernike_number_list=[1, 2, 3]
        )

        fig1 = plt.figure(figsize=(8, 6))
        gs1 = fig1.add_gridspec(3, 2)

        ax11 = surface.make_image_plot(
            figure=fig1,
            position=gs1[0:2, 0],
            cbar_min_percent=5,
            cbar_max_percent=95
        )
        ax11.set_title("M" + str(i + 1).zfill(2) + "\n" + ax11.get_title())

        ax12 = zernike_removed_surface.make_image_plot(
            figure=fig1,
            position=gs1[0:2, 1],
            cbar_min_percent=5,
            cbar_max_percent=95
        )
        ax12.set_title("zernike [1, 2, 3] removed\n" + ax12.get_title())

        ax13 = fig1.add_subplot(gs1[2, 0])
        ax13.plot(
            np.arange(CONSTS.zernike_max_degree) + 1,
            surface.zernike_value_array,
            marker="s"
        )
        ax13.grid()

        ax14 = fig1.add_subplot(gs1[2, 1])
        ax14.plot(
            np.arange(CONSTS.zernike_max_degree) + 1,
            zernike_removed_surface.zernike_value_array,
            marker="s"
        )
        ax14.grid()

        fig1.tight_layout()
        fig1.savefig(mkfolder() + "fig1_M" + str(i + 1).zfill(2) + ".png")
