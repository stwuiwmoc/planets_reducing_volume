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
        ignore_radius=5e-3,
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

    target_zernike_value = 1e-6
    torque_value_limit = 5  # [mm] 駆動量の制約なしの場合はnp.inf

    for i in range(CONSTS.zernike_max_degree):
        print(i)
        zernike_num = i + 1
        target_zernike_value_array = np.zeros(CONSTS.zernike_max_degree)
        target_zernike_value_array[i] = target_zernike_value

        target_surface = pom.ZernikeToSurface(
            constants=CONSTS,
            zernike_value_array=target_zernike_value_array)

        wh_torque = pom.ZernikeToTorque(
            constants=CONSTS,
            target_zernike_value_array=target_zernike_value_array,
            ignore_zernike_number_list=[],
            restructed_torque_value=torque_value_limit)

        wh_reproducted_surface = pom.TorqueToSurface(
            constants=CONSTS,
            torque_value_array=wh_torque.torque_value_array
        )

        result_surface = pom.Surface(
            constants=CONSTS,
            surface=target_surface.surface - wh_reproducted_surface.surface)

        # plot

        cbar_reference_zernike_value_array = np.zeros(CONSTS.zernike_max_degree)
        cbar_reference_zernike_value_array[2] = target_zernike_value
        cbar_reference_surface = pom.ZernikeToSurface(
            constants=CONSTS,
            zernike_value_array=cbar_reference_zernike_value_array
        )

        fig1 = plt.figure(figsize=(10, 15))
        gs1 = fig1.add_gridspec(10, 2)
        fig1.suptitle("zernike " + str(zernike_num))

        ax15 = target_surface.make_image_plot(
            figure=fig1,
            position=gs1[0:3, 0],
            cbar_surface=cbar_reference_surface.surface,
        )
        ax15.set_title("target surface\n" + ax15.get_title())

        ax11 = wh_reproducted_surface.make_image_plot(
            figure=fig1,
            position=gs1[0:3, 1],
            cbar_surface=cbar_reference_surface.surface,
        )
        ax11.set_title("wh reproducted surface\n" + ax11.get_title())

        ax12 = result_surface.make_image_plot(
            figure=fig1,
            position=gs1[3:6, 0],
            cbar_surface=cbar_reference_surface.surface,
            pv_digits=5,
            rms_digits=5,
        )
        ax12.set_title("residual (colorbar is same for above)\n" + ax12.get_title())

        ax16 = result_surface.make_image_plot(
            figure=fig1,
            position=gs1[3:6, 1],
            pv_digits=5,
            rms_digits=5,
        )
        ax16.set_title("residual (colorbar is fitted to pv)\n" + ax16.get_title())

        ax13_xaxis = np.arange(CONSTS.zernike_max_degree) + 1
        ax13 = fig1.add_subplot(gs1[6:8, :])
        ax13.plot(
            ax13_xaxis,
            target_zernike_value_array,
            label="target zernike",
            marker="s",
        )
        ax13.plot(
            ax13_xaxis,
            np.dot(CONSTS.operation_matrix_A, wh_torque.torque_value_array),
            label="wh_reproducted",
            marker="s",
        )
        ax13.grid()
        ax13.legend()

        ax14 = wh_torque.make_torque_plot(
            figure=fig1,
            position=gs1[8:10, :],
        )

        if torque_value_limit is not np.inf:
            # set_ylim は np.inf には対応しないため
            ax14.set_ylim(-1.1 * torque_value_limit, 1.1 * torque_value_limit)

        fig1.tight_layout()

        picname1 = mkfolder() + "zer" + str(zernike_num).zfill(2) + "_fig1.png"
        fig1.savefig(picname1)
        fig1.clf()
