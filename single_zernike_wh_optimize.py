# %%

import numpy as np
import matplotlib.pyplot as plt

import planets_optimize_myclass as pom
import importlib


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
        ignore_radius=0,
        pixel_number=1024,
        zernike_max_degree=11,
        offset_height_percent=2)

    torque_value_limit = np.inf

    for i in range(CONSTS.zernike_max_degree):
        print(i)
        zernike_num = i + 1
        target_zernike_value_array = np.zeros(CONSTS.zernike_max_degree)
        target_zernike_value_array[i] = 5e-7

        target_surface = pom.ZernikeToSurface(
            constants=CONSTS,
            zernike_value_array=target_zernike_value_array)

        wh_torque = pom.ZernikeToTorque(
            constants=CONSTS,
            target_zernike_value_array=target_zernike_value_array,
            ignore_zernike_number_list=[],
            restructed_torque_value=torque_value_limit)

        wh_reproducted_zernike = pom.TorqueToZernike(
            constants=CONSTS,
            torque_value_array=wh_torque.torque_value_array)

        wh_reproducted_surface = pom.ZernikeToSurface(
            constants=CONSTS,
            zernike_value_array=wh_reproducted_zernike.zernike_value_array)

        result_surface = pom.Surface(
            constants=CONSTS,
            surface=target_surface.surface - wh_reproducted_surface.surface)

        fig1 = plt.figure(figsize=(10, 10))
        gs1 = fig1.add_gridspec(4, 2)
        fig1.suptitle("zernike " + str(zernike_num))

        ax11 = wh_reproducted_surface.make_image_plot(fig1, gs1[0:2, 0])

        ax12 = result_surface.make_image_plot(fig1, gs1[0:2, 1], pv_digits=5, rms_digits=5)

        ax13_xaxis = np.arange(CONSTS.zernike_max_degree) + 1
        ax13 = fig1.add_subplot(gs1[2, :])
        ax13.plot(
            ax13_xaxis, target_zernike_value_array,
            label="target zernike", marker="s")
        ax13.plot(
            ax13_xaxis, wh_reproducted_zernike.zernike_value_array,
            label="wh_reproducted", marker="s")
        ax13.grid()
        ax13.legend()

        ax14 = wh_reproducted_zernike.make_torque_plot(fig1, gs1[3, :])

        fig1.tight_layout()

        picname1 = mkfolder() + "zer" + str(zernike_num).zfill(2) + "_fig1.png"
        fig1.savefig(picname1)
        fig1.clf()
