# %%

import matplotlib.pyplot as plt
import importlib
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


def make_full_torque_value_array(torque_number_list_, torque_value_array_):
    full_value_array = np.zeros(36)
    idx_array = np.array(torque_number_list_) - 1
    for i in range(len(idx_array)):
        idx = idx_array[i]
        full_value_array[idx] = torque_value_array_[i]
    return full_value_array


if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=25e-3,
        pixel_number=256,
        zernike_max_degree=11,
        offset_height_percent=0)

    ignore_zernike_number_list = [1, 2, 3, 4, 5, 6, 7, 8]

    # measurement
    ptn: int = 6

    # serials = ["0117e", "0117f", "0117g", "0117h", "0117j", "0117k", "0117l", "0117i"]
    # folderpath = "mkfolder/psm_test_kagi_data_integration/"

    serials = ["0217ym870CirB", "", "", "", "0217ym870CirD", "0217ym870CirE", "0217ym870CirF"]
    folderpath = "mkfolder/kagi_hei_to_csv/"

    input_filename_0 = "../sag_integration_code/" + folderpath + serials[0] + "_height.csv"
    input_filename_n = "../sag_integration_code/" + folderpath + serials[ptn] + "_height.csv"

    mes_raw0n = pom.CirclePathMeasurementReading(
        Constants=CONSTS,
        original_csv_fpath=input_filename_0,
        deformed_csv_fpath=input_filename_n)

    mes_zer0n = pom.CirclePathZernikeFitting(
        Constants=CONSTS,
        circle_path_radius=870e-3,
        df_diff=mes_raw0n.df_diff,
        ignore_zernike_number_list=ignore_zernike_number_list)

    # calculate
    torque_number_list = [0 + ptn, 6 + ptn, 12 + ptn, 18 + ptn, 24 + ptn, 30 + ptn]
    torque_value_array = np.array([5, -5, 5, -5, 5, -5])

    original_torque_value_array = make_full_torque_value_array(
        torque_number_list_=torque_number_list,
        torque_value_array_=torque_value_array)

    omx_deformed_zernike = pom.TorqueToZernike(
        constants=CONSTS,
        torque_value_array=original_torque_value_array)

    omx_deformed_surface = pom.ZernikeToSurface(
        constants=CONSTS,
        zernike_value_array=omx_deformed_zernike.zernike_value_array)

    omx_deformed_zernike_removed_surface = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=omx_deformed_surface.surface,
        removing_zernike_number_list=ignore_zernike_number_list)

    # plot

    fig1 = plt.figure(figsize=(7, 10))
    gs1 = fig1.add_gridspec(3, 2)

    ax11 = omx_deformed_zernike_removed_surface.make_image_plot(figure=fig1, position=gs1[0:2, :])
    ax11 = omx_deformed_zernike_removed_surface.make_torque_fulcrum_plot(ax=ax11, torque_value_array=original_torque_value_array)

    ax12 = omx_deformed_zernike_removed_surface.make_circle_path_plot(figure=fig1, position=gs1[2, :], line_label="model")
    ax12.plot(mes_zer0n.df_diff["degree"], mes_zer0n.height_removed * 1e9, label="measurement " + serials[ptn])
    ax12.legend()

    pv_omx = float(ax12.get_title()[19:25])
    pv_mes = 1e9 * (mes_zer0n.height_removed.max() - mes_zer0n.height_removed.min())

    ax12.set_title(
        "measurement pv = " + str(round(pv_mes, 2)) + " [nm] / " +
        "model pv = " + str(round(pv_omx, 2)) + " [nm]\n" +
        "measurement / model = " + str(round(pv_mes / pv_omx, 3)))

    fig1.tight_layout()
    fig1.savefig(mkfolder() + serials[ptn] + "_fig1.png")
