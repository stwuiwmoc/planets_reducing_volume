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

    mes_raw0n = pom.CirclePathMeasurementCsvReading(
        Constants=CONSTS,
        original_csv_fpath=input_filename_0,
        deformed_csv_fpath=input_filename_n)

    mes_zer0n = pom.CirclePathZernikeFitting(
        Constants=CONSTS,
        circle_path_radius=870e-3,
        df_diff=mes_raw0n.df_diff,
        ignore_zernike_number_list=ignore_zernike_number_list)

    # operation matrix (omx) calculate
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

    # fem calculate
    fem_surface_array = np.zeros((CONSTS.pixel_number, CONSTS.pixel_number))
    for i in range(len(torque_number_list)):
        fem_single_surface = pom.FemTxtToSurface(
            constants=CONSTS,
            wt_version=6,
            wt_number=torque_number_list[i])

        fem_surface_array += fem_single_surface.surface * torque_value_array[i]
        del fem_single_surface

    fem_composited_surface = pom.Surface(
        constants=CONSTS,
        surface=fem_surface_array)

    fem_composited_zernike_removed_surface = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=fem_composited_surface.surface,
        removing_zernike_number_list=[1, 2, 3, 4])

    # omx and fem comparison
    fem_omx_diff_surface = pom.Surface(
        constants=CONSTS,
        surface=fem_composited_surface.surface - omx_deformed_surface.surface)

    fem_omx_diff_zernike_removed_surface = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=fem_omx_diff_surface.surface,
        removing_zernike_number_list=ignore_zernike_number_list)

    # plot

    fig1 = plt.figure(figsize=(7, 10))
    gs1 = fig1.add_gridspec(3, 2)

    ax11 = omx_deformed_zernike_removed_surface.make_image_plot(figure=fig1, position=gs1[0:2, :])
    ax11 = omx_deformed_zernike_removed_surface.make_torque_fulcrum_plot(ax=ax11, torque_value_array=original_torque_value_array)
    ax11.set_title("operation matrix model\n" + ax11.get_title())

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

    fig2 = plt.figure(figsize=(7, 10))
    gs2 = fig2.add_gridspec(3, 2)

    ax21 = fem_composited_zernike_removed_surface.make_image_plot(figure=fig2, position=gs2[0:2, :])
    ax21 = fem_composited_zernike_removed_surface.make_torque_fulcrum_plot(ax=ax21, torque_value_array=original_torque_value_array)
    ax21.set_title("FEM model\n" + ax21.get_title())

    ax22 = fem_composited_zernike_removed_surface.make_circle_path_plot(figure=fig2, position=gs2[2, :], line_label="model")
    ax22.plot(mes_zer0n.df_diff["degree"], mes_zer0n.height_removed * 1e9, label="measurement " + serials[ptn])
    ax22.legend()

    pv_fem = float(ax22.get_title()[19:25])
    pv_mes = 1e9 * (mes_zer0n.height_removed.max() - mes_zer0n.height_removed.min())

    ax22.set_title(
        "measurement pv = " + str(round(pv_mes, 2)) + " [nm] / " +
        "model pv = " + str(round(pv_fem, 2)) + " [nm]\n" +
        "measurement / model = " + str(round(pv_mes / pv_fem, 3)))

    fig2.tight_layout()
    fig2.savefig(mkfolder() + serials[ptn] + "_fig2.png")

    fig3 = plt.figure(figsize=(15, 10))
    gs3 = fig3.add_gridspec(2, 3)

    ax31 = fem_composited_surface.make_image_plot(
        figure=fig3, position=gs3[0, 0], cbar_surface=None)
    ax31.set_title("FEM model\n" + ax31.get_title())

    ax32 = fem_composited_zernike_removed_surface.make_image_plot(
        figure=fig3, position=gs3[1, 0], cbar_surface=fem_composited_surface.surface)

    ax33 = omx_deformed_surface.make_image_plot(
        figure=fig3, position=gs3[0, 1], cbar_surface=fem_composited_surface.surface)
    ax33.set_title("operation matrix model\n" + ax33.get_title())

    ax34 = omx_deformed_zernike_removed_surface.make_image_plot(
        figure=fig3, position=gs3[1, 1], cbar_surface=fem_composited_surface.surface)

    ax35 = fem_omx_diff_surface.make_image_plot(
        figure=fig3, position=gs3[0, 2], cbar_surface=None)
    ax35.set_title("residual by ( FEM - OM )\n" + ax35.get_title())

    ax36 = fem_omx_diff_zernike_removed_surface.make_image_plot(
        figure=fig3, position=gs3[1, 2], cbar_surface=fem_omx_diff_surface.surface)

    fig3.tight_layout()
    fig3.savefig(mkfolder() + serials[ptn] + "_fig3.png")
