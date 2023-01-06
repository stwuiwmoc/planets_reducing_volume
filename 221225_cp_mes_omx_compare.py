# %%

import importlib

import matplotlib.pyplot as plt
import numpy as np

import cp_mes_omx_compare as cp
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

    ignore_zernike_number_list = [1, 2, 3, 4, 7, 8, 11]
    angle_division_number = 360

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=25e-3,
        pixel_number=256,
        zernike_max_degree=11,
        offset_height_percent=0,
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
    # 測定結果の処理
    importlib.reload(pom)
    mesde_filepath = "raw_data/221225xcmM_R850e-R850d.hei.txt"

    mesde = pom.CirclePathMeasurementTxtReading(
        Constants=CONSTS,
        original_txt_fpath=mesde_filepath,
        deformed_txt_fpath=None
    )

    mesde_zerfit = pom.CirclePathZernikeFitting(
        Constants=CONSTS,
        circle_path_radius=mesde.circle_path_radius,
        degree_array=mesde.df_diff["degree"],
        unprocessed_height_array=mesde.df_diff["height"],
        ignore_zernike_number_list=ignore_zernike_number_list
    )

    # 作用行列の処理
    full_torque_value_array = cp.make_full_torque_value_array(
        torque_number_list_=[
            3, 9, 15, 21, 27, 33,
            4, 10, 16, 22, 28, 34],
        torque_value_array_=[
            +5, -5, +5, -5, +5, -5,
            +5, -5, +5, -5, +5, -5]
    )

    omx_deformed_surface = pom.TorqueToSurface(
        constants=CONSTS,
        torque_value_array=full_torque_value_array
    )

    omx_deformed_zernike_removed_surface = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=omx_deformed_surface.surface,
        removing_zernike_number_list=ignore_zernike_number_list
    )

    omx_angle_array = omx_deformed_surface.calc_circle_path_height(
        radius=mesde_zerfit.circle_path_radius,
        angle_division_number=angle_division_number
    )[0]

    omx_deformed_height_array = omx_deformed_surface.calc_circle_path_height(
        radius=mesde_zerfit.circle_path_radius,
        angle_division_number=angle_division_number
    )[1]

    omx_deformed_zernike_removed_height_array = omx_deformed_zernike_removed_surface.calc_circle_path_height(
        radius=mesde_zerfit.circle_path_radius,
        angle_division_number=angle_division_number
    )[1]

    # 半径一定のzernike係数ベクトルに変換
    omx_r_const_zernike_value_array = cp.convert_to_r_const_zernike_value_array(
        zernike_value_array_=omx_deformed_surface.zernike_value_array
    )

    # 測定結果のplot
    fig1 = plt.figure(figsize=(10, 12))
    gs1 = fig1.add_gridspec(3, 1)

    ax11 = fig1.add_subplot(gs1[0, 0])
    ax11.plot(
        mesde_zerfit.degree_array,
        mesde_zerfit.unprocessed_height_array,
        label="c - d",
        color="blue"
    )
    ax11.plot(
        mesde_zerfit.degree_array,
        mesde_zerfit.removing_zernike_height_array,
        label="zernike fit",
        color="blue",
        linestyle=":"
    )
    ax11.grid()
    ax11.legend()
    ax11.set_ylabel("height diff [m]")

    ax12 = fig1.add_subplot(gs1[1, 0])
    ax12.plot(
        mesde_zerfit.degree_array,
        mesde_zerfit.zernike_removed_height_array,
        color="blue",
        linestyle="--"
    )
    ax12.grid()
    ax12.set_ylabel("zernike removed\nheight diff [m]")

    ax13 = fig1.add_subplot(gs1[2, 0])
    ax13.plot(
        mesde_zerfit.r_const_zernike_number_meaning_list,
        mesde_zerfit.r_const_zernike_polynomial_array,
        marker="s"
    )
    ax13.grid()
    ax13.set_xlabel("zernike number")
    ax13.set_ylabel("zernike polynomial value\nfor r-const [m]")

    fig1.tight_layout()

    # 作用行列結果のplot
    fig2 = plt.figure(figsize=(12, 10))
    gs2 = fig2.add_gridspec(3, 4)

    ax21 = omx_deformed_surface.make_image_plot(
        figure=fig2,
        position=gs2[1:3, 0:2]
    )

    ax22 = fig2.add_subplot(gs2[0, 2:4])
    ax22.plot(
        omx_angle_array,
        omx_deformed_height_array,
        label="WH deforming",
        color="red"
    )

    ax22.plot(
        omx_angle_array,
        omx_deformed_height_array - omx_deformed_zernike_removed_height_array,
        label="zernike fit",
        color="red",
        linestyle=":"
    )
    ax22.legend()
    ax22.grid()
    ax22.set_ylabel("height diff [m]")

    ax23 = fig2.add_subplot(gs2[1, 2:4])
    ax23.plot(
        omx_angle_array,
        omx_deformed_zernike_removed_height_array,
        color="red",
        linestyle="--"
    )
    ax23.grid()
    ax23.set_ylabel("zernike removed\nheight diff [m]")

    ax24 = omx_deformed_surface.make_zernike_value_plot(
        figure=fig2,
        position=gs2[2, 2:4]
    )

    fig2.tight_layout()

    # 測定と作用行列結果の比較プロット

    fig3 = plt.figure(figsize=(12, 10))
    gs3 = fig3.add_gridspec(3, 2)

    ax31 = fig3.add_subplot(gs3[0, 0])

    ax31.plot(
        mesde_zerfit.r_const_zernike_number_meaning_list,
        mesde_zerfit.r_const_zernike_polynomial_array,
        marker="s",
        label="measurement",
        color="blue"
    )
    ax31.plot(
        mesde_zerfit.r_const_zernike_number_meaning_list,
        omx_r_const_zernike_value_array,
        marker="s",
        label="operation matrix",
        color="red"
    )
    ax31.legend()
    ax31.grid()
    ax31.set_xlabel("zernike number")
    ax31.set_ylabel("zernike polynomial value\nfor r-const [m]")

    ax32 = fig3.add_subplot(gs3[1, 0])
    ax32.set_title("operation matrix")
    ax32.plot(
        omx_angle_array,
        omx_deformed_height_array,
        label="WH deforming",
        color="red"
    )
    ax32.plot(
        omx_angle_array,
        omx_deformed_height_array - omx_deformed_zernike_removed_height_array,
        label="zernike fit",
        color="red",
        linestyle=":"
    )
    ax32.plot(
        omx_angle_array,
        omx_deformed_zernike_removed_height_array,
        label="deforming - fit",
        color="red",
        linestyle="--"
    )
    ax32.legend()
    ax32.grid()
    ax32.set_ylabel("height_diff [m]")

    ax33 = fig3.add_subplot(gs3[1, 1])
    ax33.set_title("measurement")
    ax33.plot(
        mesde_zerfit.degree_array,
        mesde_zerfit.unprocessed_height_array,
        label="WH deforming\n(deformed - original)",
        color="blue",
    )
    ax33.plot(
        mesde_zerfit.degree_array,
        mesde_zerfit.removing_zernike_height_array,
        label="zernike fit",
        color="blue",
        linestyle=":"
    )
    ax33.plot(
        mesde_zerfit.degree_array,
        mesde_zerfit.zernike_removed_height_array,
        label="deforming - fit",
        color="blue",
        linestyle="--"
    )
    ax33.legend()
    ax33.grid()
    ax33.set_ylabel("height_diff [m]")

    ax35 = fig3.add_subplot(gs3[2, 0])
    ax35.set_title("compare")
    ax35.plot(
        mesde_zerfit.degree_array,
        mesde_zerfit.unprocessed_height_array,
        color="blue",
        label="measurement"
    )
    ax35.plot(
        omx_angle_array,
        omx_deformed_height_array,
        label="operation matrix",
        color="red"
    )
    ax35.grid()
    ax35.set_ylabel("height diff [m]")
    ax35.legend()

    ax34 = fig3.add_subplot(gs3[2, 1])
    ax34.set_title("compare (zernike removed)")
    ax34.plot(
        mesde_zerfit.degree_array,
        mesde_zerfit.zernike_removed_height_array,
        color="blue",
        label="measurement",
        linestyle="--"
    )
    ax34.plot(
        omx_angle_array,
        omx_deformed_zernike_removed_height_array,
        label="operation matrix",
        color="red",
        linestyle="--"
    )
    ax34.grid()
    ax34.set_ylabel("zernike removed\nheight diff [m]")
    ax34.legend()

    fig3.tight_layout()
