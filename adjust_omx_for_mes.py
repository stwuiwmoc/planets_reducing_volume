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


def make_adjusted_operation_matrix(
        operation_matrix: np.ndarray,
        magnification_for_remainder_1: float,
        magnification_for_remainder_2: float,
        magnification_for_remainder_3: float,
        magnification_for_remainder_4: float,
        magnification_for_remainder_5: float,
        magnification_for_remainder_6: float) -> np.ndarray:
    """モデルから作った作用行列に倍率をかけて測定結果を再現する作用行列にする

    Parameters
    ----------
    operation_matrix : np.ndarray
        Constants.operation_matrix
    magnification_for_remainder_1 : float
        act [ 1,  7, 13, 19, 25, 31] にかける倍率
    magnification_for_remainder_2 : float
        act [ 2,  8, 14, 20, 26, 32] にかける倍率
    magnification_for_remainde_3 : float
        act [ 3,  9, 15, 21, 27, 33] にかける倍率
    magnification_for_remainder_4 : float
        act [ 4, 10, 16, 22, 28, 34] にかける倍率
    magnification_for_remainder_5 : float
        act [ 5, 11, 17, 23, 29, 35] にかける倍率
    magnification_for_remainder_6 : float
        act [ 6, 12, 18, 24, 30, 36] にかける倍率

    Returns
    -------
    np.ndarray
        測定結果を再現する作用行列
    """

    adjusted_operation_matrix = operation_matrix.copy()

    magnification_list = [
        magnification_for_remainder_1,
        magnification_for_remainder_2,
        magnification_for_remainder_3,
        magnification_for_remainder_4,
        magnification_for_remainder_5,
        magnification_for_remainder_6,
    ]

    for i in range(6):
        magnification = magnification_list[i]

        adjusted_operation_matrix[:, 0 + i] = magnification * adjusted_operation_matrix[:, 0 + i]
        adjusted_operation_matrix[:, 6 + i] = magnification * adjusted_operation_matrix[:, 6 + i]
        adjusted_operation_matrix[:, 12 + i] = magnification * adjusted_operation_matrix[:, 12 + i]
        adjusted_operation_matrix[:, 18 + i] = magnification * adjusted_operation_matrix[:, 18 + i]
        adjusted_operation_matrix[:, 24 + i] = magnification * adjusted_operation_matrix[:, 24 + i]
        adjusted_operation_matrix[:, 30 + i] = magnification * adjusted_operation_matrix[:, 30 + i]

    return adjusted_operation_matrix


if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=25e-3,
        pixel_number=256,
        zernike_max_degree=11,
        offset_height_percent=0)

    CONSTS.operation_matrix = make_adjusted_operation_matrix(
        operation_matrix=CONSTS.operation_matrix,
        magnification_for_remainder_1=0.75955,
        magnification_for_remainder_2=3.42924,
        magnification_for_remainder_3=0.29211,
        magnification_for_remainder_4=0.50750,
        magnification_for_remainder_5=2.64606,
        magnification_for_remainder_6=0.86704
    )

    mes_n_serial = "F"  # F, G, H, I, J, K, L から選択
    ignore_zernike_number_list = [1, 2, 3, 4, 5, 6, 7, 8, 11]

    # 測定結果の処理
    mes_0_filepath = "raw_data/220117xrmEAi.v5.60.hei.txt"
    mes_n_filepath = "raw_data/220117xrm" + mes_n_serial + "Ei.v5.60.hei.txt"

    mes0n = pom.CirclePathMeasurementTxtReading(
        Constants=CONSTS,
        original_txt_fpath=mes_0_filepath,
        deformed_txt_fpath=mes_n_filepath
    )

    mes0n_zerfit = pom.CirclePathZernikeFitting(
        Constants=CONSTS,
        circle_path_radius=mes0n.circle_path_radius,
        degree_array=mes0n.df_diff["degree"],
        unprocessed_height_array=mes0n.df_diff["height"],
        ignore_zernike_number_list=ignore_zernike_number_list
    )

    # 作用行列による計算
    full_torque_value_array = cp.make_full_torque_value_array_for_mes(mes_n_serial_=mes_n_serial)
    angle_division_number = 360

    omx_deformed_zernike = pom.TorqueToZernike(
        constants=CONSTS,
        torque_value_array=full_torque_value_array)

    omx_deformed_surface = pom.ZernikeToSurface(
        constants=CONSTS,
        zernike_value_array=omx_deformed_zernike.zernike_value_array)

    omx_deformed_zernike_removed_surface = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=omx_deformed_surface.surface,
        removing_zernike_number_list=ignore_zernike_number_list)

    omx_angle_array = omx_deformed_surface.calc_circle_path_height(
        radius=mes0n_zerfit.circle_path_radius,
        angle_division_number=angle_division_number
    )[0]

    omx_deformed_height_array = omx_deformed_surface.calc_circle_path_height(
        radius=mes0n_zerfit.circle_path_radius,
        angle_division_number=angle_division_number
    )[1]

    omx_deformed_zernike_removed_height_array = omx_deformed_zernike_removed_surface.calc_circle_path_height(
        radius=mes0n.circle_path_radius,
        angle_division_number=angle_division_number
    )[1]

    # 半径一定のzernike係数ベクトルに変換
    omx_r_const_zernike_value_array = cp.convert_to_r_const_zernike_value_array(
        zernike_value_array_=omx_deformed_surface.zernike_value_array
    )

    # 作用行列を現実に一致させるための倍率の導出
    zernike_09_magnification = mes0n_zerfit.r_const_zernike_polynomial_array[5] / omx_r_const_zernike_value_array[5]
    zernike_10_magnification = mes0n_zerfit.r_const_zernike_polynomial_array[6] / omx_r_const_zernike_value_array[6]
    total_magnification = \
        np.sqrt(mes0n_zerfit.r_const_zernike_polynomial_array[5]**2 + mes0n_zerfit.r_const_zernike_polynomial_array[6]**2) \
        / np.sqrt(omx_r_const_zernike_value_array[5]**2 + omx_r_const_zernike_value_array[6]**2)

    # plot
    # 測定と作用行列の比較プロット

    fig1 = plt.figure(figsize=(12, 10))
    gs1 = fig1.add_gridspec(3, 2)

    ax11 = omx_deformed_zernike.make_torque_plot(
        figure=fig1,
        position=gs1[0, 0]
    )

    ax12 = fig1.add_subplot(gs1[1, 0])
    ax12.plot(
        mes0n_zerfit.r_const_zernike_number_meaning_list,
        mes0n_zerfit.r_const_zernike_polynomial_array,
        marker="s",
        label="measurement",
        color="blue"
    )
    ax12.plot(
        mes0n_zerfit.r_const_zernike_number_meaning_list,
        omx_r_const_zernike_value_array,
        marker="s",
        label="operation matrix",
        color="red"
    )
    ax12.legend()
    ax12.grid()
    ax12.set_xlabel("zernike number")
    ax12.set_ylabel("zernike polynomial value\nfor r-const [m]")

    ax13 = fig1.add_subplot(gs1[0, 1])
    ax13.set_title("operation matrix model")
    ax13.plot(
        omx_angle_array,
        omx_deformed_height_array,
        label="WH deforming",
        color="red"
    )
    ax13.plot(
        omx_angle_array,
        omx_deformed_height_array - omx_deformed_zernike_removed_height_array,
        label="zernike fit",
        color="red",
        linestyle=":"
    )
    ax13.plot(
        omx_angle_array,
        omx_deformed_zernike_removed_height_array,
        label="deforming - fit",
        color="red",
        linestyle="--"
    )
    ax13.legend()
    ax13.grid()
    ax13.set_ylabel("height_diff [m]")

    ax14 = fig1.add_subplot(gs1[1, 1])
    ax14.set_title("measurement")
    ax14.plot(
        mes0n_zerfit.degree_array,
        mes0n_zerfit.unprocessed_height_array,
        label="WH deforming\n(deformed - original)",
        color="blue",
    )
    ax14.plot(
        mes0n_zerfit.degree_array,
        mes0n_zerfit.removing_zernike_height_array,
        label="zernike fit",
        color="blue",
        linestyle=":"
    )
    ax14.plot(
        mes0n_zerfit.degree_array,
        mes0n_zerfit.zernike_removed_height_array,
        label="deforming - fit",
        color="blue",
        linestyle="--"
    )
    ax14.legend()
    ax14.grid()
    ax14.set_ylabel("height_diff [m]")

    ax15 = fig1.add_subplot(gs1[2, 1])
    ax15.plot(
        mes0n_zerfit.degree_array,
        mes0n_zerfit.zernike_removed_height_array,
        color="blue",
        label="measurement",
        linestyle="--"
    )
    ax15.plot(
        omx_angle_array,
        omx_deformed_zernike_removed_height_array,
        label="operation matrix",
        color="red",
        linestyle="--"
    )
    ax15.grid()
    ax15.set_ylabel("zernike removed\nheight diff [m]")
    ax15.legend()

    parameter_list_1 = [
        pom.get_latest_commit_datetime(),
        ["Have some change", "from above commit", pom.have_some_change_in_git_status()],
        ["original", "", mes_0_filepath],
        ["deformed", "", mes_n_filepath],
        ["Z09 : mes / omx", round(zernike_09_magnification, 5), ""],
        ["Z10 : mes / omx", round(zernike_10_magnification, 5), ""],
        ["Z09 & Z10 : mes / omx", round(total_magnification, 5), ""]
    ]

    ax16 = pom.plot_parameter_table(
        fig=fig1,
        position=gs1[2, 0],
        parameter_table=parameter_list_1,
        fontsize=10
    )

    fig1.tight_layout()
