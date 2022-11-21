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


def make_full_torque_value_array_for_mes(mes_n_serial_: str) -> np.ndarray:
    """1/17 F～L の円環パス測定それぞれでのWH駆動量ベクトルを出力

    Parameters
    ----------
    mes_n_serial_ : str
        F, F, G, H, I, J, K, L のどれかの文字 \n

    Returns
    -------
    np.ndarray
        [mm] F～Lそれぞれの測定を行った時の駆動量ベクトル、len() = 36 \n
        Fの場合、以下のリンク先の1117fを測定した時のトルクベクトルを出力する \n
        https://www.dropbox.com/scl/fi/306tibn2etnpiwv9z5xcl/planets%E6%B8%AC%E5%AE%9A%E8%A8%98%E9%8C%B2.gsheet?dl=0&rlkey=oxj5n8wfy8jg02g3o8171ubbw#gid=641027351
    """

    if mes_n_serial_ == "F":
        original_torque_value_array_ = make_full_torque_value_array(
            [1, 7, 13, 19, 25, 31],
            [5, -5, 5, -5, 5, -5])

    elif mes_n_serial_ == "G":
        original_torque_value_array_ = make_full_torque_value_array(
            [2, 8, 14, 20, 26, 32],
            [5, -5, 5, -5, 5, -5])

    elif mes_n_serial_ == "H":
        original_torque_value_array_ = make_full_torque_value_array(
            [3, 9, 15, 21, 27, 33],
            [5, -5, 5, -5, 5, -5])

    elif mes_n_serial_ == "I":
        original_torque_value_array_ = make_full_torque_value_array(
            [4, 10, 15, 21, 28, 34],
            [5, -5, 5, -5, 5, -5])

    elif mes_n_serial_ == "J":
        original_torque_value_array_ = make_full_torque_value_array(
            [4, 10, 16, 22, 28, 34],
            [5, -5, 5, -5, 5, -5])

    elif mes_n_serial_ == "K":
        original_torque_value_array_ = make_full_torque_value_array(
            [5, 11, 17, 23, 29, 35],
            [5, -5, 5, -5, 5, -5])

    elif mes_n_serial_ == "L":
        original_torque_value_array_ = make_full_torque_value_array(
            [6, 12, 18, 24, 30, 36],
            [5, -5, 5, -5, 5, -5])

    else:
        print("Error! input must be any of F, F, G, H, I, J, K, L")

    return original_torque_value_array_


if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=25e-3,
        pixel_number=256,
        zernike_max_degree=11,
        offset_height_percent=0)

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
    full_torque_value_array = make_full_torque_value_array_for_mes(mes_n_serial_=mes_n_serial)
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

    # plot
    fig1 = plt.figure(figsize=(10, 10))
    gs1 = fig1.add_gridspec(4, 1)

    ax11 = fig1.add_subplot(gs1[0, 0])
    ax11.plot(
        mes0n.df_raw_original["degree"],
        mes0n.df_raw_original["height"],
        label="original"
    )
    ax11.plot(
        mes0n.df_raw_deformed["degree"],
        mes0n.df_raw_deformed["height"],
        label="deformed"
    )
    ax11.legend()
    ax11.grid()
    ax11.set_ylabel("height [m]")

    ax12 = fig1.add_subplot(gs1[1, 0])
    ax12.plot(
        mes0n_zerfit.degree_array,
        mes0n_zerfit.unprocessed_height_array,
        label="deformed - original",
        color="blue"
    )
    ax12.plot(
        mes0n_zerfit.degree_array,
        mes0n_zerfit.removing_zernike_height_array,
        label="zernike fit",
        color="blue",
        linestyle=":"
    )
    ax12.legend()
    ax12.grid()
    ax12.set_ylabel("height_diff [m]")

    ax13 = fig1.add_subplot(gs1[2, 0])
    ax13.plot(
        mes0n_zerfit.degree_array,
        mes0n_zerfit.zernike_removed_height_array,
        color="blue",
        linestyle="--"
    )
    ax13.grid()
    ax13.set_ylabel("zernike removed\nheight diff [m]")

    ax14 = fig1.add_subplot(gs1[3, 0])
    ax14.plot(
        mes0n_zerfit.r_const_zernike_number_meaning_list,
        mes0n_zerfit.r_const_zernike_polynomial_array,
        marker="s"
    )
    ax14.grid()
    ax14.set_xlabel("zernike number")
    ax14.set_ylabel("zernike polynomial value\nfor r-const [m]")

    fig1.tight_layout()

    fig2 = plt.figure(figsize=(12, 10))
    gs2 = fig2.add_gridspec(3, 4)

    ax21 = omx_deformed_zernike.make_torque_plot(
        figure=fig2,
        position=gs2[0, 0:2]
    )

    ax22 = omx_deformed_surface.make_image_plot(
        figure=fig2,
        position=gs2[1:3, 0:2]
    )

    ax23 = fig2.add_subplot(gs2[0, 2:4])
    ax23.set_title("operation matrix model")
    ax23.plot(
        omx_angle_array,
        omx_deformed_height_array,
        label="WH deforming",
        color="red"
    )
    ax23.plot(
        omx_angle_array,
        omx_deformed_height_array - omx_deformed_zernike_removed_height_array,
        label="zernike fit",
        color="red",
        linestyle=":"
    )
    ax23.legend()
    ax23.grid()
    ax23.set_ylabel("height_diff [m]")

    ax24 = fig2.add_subplot(gs2[1, 2:4])
    ax24.set_title("measurement")
    ax24.plot(
        mes0n_zerfit.degree_array,
        mes0n_zerfit.unprocessed_height_array,
        label="WH deforming\n(deformed - original)",
        color="blue",
    )
    ax24.plot(
        mes0n_zerfit.degree_array,
        mes0n_zerfit.removing_zernike_height_array,
        label="zernike fit",
        color="blue",
        linestyle=":"
    )
    ax24.legend()
    ax24.grid()
    ax24.set_ylabel("height_diff [m]")

    ax25 = fig2.add_subplot(gs2[2, 2:4])
    ax25.plot(
        mes0n_zerfit.degree_array,
        mes0n_zerfit.zernike_removed_height_array,
        color="blue",
        label="measurement",
        linestyle="--"
    )
    ax25.plot(
        omx_angle_array,
        omx_deformed_zernike_removed_height_array,
        label="operation matrix",
        color="red",
        linestyle="--"
    )
    ax25.grid()
    ax25.set_ylabel("zernike removed\nheight diff [m]")
    ax25.legend()

    fig2.tight_layout()
