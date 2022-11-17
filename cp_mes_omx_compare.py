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

    mes_0_filepath = "raw_data/220117xrmEAi.v5.60.hei.txt"
    mes_n_filepath = "raw_data/220117xrm" + mes_n_serial + "Ei.v5.60.hei.txt"
