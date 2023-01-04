# %%

import time

import numpy as np
import pandas as pd


def mkfolder(suffix=""):
    """
    Parameters
    ----------
    suffix : str, optional
        The default is "".

    Returns
    -------
    str ( script name + suffix )
    """
    import os
    filename = os.path.basename(__file__)
    filename = filename.replace(".py", "") + suffix
    folder = "mkfolder/" + filename + "/"
    os.makedirs(folder, exist_ok=True)
    return folder


def read(filename):
    # skip行数の設定
    skip = 0
    with open(filename) as f:
        while True:
            line = f.readline()
            if line[0] == '%':
                skip += 1
            else:
                break
            # エラーの処理
            if len(line) == 0:
                break

    # データの読み出し
    df = pd.read_csv(
        filename,
        sep='\\s+',
        skiprows=skip,
        header=None)  # \s+...スペース数に関わらず区切る
    df.columns = ["x", "y", "z", "color", "dx", "dy", "dz"]
    df = df * 10**3  # [m] -> [mm]
    return df


if __name__ == "__main__":
    df0 = read("raw_data/Fxx/PM3.5_36ptAxWT06_F00.smesh.txt")

    fem_x_array = df0["x"].values * 1e-3  # [mm] -> [m] 換算
    fem_y_array = df0["y"].values * 1e-3  # [mm] -> [m] 換算
    h_0_array = df0["dz"].values * 1e-3  # [mm] -> [m] 換算
    fem_xy_matrix = np.stack([fem_x_array, fem_y_array], axis=1)

    np.savetxt(
        fname=mkfolder() + "fem_xy_matrix.csv",
        X=fem_xy_matrix,
        delimiter=","
    )

    file_num = 36
    data_length = len(df0)
    h_n_minus_h_0_matrix = np.zeros((data_length, file_num))

    for i in range(0, file_num):
        start = time.time()
        num = str(i + 1).zfill(2)

        data_fname = "raw_data/Fxx/PM3.5_36ptAxWT06_F" + num + ".smesh.txt"
        dfxx = read(data_fname)
        h_n_array = dfxx["dz"].values * 1e-3  # [mm] -> [m] 換算
        h_n_minus_h_0_matrix[:, i] = h_n_array - h_0_array

    np.savetxt(
        fname=mkfolder() + "h_n_minus_h_0_matrix.csv",
        X=h_n_minus_h_0_matrix,
        delimiter=",")
