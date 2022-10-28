# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 22:07:13 2021

@author: swimc
"""
# %%

import numpy as np


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
    raw = np.genfromtxt(
        "raw_data/digitFig01.csv",
        delimiter=",",
        encoding="utf-8_sig")

    reshaped = raw.reshape((1024, 1026))

    # 縦横の5.52182だけが並ぶ行と列を削除して1023*1023に成形
    deleted = np.delete(reshaped, 0, 0)
    deleted = np.delete(deleted, [0, 1024, 1025], 1)

    inner_nan_removed = np.where(
        np.isnan(deleted),
        np.nanmean(deleted),
        deleted)

    pixels_in_mm = 1e-3 * inner_nan_removed
    save_fname = mkfolder() + "exelis_reshaped_mm.csv"
    np.savetxt(save_fname, pixels_in_mm, delimiter=",", encoding="cp932", fmt="%.10f")
