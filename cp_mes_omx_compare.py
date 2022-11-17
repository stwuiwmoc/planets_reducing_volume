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
