# %%
import importlib

import matplotlib.pyplot as plt

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

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=25e-3,
        pixel_number=1024,
        zernike_max_degree=10,
        offset_height_percent=2)

    original_filepath = "raw_data/220117xrmEAi.v5.40.hei.txt"
    deformed_filepath = "raw_data/220117xrmFEi.v5.40.hei.txt"

    ignore_zernike_number_list = [1, 2, 3, 4, 5, 6, 7, 8, 11]

    measurement_diff = pom.CirclePathMeasurementTxtReading(
        Constants=CONSTS,
        original_txt_fpath=original_filepath,
        deformed_txt_fpath=deformed_filepath
    )

    residual = pom.CirclePathZernikeFitting(
        Constants=CONSTS,
        circle_path_radius=870e-3,
        degree_array=measurement_diff.df_diff["degree"].values,
        unprocessed_height_array=measurement_diff.df_diff["height"].values,
        ignore_zernike_number_list=ignore_zernike_number_list
    )

    # plot

    fig1 = plt.figure(figsize=(10, 7))
    gs1 = fig1.add_gridspec(2, 1)

    ax11 = fig1.add_subplot(gs1[0, 0])
    ax11.plot(
        measurement_diff.df_diff["degree"],
        measurement_diff.df_diff["height"],
        label="raw")
    ax11.plot(
        residual.degree_array,
        residual.removing_zernike_height_array,
        label="zernike fit")

    ax11.set_xlabel("degree")
    ax11.set_ylabel("height [m]")
    ax11.grid()
    ax11.legend()

    ax12 = fig1.add_subplot(gs1[1, 0])
    ax12.plot(residual.degree_array, residual.zernike_removed_height_array)
    ax12.set_xlabel("degree")
    ax12.set_ylabel("zernike removed height [m]")
    ax12.grid()

    fig1.tight_layout()
