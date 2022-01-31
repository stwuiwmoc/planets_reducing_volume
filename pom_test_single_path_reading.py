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

    CONSTS = pom.Constants(physical_radius=925e-3,
                           ignore_radius=25e-3,
                           pixel_number=1024,
                           zernike_max_degree=10,
                           offset_height_percent=2)

    serial_n = "0117l"
    input_filename_n = "../sag_integration_code/" + "mkfolder/psm_test_kagi_data_integration/" + serial_n + "_height.csv"

    mes = pom.CirclePathMeasurementReading(Constants=CONSTS,
                                           original_csv_fpath=input_filename_n,
                                           deformed_csv_fpath="")

    res = pom.CirclePathZernikeFitting(Constants=CONSTS,
                                       circle_path_radius=870e-3,
                                       df_diff=mes.df_diff,
                                       ignore_zernike_number_list=[1, 2, 3, 4, 5, 6, 7, 8, 11])

    circumference = -res.df_diff["radian"] * res.circle_path_radius * 1e3

    fig1 = plt.figure(figsize=(10, 5))
    gs1 = fig1.add_gridspec(2, 1)
    fig1.suptitle(input_filename_n[-16:-11])

    ax11 = fig1.add_subplot(gs1[0, 0])
    ax11.plot(circumference, res.df_diff["height"], label="height")
    ax11.plot(circumference, res.removing_zernike, label="zernike_fitting")
    ax11.grid()
    ax11.legend()
    ax11.set_ylabel("height [m]")

    ax12 = fig1.add_subplot(gs1[1, 0])
    ax12.plot(circumference, res.height_removed)
    ax12.grid()
    ax12.set_xlabel("robot-arm scanning coordinate [mm]")

    fig1.tight_layout()
    fig1.savefig(mkfolder() + serial_n + "_fig1.png")
