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

    serial_0 = "0117e"
    serial_n = "0117l"
    input_filename_0 = "../sag_integration_code/" + "mkfolder/psm_test_kagi_data_integration/" + serial_0 + "_height.csv"
    input_filename_n = "../sag_integration_code/" + "mkfolder/psm_test_kagi_data_integration/" + serial_n + "_height.csv"

    mes0n = pom.CirclePathMeasurementReading(Constants=CONSTS,
                                             original_csv_fpath=input_filename_0,
                                             deformed_csv_fpath=input_filename_n)

    res0n = pom.CirclePathZernikeFitting(Constants=CONSTS,
                                         circle_path_radius=870e-3,
                                         df_diff=mes0n.df_diff,
                                         ignore_zernike_number_list=[1, 2, 3, 4, 5, 6, 7, 8, 11])

    fig1 = plt.figure(figsize=(7, 7))
    gs1 = fig1.add_gridspec(2, 1)

    ax11 = fig1.add_subplot(gs1[0, 0])
    ax11.plot(res0n.df_diff["degree"], res0n.df_diff["height"], label="height_transform")
    ax11.plot(res0n.df_diff["degree"], res0n.removing_zernike, label="zernike_fitting")
    ax11.grid()
    ax11.legend()

    ax12 = fig1.add_subplot(gs1[1, 0])
    ax12.plot(res0n.df_diff["degree"], res0n.height_removed)
    ax12.grid()
    ax12.set_xlabel("robot-arm theta [deg]")

    fig1.tight_layout()

    fig2 = plt.figure(figsize=(7, 4))
    gs2 = fig2.add_gridspec(1, 1)
    fig2.suptitle(serial_n + " - " + serial_0)

    ax21 = fig2.add_subplot(gs2[0, 0])
    ax21.plot(res0n.df_diff["degree"], res0n.height_removed * 1e9)
    ax21.grid()
    ax21.set_xlabel("robot-arm theta [deg]")
    ax21.set_ylabel("height_deform [nm]")

    fig2.savefig(mkfolder() + serial_n + "_fig2.png")
