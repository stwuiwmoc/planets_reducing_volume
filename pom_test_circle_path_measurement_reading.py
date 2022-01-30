# %%
import importlib
import matplotlib.pyplot as plt

import planets_optimize_myclass as pom

if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(physical_radius=925e-3,
                           ignore_radius=25e-3,
                           pixel_number=1024,
                           zernike_max_degree=10,
                           offset_height_percent=2)

    input_filename0 = "../sag_integration_code/" + "mkfolder/psm_test_kagi_data_integration/0117e_height.csv"
    input_filename3 = "../sag_integration_code/" + "mkfolder/psm_test_kagi_data_integration/0117h_height.csv"

    mes03 = pom.CirclePathMeasurementReading(Constants=CONSTS,
                                             original_csv_fpath=input_filename0,
                                             deformed_csv_fpath=input_filename3)

    res03 = pom.CirclePathZernikeFitting(Constants=CONSTS,
                                         circle_path_radius=870e-3,
                                         df_diff=mes03.df_diff,
                                         ignore_zernike_number_list=[1, 2, 3, 4, 5, 6, 7, 8, 11])

    fig1 = plt.figure(figsize=(7, 7))
    gs1 = fig1.add_gridspec(2, 1)

    ax11 = fig1.add_subplot(gs1[0, 0])
    ax11.plot(res03.df_diff["degree"], res03.df_diff["height"], label="height_diff")
    ax11.plot(res03.df_diff["degree"], res03.removing_zernike, label="zernike_fitting")
    ax11.grid()
    ax11.legend()

    ax12 = fig1.add_subplot(gs1[1, 0])
    ax12.plot(res03.df_diff["degree"], res03.height_removed)
    ax12.grid()
    ax12.set_xlabel("robot-arm theta [deg]")

    fig1.tight_layout()
