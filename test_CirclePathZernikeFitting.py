# %%
import importlib

import matplotlib.pyplot as plt

import planets_optimize_myclass as pom

if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=25e-3,
        pixel_number=1024,
        zernike_max_degree=10,
        offset_height_percent=2)

    original_fpath = "raw_data/220117xrmEAi.v5.60.hei.txt"
    deformed_fpath = "raw_data/220117xrmFEi.v5.60.hei.txt"

    ignore_zernike_number_list = [1, 2, 3, 4, 5, 6, 7, 8, 11]

    measurement = pom.CirclePathMeasurementTxtReading(
        Constants=CONSTS,
        original_txt_fpath=original_fpath,
        deformed_txt_fpath=deformed_fpath
    )

    removed = pom.CirclePathZernikeFitting(
        Constants=CONSTS,
        circle_path_radius=measurement.circle_path_radius,
        degree_array=measurement.df_diff["degree"],
        unprocessed_height_array=measurement.df_diff["height"],
        ignore_zernike_number_list=ignore_zernike_number_list
    )

    fig1 = plt.figure(figsize=(10, 10))
    gs1 = fig1.add_gridspec(3, 1)

    # WH変形前後の差分と、そこに対するzernike fitting
    ax11 = fig1.add_subplot(gs1[0, 0])
    ax11.plot(
        measurement.df_diff["degree"],
        measurement.df_diff["height"],
        label="WH deform height")
    ax11.plot(
        removed.degree_array,
        removed.removing_zernike_height_array,
        label="zernike fit"
    )
    ax11.grid()
    ax11.set_ylabel("height [m]")
    ax11.legend()

    ax12 = fig1.add_subplot(gs1[1, 0])
    ax12.plot(
        removed.degree_array,
        removed.zernike_removed_height_array
    )
    ax12.grid()
    ax12.set_xlabel("angle [deg]")
    ax12.set_ylabel("height [m]")

    ax13 = fig1.add_subplot(gs1[2, 0])
    ax13.plot(
        removed.r_const_zernike_number_meaning_list,
        removed.r_const_zernike_polynomial_array
    )
    ax13.grid()
    ax13.set_xlabel("zernike number in radius const")
    ax13.set_ylabel("zernike polynomial value [m]")

    fig1.tight_layout()
