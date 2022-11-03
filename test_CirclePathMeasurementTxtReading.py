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

    original_filepath = "raw_data/220117xrmEAi.v5.41.hei.txt"
    deformed_filepath = "raw_data/220117xrmFEi.v5.41.hei.txt"

    diff = pom.CirclePathMeasurementTxtReading(
        Constants=CONSTS,
        original_txt_fpath=original_filepath,
        deformed_txt_fpath=deformed_filepath
    )

    fig1 = plt.figure()
    gs1 = fig1.add_gridspec(1, 1)

    ax11 = fig1.add_subplot(gs1[0, 0])
    ax11.plot(diff.df_diff["degree"], diff.df_diff["height"])

    ax11.set_xlabel("degree")
    ax11.set_ylabel("height [m]")
    ax11.grid()

    fig1.tight_layout()
