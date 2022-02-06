# %%

import matplotlib.pyplot as plt
import importlib

import planets_optimize_myclass as pom

if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=25e-3,
        pixel_number=1024,
        zernike_max_degree=10,
        offset_height_percent=2
    )

    mes0 = pom.StitchedCsvToSurface(
        constants=CONSTS,
        original_stitched_csv_fpath="mkfolder/stitch2mesh/zer03_0923xm130_1007ym830.hei.v2_dense.csv",
        deformed_stitched_csv_fpath=""
    )

    mes8 = pom.StitchedCsvToSurface(
        constants=CONSTS,
        original_stitched_csv_fpath="mkfolder/stitch2mesh/zer03_0131xm130allcc_0201cirAll.v4.8.hei_dense.csv",
        deformed_stitched_csv_fpath=""
    )

    zer0 = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=mes0.surface,
        removing_zernike_number_list=[1, 2, 3, 4, 5, 6]
    )

    zer8 = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=mes8.surface,
        removing_zernike_number_list=[1, 2, 3, 4, 5, 6]
    )

    filtered0 = pom.FilteredSurface(
        constants=CONSTS,
        inputed_surface=zer0.surface,
        filter_parameter=50
    )

    filtered8 = pom.FilteredSurface(
        constants=CONSTS,
        inputed_surface=zer8.surface,
        filter_parameter=50
    )

    cbar_min_percent_ = 35
    cbar_max_percent_ = 90

    fig1 = plt.figure(figsize=(10, 14))
    gs1 = fig1.add_gridspec(3, 2)

    ax11 = mes0.make_image_plot(
        fig1, gs1[0, 0],
        None, cbar_min_percent_, cbar_max_percent_)
    ax12 = mes8.make_image_plot(
        fig1, gs1[0, 1],
        mes0.surface, cbar_min_percent_, cbar_max_percent_)
    ax15 = zer0.make_image_plot(
        fig1, gs1[1, 0],
        None, cbar_min_percent_, cbar_max_percent_)
    ax16 = zer8.make_image_plot(
        fig1, gs1[1, 1],
        zer0.surface, cbar_min_percent_, cbar_max_percent_)
    ax17 = filtered0.make_3d_plot(
        fig1, gs1[2, 0],
        None, cbar_min_percent_, cbar_max_percent_)
    ax18 = filtered8.make_3d_plot(
        fig1, gs1[2, 1],
        filtered0.surface, cbar_min_percent_, cbar_max_percent_)

    fig1.tight_layout()
