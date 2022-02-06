# %%

import matplotlib.pyplot as plt
import importlib

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
        offset_height_percent=2
    )

    mesE = pom.StitchedCsvToSurface(
        constants=CONSTS,
        original_stitched_csv_fpath="mkfolder/stitch2mesh/zer03_exelis.csv",
        deformed_stitched_csv_fpath=""
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

    zerE = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=mesE.surface,
        removing_zernike_number_list=[1, 2, 3, 4, 5, 6]
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

    filterE = pom.FilteredSurface(
        constants=CONSTS,
        inputed_surface=zerE.surface,
        filter_parameter=50
    )

    filter0 = pom.FilteredSurface(
        constants=CONSTS,
        inputed_surface=zer0.surface,
        filter_parameter=50
    )

    filter8 = pom.FilteredSurface(
        constants=CONSTS,
        inputed_surface=zer8.surface,
        filter_parameter=50
    )

    cbar_min_percent_ = 35
    cbar_max_percent_ = 90

    fig1 = plt.figure(figsize=(14, 14))
    gs1 = fig1.add_gridspec(3, 3)

    ax19 = mesE.make_image_plot(
        fig1, gs1[0, 0],
        mes0.surface, cbar_min_percent_, cbar_max_percent_)
    ax19.set_title("exelis\n" + ax19.get_title())
    ax11 = mes0.make_image_plot(
        fig1, gs1[0, 1],
        None, cbar_min_percent_, cbar_max_percent_)
    ax11.set_title("before_polish\n" + ax11.get_title())
    ax12 = mes8.make_image_plot(
        fig1, gs1[0, 2],
        mes0.surface, cbar_min_percent_, cbar_max_percent_)
    ax12.set_title("after 8th polish\n" + ax12.get_title())

    ax110 = zerE.make_image_plot(
        fig1, gs1[1, 0],
        zer0.surface, cbar_min_percent_, cbar_max_percent_)
    ax15 = zer0.make_image_plot(
        fig1, gs1[1, 1],
        None, cbar_min_percent_, cbar_max_percent_)
    ax16 = zer8.make_image_plot(
        fig1, gs1[1, 2],
        zer0.surface, cbar_min_percent_, cbar_max_percent_)

    ax111 = filterE.make_3d_plot(
        fig1, gs1[2, 0],
        filter0.surface, cbar_min_percent_, cbar_max_percent_)
    ax17 = filter0.make_3d_plot(
        fig1, gs1[2, 1],
        None, cbar_min_percent_, cbar_max_percent_)
    ax18 = filter8.make_3d_plot(
        fig1, gs1[2, 2],
        filter0.surface, cbar_min_percent_, cbar_max_percent_)

    fig1.tight_layout()
    fig1.savefig(mkfolder() + "fig1.png")
