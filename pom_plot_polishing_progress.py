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
        ignore_radius=50e-3,
        pixel_number=1024,
        zernike_max_degree=11,
        offset_height_percent=2
    )

    cbar_min_percent_ = 40
    cbar_max_percent_ = 90
    zaxis_bottom_percent_ = 20

    mesE = pom.ExelisCsvToSurface(
        constants=CONSTS)

    mes0 = pom.KagiStitchToSurface(
        constants=CONSTS,
        txt_fpath="raw_data/0923xm130_1007ym830.hei.v2_dense.txt")

    mesN = pom.KagiStitchToSurface(
        constants=CONSTS,
        txt_fpath="raw_data/0228xm130_3degintAll_0301ym870-510CirAll_0228ykagomeACc.v4.6.hei_dense.txt"
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

    zerN = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=mesN.surface,
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

    filterN = pom.FilteredSurface(
        constants=CONSTS,
        inputed_surface=zerN.surface,
        filter_parameter=50
    )

    fig1 = plt.figure(figsize=(14, 14))
    gs1 = fig1.add_gridspec(3, 3)

    ax19 = mesE.make_image_plot(
        fig1, gs1[2, 0],
        mes0.surface, cbar_min_percent_, cbar_max_percent_)
    ax11 = mes0.make_image_plot(
        fig1, gs1[2, 1],
        None, cbar_min_percent_, cbar_max_percent_)
    ax12 = mesN.make_image_plot(
        fig1, gs1[2, 2],
        mes0.surface, cbar_min_percent_, cbar_max_percent_)

    ax110 = zerE.make_image_plot(
        fig1, gs1[1, 0],
        zer0.surface, cbar_min_percent_, cbar_max_percent_)
    ax15 = zer0.make_image_plot(
        fig1, gs1[1, 1],
        None, cbar_min_percent_, cbar_max_percent_)
    ax16 = zerN.make_image_plot(
        fig1, gs1[1, 2],
        zer0.surface, cbar_min_percent_, cbar_max_percent_)

    ax111 = filterE.make_3d_plot(
        fig1, gs1[0, 0],
        zer0.surface, cbar_min_percent_, cbar_max_percent_, zaxis_bottom_percent_)
    ax111.set_title("Exelis")

    ax17 = filter0.make_3d_plot(
        fig1, gs1[0, 1],
        zer0.surface, cbar_min_percent_, cbar_max_percent_, zaxis_bottom_percent_)
    ax17.set_title("Before Polish")

    ax18 = filterN.make_3d_plot(
        fig1, gs1[0, 2],
        zer0.surface, cbar_min_percent_, cbar_max_percent_, zaxis_bottom_percent_)
    ax18.set_title("After Polish")

    fig1.tight_layout()
    fig1.savefig(mkfolder() + "fig1.png")

    fig2 = plt.figure(figsize=(7, 10))
    gs2 = fig2.add_gridspec(3, 1)

    ax21 = mesE.make_zernike_value_plot(figure=fig2, position=gs2[0, 0])
    ax21.set_ylim(mesE.zernike_value_array.min(), mesE.zernike_value_array.max())

    ax22 = mes0.make_zernike_value_plot(figure=fig2, position=gs2[1, 0])
    ax22.set_ylim(mesE.zernike_value_array.min(), mesE.zernike_value_array.max())

    ax23 = mesN.make_zernike_value_plot(figure=fig2, position=gs2[2, 0])
    ax23.set_ylim(mesE.zernike_value_array.min(), mesE.zernike_value_array.max())

    fig2.tight_layout()
    fig2.savefig(mkfolder() + "fig2.png")
