# %%
import importlib

import matplotlib.pyplot as plt
import numpy as np

import planets_optimize_myclass as pom

if __name__ == "__main__":
    importlib.reload(pom)

    CONSTS = pom.Constants(
        physical_radius=925e-3,
        ignore_radius=25e-3,
        pixel_number=1024,
        zernike_max_degree=11,
        offset_height_percent=2,
        alpha_array=np.array([
            -3.7, +3.7, -2.0, -2.0, +3.8, -3.8,
            -3.7, +3.7, -2.0, -2.0, +3.8, -3.8,
            -3.7, +3.7, -2.0, -2.0, +3.8, -3.8,
            -3.7, +3.7, -2.0, -2.0, +3.8, -3.8,
            -3.7, +3.7, -2.0, -2.0, +3.8, -3.8,
            -3.7, +3.7, -2.0, -2.0, +3.8, -3.8,
        ])
    )

    # %%
    importlib.reload(pom)
    torque_num_a = 6
    torque_num_b = 12

    torque_value_array_a = np.zeros(36)
    torque_value_array_a[torque_num_a - 1] = +1

    torque_value_array_b = np.zeros(36)
    torque_value_array_b[torque_num_b - 1] = +1

    surface_a = pom.TorqueToSurface(
        constants=CONSTS,
        torque_value_array=torque_value_array_a
    )

    surface_b = pom.TorqueToSurface(
        constants=CONSTS,
        torque_value_array=torque_value_array_b
    )

    zernike_removed_surface_a = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=surface_a.surface,
        removing_zernike_number_list=[1, 2, 3]
    )

    zernike_removed_surface_b = pom.ZernikeRemovedSurface(
        constants=CONSTS,
        inputed_surface=surface_b.surface,
        removing_zernike_number_list=[1, 2, 3]
    )

    surface_ab = pom.Surface(
        constants=CONSTS,
        surface=surface_a.surface - surface_b.surface
    )

    zernike_removed_surface_ab = pom.Surface(
        constants=CONSTS,
        surface=zernike_removed_surface_a.surface - zernike_removed_surface_b.surface
    )

    # plot
    fig1 = plt.figure(figsize=(10, 10))
    gs1 = fig1.add_gridspec(3, 3)

    # zernike除去前の形状プロット
    ax11 = surface_a.make_image_plot(
        figure=fig1, position=gs1[0, 0]
    )
    ax11.set_title("M" + str(torque_num_a).zfill(2) + "\n" + ax11.get_title())

    ax12 = surface_b.make_image_plot(
        figure=fig1, position=gs1[1, 0]
    )
    ax12.set_title("M" + str(torque_num_b).zfill(2) + "\n" + ax12.get_title())

    ax13 = surface_ab.make_image_plot(
        figure=fig1, position=gs1[2, 0]
    )
    ax13.set_title("M" + str(torque_num_a).zfill(2) + " - M" + str(torque_num_b).zfill(2) + "\n" + ax13.get_title())

    # zernike除去後の形状プロット
    ax14 = zernike_removed_surface_a.make_image_plot(
        figure=fig1, position=gs1[0, 1]
    )
    ax14.set_title("M" + str(torque_num_a).zfill(2) + " z1, 2, 3 removed\n" + ax14.get_title())

    ax15 = zernike_removed_surface_b.make_image_plot(
        figure=fig1, position=gs1[1, 1]
    )
    ax15.set_title("M" + str(torque_num_b).zfill(2) + " z1, 2, 3 removed\n" + ax15.get_title())

    ax16 = zernike_removed_surface_ab.make_image_plot(
        figure=fig1, position=gs1[2, 1]
    )
    ax16.set_title("z1, 2, 3 removed\n" + ax16.get_title())

    # zernike係数ベクトルのプロット
    ax17 = fig1.add_subplot(gs1[0, 2])
    ax17.plot(
        np.arange(CONSTS.zernike_max_degree) + 1,
        surface_a.zernike_value_array,
        marker="s"
    )
    ax17.grid()

    ax18 = fig1.add_subplot(gs1[1, 2])
    ax18.plot(
        np.arange(CONSTS.zernike_max_degree) + 1,
        surface_b.zernike_value_array,
        marker="s"
    )
    ax18.grid()

    ax19 = fig1.add_subplot(gs1[2, 2])
    ax19.plot(
        np.arange(CONSTS.zernike_max_degree) + 1,
        surface_ab.zernike_value_array,
        marker="s"
    )
    ax19.grid()

    fig1.tight_layout()
