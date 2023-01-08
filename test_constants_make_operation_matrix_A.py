# %%
import importlib

import matplotlib.pyplot as plt
import numpy as np

import planets_optimize_myclass as pom


def make_zernike_plot(fig, position, zernike_array):
    ax = fig.add_subplot(position)
    ax.plot(
        np.arange(0, len(zernike_array)) + 1,
        zernike_array,
        marker="s"
    )
    ax.grid()
    return ax


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
    torque_num_a = 6  # 6, 18, 30 のどれか
    torque_num_b = 12  # 12, 24, 36 のどれか
    zernike_num = 6  # 6, 17, 28 のどれか

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

    a_n_array = CONSTS.operation_matrix_A[:, zernike_num - 1]
    zernike_removed_a_n_array = a_n_array.copy()
    zernike_removed_a_n_array[0] = zernike_removed_a_n_array[1] = zernike_removed_a_n_array[2] = 0

    fig1 = plt.figure(figsize=(10, 10))
    gs1 = fig1.add_gridspec(4, 2)

    ax11 = make_zernike_plot(fig1, gs1[0, 0], surface_a.zernike_value_array)
    ax11.set_title("M" + str(torque_num_a).zfill(2))

    ax12 = make_zernike_plot(fig1, gs1[1, 0], surface_b.zernike_value_array)
    ax12.set_title("M" + str(torque_num_b).zfill(2))

    ax13 = make_zernike_plot(fig1, gs1[2, 0], surface_ab.zernike_value_array)
    ax13.set_title("M" + str(torque_num_a).zfill(2) + " - M" + str(torque_num_b).zfill(2))

    ax14 = make_zernike_plot(fig1, gs1[3, 0], a_n_array)
    ax14.set_title("a_" + str(zernike_num) + " from CONSTS.operation_matrix_A")

    ax15 = make_zernike_plot(fig1, gs1[0, 1], zernike_removed_surface_a.zernike_value_array)

    ax16 = make_zernike_plot(fig1, gs1[1, 1], zernike_removed_surface_b.zernike_value_array)

    ax17 = make_zernike_plot(fig1, gs1[2, 1], zernike_removed_surface_ab.zernike_value_array)

    ax18 = make_zernike_plot(fig1, gs1[3, 1], zernike_removed_a_n_array)

    fig1.tight_layout()
