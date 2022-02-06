# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 14:59:35 2021

@author: swimc
"""
# %%
from numpy import ndarray
import numpy as np
import proper as pr
import matplotlib.pyplot as plt
import scipy as sp
from scipy import optimize
import cv2
import PIL
import time
import pandas as pd

from matplotlib import cm
from matplotlib.colors import Normalize
import mpl_toolkits.axes_grid1


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


def mkhelp(instance):
    import inspect
    attr_list = list(instance.__dict__.keys())
    for attr in attr_list:
        if attr.startswith("_"):
            continue
        print(attr)
    for method in inspect.getmembers(instance, inspect.ismethod):
        if method[0].startswith("_"):
            continue
        print(method[0] + "()")


def make_meshgrid(x_min, x_max, y_min, y_max, pixel_number):
    x_array = np.linspace(x_min, x_max, pixel_number)
    y_array = np.linspace(y_min, y_max, pixel_number)
    return np.meshgrid(x_array, y_array)


def make_remaining_matrix(matrix, ignore_zernike_number_list):
    idx_array = np.array(ignore_zernike_number_list) - 1
    remaining_matrix = np.delete(arr=matrix, obj=idx_array, axis=0)
    return remaining_matrix


def oap_calculation(radius_of_curvature, off_axis_distance, clocking_angle_rad, x_mesh, y_mesh):
    """

    Parameters
    ----------
    radius_of_curvature : float [m]
        曲率半径
    off_axis_distance : floot [m]
        軸外し距離 (off-axis)
    clocking_angle_rad : float [rad]
        回転角

    x_mesh : 2D-mesh-array [m]
        軸外し方向。 off_axis_distanceを増加させる方向をx_meshの正の向きとして定義
    y_mesh : 2D-mesh-array [m]

    Returns
    -------
    oap_height : 2D-mesh-array [m]
        input された clocking_angle_rad, radius_of_curvature, off_axis_distance での切り取り放物面
    """
    phi = clocking_angle_rad

    # [m] --> [mm] に変換
    p = 1e3 * radius_of_curvature
    aqx = 1e3 * off_axis_distance
    cx = 1e3 * x_mesh
    cy = 1e3 * y_mesh

    a = 1 / (2 * p)
    aqz = a * (aqx ** 2 + 0 ** 2)
    theta = np.arctan(2 * a * aqx)

    D = -4 * a**2 * cx**2 * np.sin(phi)**2 * np.sin(theta)**2 - 8 * a**2 * cx * cy * np.sin(phi) * np.sin(theta)**2 * np.cos(phi) + 4 * a**2 * cy**2 * np.sin(phi)**2 * np.sin(theta)**2 - 4 * a**2 * cy**2 * np.sin(theta)**2 + 2 * a * aqx * np.sin(2 * theta) + 4 * a * aqz * np.sin(theta)**2 + 4 * a * cx * np.sin(theta) * np.cos(phi) - 4 * a * cy * np.sin(phi) * np.sin(theta) - np.sin(theta)**2 + 1

    cz1 = (4 * a * aqx * np.sin(theta) - a * cx * (np.sin(phi - 2 * theta) - np.sin(phi + 2 * theta)) - a * cy * (np.cos(phi - 2 * theta) - np.cos(phi + 2 * theta)) - 2 * np.sqrt(D) + 2 * np.cos(theta)) / (4 * a * np.sin(theta)**2)
    oap_height = cz1 * 1e-3  # [mm] --> [m] に変換
    return oap_height


def zernike_term_calculation(term_number: int,
                             radius: ndarray,
                             theta: ndarray) -> ndarray:
    """zernike_term_calculation
    zernikeの各項を計算する（係数は入ってない）

    Parameters
    ----------
    term_number : int
        計算したいzernikeの項の番号（1~11）
    radius : ndarray
        半径 [m]
    theta : ndarray
        角度 [rad]

    Returns
    -------
    ndarray
        zernikeのどれか一つの項の計算結果（係数は入ってない）
    """
    if term_number == 1:
        return 1
    if term_number == 2:
        return 2 * (radius) * np.cos(theta)
    if term_number == 3:
        return 2 * (radius) * np.sin(theta)
    if term_number == 4:
        return np.sqrt(3) * (2.0 * radius ** 2 - 1.0)
    if term_number == 5:
        return np.sqrt(6) * (radius ** 2) * np.sin(2 * theta)
    if term_number == 6:
        return np.sqrt(6) * (radius ** 2) * np.cos(2 * theta)
    if term_number == 7:
        return np.sqrt(8) * (3.0 * radius ** 3 - 2.0 * radius) * np.sin(theta)
    if term_number == 8:
        return np.sqrt(8) * (3.0 * radius ** 3 - 2.0 * radius) * np.cos(theta)
    if term_number == 9:
        return np.sqrt(8) * (radius ** 3) * np.sin(3 * theta)
    if term_number == 10:
        return np.sqrt(8) * (radius ** 3) * np.cos(3 * theta)
    if term_number == 11:
        return np.sqrt(5) * (6.0 * radius ** 4 - 6.0 * radius ** 2 + 1.0)
    else:
        print("term_number must be 1~11")
        return


def zernike_polynomial_calculation(coef: list[float],
                                   radius: ndarray,
                                   theta: ndarray) -> ndarray:
    """
    zernike多項式のarrayを計算する（係数込み）

    Parameters
    ----------
    coef : list[float]
        zernike係数ベクトル
    radius : ndarray
        半径 [m]
    theta : ndarray
        角度 [rad]

    Returns
    -------
    ndarray
        zernike多項式のarray
    """

    zernike1 = coef[0] * zernike_term_calculation(1, radius, theta)
    zernike2 = coef[1] * zernike_term_calculation(2, radius, theta)
    zernike3 = coef[2] * zernike_term_calculation(3, radius, theta)
    zernike4 = coef[3] * zernike_term_calculation(4, radius, theta)
    zernike5 = coef[4] * zernike_term_calculation(5, radius, theta)
    zernike6 = coef[5] * zernike_term_calculation(6, radius, theta)
    zernike7 = coef[6] * zernike_term_calculation(7, radius, theta)
    zernike8 = coef[7] * zernike_term_calculation(8, radius, theta)
    zernike9 = coef[8] * zernike_term_calculation(9, radius, theta)
    zernike10 = coef[9] * zernike_term_calculation(10, radius, theta)
    zernike11 = coef[10] * zernike_term_calculation(11, radius, theta)

    zernike_polynomial = zernike1 + zernike2 + zernike3 + zernike4 + zernike5 + zernike6 + zernike7 + zernike8 + zernike9 + zernike10 + zernike11
    return zernike_polynomial


class Constants:

    def __init__(self, physical_radius, ignore_radius,
                 pixel_number, zernike_max_degree, offset_height_percent):
        """
        class : constants

        Parameters
        ----------
        physical_radius : float
            [m] physical radius of the mirror
        ignore_radius : float
            [m] outer mask
        pixel_number : int
            vertical and horizontal pixel number
        zernike_max_degree : int
            max zernike number
        offset_height_percent : float
            ignore height in percent. if you set 2, the lower 2% is ignored in volume calculation

        Returns
        -------
        None.

        """
        self.physical_radius = physical_radius
        self.ignore_radius = ignore_radius
        self.pixel_number = pixel_number
        self.offset_height_percent = offset_height_percent

        self.varid_radius = physical_radius - ignore_radius
        self.xx, self.yy = make_meshgrid(-physical_radius, physical_radius,
                                         -physical_radius, physical_radius,
                                         pixel_number)
        self.tf = np.where(self.xx ** 2 + self.yy ** 2 <= self.varid_radius**2, True, False)
        self.mask = np.where(self.tf, 1, np.nan)
        self.zernike_max_degree = zernike_max_degree
        self.operation_matrix = np.genfromtxt("raw_data/WT06_zer10_operation_matrix[m].csv", delimiter=",").T

    def h(self):
        mkhelp(self)


class Surface:
    def __init__(self, constants, surface):
        self.consts = constants
        self.surface = surface

        self.pv = self._pv_calculation()
        self.rms = self._rms_calculation()
        self.volume = self._volume_calculation()[0]
        self.offset_height_value = self._volume_calculation()[1]

    def h(self):
        mkhelp(self)

    def _pv_calculation(self):
        array2d = self.surface
        peak = np.nanmax(array2d)
        valley = np.nanmin(array2d)
        pv = peak - valley
        return pv

    def _rms_calculation(self):
        array2d = self.surface
        sigma = np.nansum(array2d**2)
        data_count = np.sum(~np.isnan(array2d))
        rms = np.sqrt(sigma / data_count)
        return rms

    def _volume_calculation(self):
        # 下位○%は体積計算で無視するために、下位○%の閾値を計算
        sorted_surface = np.sort(self.surface.flatten())  # np.nanは最大値の後に並ぶ
        value_count = np.sum(~np.isnan(self.surface))
        offset_height_idx = int(
            value_count *
            self.consts.offset_height_percent /
            100)
        offset_height_value = sorted_surface[offset_height_idx]

        lower_ignored_surface = np.where(self.surface >= offset_height_value,
                                         self.surface - offset_height_value,
                                         0)

        # 1pixelあたりの単位面積を計算
        physical_diameter = 2 * self.consts.physical_radius
        unit_pixel_area = (physical_diameter / self.consts.pixel_number)**2

        # 1 pixelあたりの体積を計算
        unit_pixel_volume = unit_pixel_area * lower_ignored_surface
        volume_in_m3 = unit_pixel_volume.sum()
        return (volume_in_m3, offset_height_value)

    def make_image_plot(self, figure, position=111,
                        color_scale=False, cbar_min_percent=0, cbar_max_percent=100,
                        pv_digits=2, rms_digits=2):

        cmap = cm.jet
        fontsize = 15
        title = "pv = " + str(round(self.pv * 1e6, pv_digits)) + " [um]" + "\n" + "RMS = " + str(round(self.rms * 1e6, rms_digits)) + " [um]"

        cbar_min = np.nanmin(self.surface) + self.pv * cbar_min_percent / 100
        cbar_max = np.nanmin(self.surface) + self.pv * cbar_max_percent / 100

        extent = [-self.consts.physical_radius, self.consts.physical_radius,
                  -self.consts.physical_radius, self.consts.physical_radius]

        ax = figure.add_subplot(position)
        ax.imshow(self.surface, interpolation="nearest", cmap=cmap,
                  vmin=cbar_min, vmax=cbar_max, origin="lower", extent=extent)
        ax.set_title(title, fontsize=fontsize)

        divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
        cax = divider.append_axes("right", "5%", pad="3%")

        norm = Normalize(vmin=cbar_min, vmax=cbar_max)
        cbar_title = "[m]"
        mappable = cm.ScalarMappable(norm=norm, cmap=cmap)

        cbar = figure.colorbar(mappable, ax=ax, cax=cax)
        cbar.set_label(cbar_title, fontsize=fontsize)
        return ax

    def make_circle_path_plot(self,
                              figure,
                              position=111,
                              radius=0.870,
                              height_magn=1e9,
                              height_unit_str="[nm]"):

        fontsize = 15
        varid_radius_pixel_number = int(self.consts.varid_radius / self.consts.physical_radius * self.consts.pixel_number / 2)
        measurement_radius_idx = int(radius * 1e3)

        image = np.where(self.consts.tf, self.surface, 0)
        flags = cv2.INTER_CUBIC + cv2.WARP_FILL_OUTLIERS + cv2.WARP_POLAR_LINEAR

        linear_polar_image = cv2.warpPolar(src=image,
                                           dsize=(int(self.consts.varid_radius * 1e3), 360),
                                           center=(self.consts.pixel_number / 2, self.consts.pixel_number / 2),
                                           maxRadius=varid_radius_pixel_number,
                                           flags=flags)

        circle_path_line = height_magn * linear_polar_image[:, measurement_radius_idx]

        geometric_degree = np.arange(360)  # CCW、pomでのx軸+方向を0degとする角度
        robot_degree = geometric_degree - 90 - 10.814
        robot_degree = np.where(robot_degree > 0,
                                robot_degree,
                                robot_degree + 360)

        height_pv = np.nanmax(circle_path_line) - np.nanmin(circle_path_line)
        height_pv_str = str(round(height_pv, 2))

        ax = figure.add_subplot(position)
        ax.plot(robot_degree, circle_path_line)
        ax.grid()

        ax.set_title("circle_path ( pv = " + height_pv_str + " " + height_unit_str + " )",
                     fontsize=fontsize)
        ax.set_xlabel("degree", fontsize=fontsize)
        ax.set_ylabel("heignt " + height_unit_str + "\nat R=" + str(measurement_radius_idx),
                      fontsize=fontsize)
        return ax

    def _make_masked_zernike_surface(self, zernike_number_list, zernike_value_array):
        optical_wavelength = 500e-9

        wavestruct = pr.prop_begin(beam_diameter=2 * self.consts.varid_radius,
                                   lamda=optical_wavelength,
                                   grid_n=self.consts.pixel_number,
                                   beam_diam_fraction=1)

        wfe = pr.prop_zernikes(wavestruct,
                               zernike_number_list,
                               zernike_value_array)

        masked_wfe = self.consts.mask * wfe
        return masked_wfe


class ZernikeRemovedSurface(Surface):
    def __init__(self, constants, inputed_surface, removing_zernike_number_list):
        self.consts = constants
        self.removing_zernike_number_list = removing_zernike_number_list

        full_zernike_number_list = [i + 1 for i in range(self.consts.zernike_max_degree)]

        self.inputed_surface = inputed_surface
        self.zernike_value_array = self._zernike_value_array_calculation(self.inputed_surface)

        ignore_zernike_number_list = make_remaining_matrix(full_zernike_number_list,
                                                           self.removing_zernike_number_list)
        self.removing_zernike_value_array = make_remaining_matrix(self.zernike_value_array,
                                                                  ignore_zernike_number_list)

        self.removing_surface = super()._make_masked_zernike_surface(self.removing_zernike_number_list,
                                                                     self.removing_zernike_value_array)
        self.surface = self.inputed_surface - self.removing_surface

        self.pv = super()._pv_calculation()
        self.rms = super()._rms_calculation()
        self.volume = super()._volume_calculation()[0]
        self.offset_height_value = super()._volume_calculation()[1]

    def h(self):
        mkhelp(self)

    def _zernike_value_array_calculation(self, surface):
        surface_without_nan = np.where(self.consts.tf,
                                       surface,
                                       0)

        zernike_value_array = pr.prop_fit_zernikes(wavefront0=surface_without_nan,
                                                   pupil0=self.consts.tf,
                                                   pupilradius0=self.consts.pixel_number // 2,
                                                   nzer=self.consts.zernike_max_degree,
                                                   xc=self.consts.pixel_number // 2,
                                                   yc=self.consts.pixel_number // 2)
        return zernike_value_array


class ZernikeToSurface(Surface):

    def __init__(self, constants, zernike_number_list, zernike_value_array):
        """
        class : 2d surface from zernike values

        Parameters
        ----------
        constants : object
            instance by class : Constants
        zernike_number_list : 1d-list of int
            zernike number which to use (1 means piston)
        zernike_value_array : 1d-array of float
            [m] value of zernike coefficient coresponding to zernike_number_list
        offset_height_percent : float
            ignore height in percent. if you set 2, the lower 2% is ignored in self._volume_calculation()

        Returns
        -------
        None.

        """

        self.consts = constants
        self.zernike_number_list = zernike_number_list
        self.zernike_value_array = zernike_value_array

        self.surface = super()._make_masked_zernike_surface(self.zernike_number_list,
                                                            self.zernike_value_array)
        self.pv = super()._pv_calculation()
        self.rms = super()._rms_calculation()
        self.volume = super()._volume_calculation()[0]
        self.offset_height_value = super()._volume_calculation()[1]

    def h(self):
        mkhelp(self)


class StitchedCsvToSurface(Surface):
    def __init__(
            self,
            constants,
            original_stitched_csv_fpath,
            deformed_stitched_csv_fpath):
        """
        class : 2d-surface from measured (stitched) 2d-csv data

        Parameters
        ----------
        constants : TYPE
            instance by class : Constants
        original_stitched_csv_fpath : str
            変形前の測定データ[mm]のcsvパス
        deformed_stitched_csv_fpath : TYPE
            変形後の測定データ[mm]のcsvパス
            "" ならoriginalだけで計算
        offset_height_percent : float
            ignore height in percent. if you set 2, the lower 2% is ignored in self._volume_calculation()

        Returns
        -------
        None.

        """
        self.consts = constants

        if deformed_stitched_csv_fpath == "":
            self.surface = self.__read_csv_to_masked_surface(original_stitched_csv_fpath)

        else:
            self.original_masked_surface = self.__read_csv_to_masked_surface(
                original_stitched_csv_fpath)
            self.deformed_masked_surface = self.__read_csv_to_masked_surface(
                deformed_stitched_csv_fpath)
            self.surface = self.deformed_masked_surface - self.original_masked_surface

        self.pv = super()._pv_calculation()
        self.rms = super()._rms_calculation()
        self.volume = super()._volume_calculation()[0]
        self.offset_height_value = super()._volume_calculation()[1]

    def __read_csv_to_masked_surface(self, filepath):
        raw = np.loadtxt(filepath)
        raw_zero_fill = np.where(np.isnan(raw), 0, raw)
        image = PIL.Image.fromarray(raw_zero_fill)
        img_resize = image.resize(
            size=(
                self.consts.pixel_number,
                self.consts.pixel_number))

        masked_surface = self.consts.mask * np.array(img_resize) * 1e-3
        return masked_surface


class FilteredSurface(Surface):
    def __init__(self, constants, inputed_surface, filter_parameter):
        self.consts = constants
        self.inputed_surface = inputed_surface
        self.filter_parameter = filter_parameter

        self.surface = self.__smoothing_filter()
        self.pv = super()._pv_calculation()
        self.rms = super()._rms_calculation()
        self.volume = super()._volume_calculation()[0]
        self.offset_height_value = super()._volume_calculation()[1]

    def h(self):
        mkhelp(self)

    def __smoothing_filter(self):
        surface_without_nan = np.where(self.consts.tf,
                                       self.inputed_surface,
                                       0)
        filtered_surface = sp.ndimage.filters.uniform_filter(surface_without_nan,
                                                             size=self.filter_parameter)

        masked_filtered_surface = self.consts.mask * filtered_surface
        return masked_filtered_surface


class OapSurface(Surface):
    def __init__(self, constants, radius_of_curvature, off_axis_distance, clocking_angle_rad):
        self.consts = constants
        self.radius_of_curvature = radius_of_curvature
        self.off_axis_distance = off_axis_distance
        self.clocking_angle_rad = clocking_angle_rad

        oap = oap_calculation(clocking_angle_rad=self.clocking_angle_rad,
                              radius_of_curvature=self.radius_of_curvature,
                              off_axis_distance=self.off_axis_distance,
                              x_mesh=self.consts.xx,
                              y_mesh=self.consts.yy)
        self.surface = self.consts.mask * oap

        self.pv = super()._pv_calculation()
        self.rms = super()._rms_calculation()
        self.volume = super()._volume_calculation()[0]
        self.offset_height_value = super()._volume_calculation()[1]

    def h(self):
        mkhelp(self)


class ZernikeToTorque:
    def __init__(
            self,
            constants,
            target_zernike_number_list,
            target_zernike_value_array,
            ignore_zernike_number_list):
        """
        class : torque values in order to reproduct zernike values

        Parameters
        ----------
        constants : TYPE
            instance by class : Constants
        target_zernike_number_list : 1d-list of int
            zernike number which to use (1 means piston)
        target_zernike_value_array : 1d-array of float
            [m] value of zernike coefficient coresponding to target_zernike_number_list
        ignore_zernike_number_list : 1d-list of int
            zernike number which is not used in WH reproduction

        Returns
        -------
        None.

        """

        self.consts = constants
        self.target_zernike_number_list = target_zernike_number_list
        self.target_zernike_value_array = target_zernike_value_array
        self.ignore_zernike_number_list = ignore_zernike_number_list

        self.full_zernike_value_array = self.__make_full_zernike_value_array()
        self.remaining_operation_matrix = make_remaining_matrix(self.consts.operation_matrix, self.ignore_zernike_number_list)
        self.remaining_zernike_value_array = make_remaining_matrix(self.full_zernike_value_array, self.ignore_zernike_number_list)
        self.remaining_zernike_number_list = make_remaining_matrix(1 + np.arange(self.consts.zernike_max_degree), self.ignore_zernike_number_list)

        self.torque_value_array = self.__make_torque_value_array()

    def h(self):
        mkhelp(self)

    def __make_full_zernike_value_array(self):
        target_zernike_number_idx_array = np.array(self.target_zernike_number_list) - 1

        full_zernike_value_array = np.zeros(self.consts.zernike_max_degree)
        for i in range(len(self.target_zernike_value_array)):
            target_zernike_number_idx = target_zernike_number_idx_array[i]
            full_zernike_value_array[target_zernike_number_idx] = self.target_zernike_value_array[i]
        return full_zernike_value_array

    def __make_torque_value_array(self):
        inverse_operation_matrix = sp.linalg.pinv(1e9 * self.remaining_operation_matrix) * 1e9

        torque_value_array = np.dot(inverse_operation_matrix, self.remaining_zernike_value_array)
        return torque_value_array


class TorqueToZernike:

    def __init__(self, constants, torque_value_array, restructed_torque_value, ignore_zernike_number_list):
        """
        class : calculate zernike values deformed by torque values

        Parameters
        ----------
        constants : object
            instance by class : Constants
        torque_value_array : 1d-array of float
            [mm] torque values (array size is 36)
        restructed_torque_value : float
            [mm] mechanical limit of WH
                attr with prefix "restructed_" is used the restructed_torque_value
        ignore_zernike_number_list : 1d-list of int
            zernike number which is ignored in calculating of opreation matrix and inversed matrix

        Returns
        -------
        None.

        """

        self.consts = constants
        self.torque_value_array = torque_value_array
        self.restructed_torque_value = abs(restructed_torque_value)
        self.ignore_zernike_number_list = ignore_zernike_number_list

        self.remaining_operation_matrix = make_remaining_matrix(self.consts.operation_matrix,
                                                                self.ignore_zernike_number_list)
        self.restructed_torque_value_array = self.__make_restructed_torque_value_array()

        self.remaining_reproducted_zernike_value_array = self.__make_reproducted_zernike_value_array(self.torque_value_array)
        self.remaining_reproducted_restructed_zernike_value_array = self.__make_reproducted_zernike_value_array(self.restructed_torque_value_array)
        self.remaining_zernike_number_list = make_remaining_matrix(1 + np.arange(self.consts.zernike_max_degree),
                                                                   self.ignore_zernike_number_list)

    def h(self):
        mkhelp(self)

    def __make_restructed_torque_value_array(self):
        only_max_restructed_torque_value_array = np.where(self.torque_value_array < self.restructed_torque_value,
                                                          self.torque_value_array,
                                                          self.restructed_torque_value)

        restructed_torque_value_array = np.where(only_max_restructed_torque_value_array > -self.restructed_torque_value,
                                                 only_max_restructed_torque_value_array,
                                                 -self.restructed_torque_value)

        return restructed_torque_value_array

    def __make_reproducted_zernike_value_array(self, using_torque_value_array):
        remaining_reproducted_zernike_value_array = np.dot(
            self.remaining_operation_matrix, using_torque_value_array)
        return remaining_reproducted_zernike_value_array

    def make_torque_plot(self, figure=plt.figure(), position=111):
        fontsize = 15
        ax_title = "WH (black : 1-12, green : 13-24, violet : 25-36)"

        x = np.arange(1, 13)
        x_str = []

        for i in range(12):
            text = ""
            for j in range(3):
                num = str(12 * j + i + 1).zfill(2)
                text = text + "\n" + num
            x_str.append(text)

        torque = self.torque_value_array
        ax = figure.add_subplot(position)
        ax.plot(x, torque[0:12], color="black", marker="s", linewidth=1)
        ax.plot(x, torque[12:24], color="green", marker="o", linewidth=1)
        ax.plot(x, torque[24:36], color="darkviolet", marker="^", linewidth=1)

        ax.set_title(ax_title, fontsize=fontsize)
        ax.grid()
        ax.set_xlabel("Motor Number", fontsize=fontsize)
        ax.set_xticks(x)
        ax.set_xticklabels(x_str)
        ax.set_ylabel("Motor drive \namount [mm]", fontsize=fontsize)

        ax.hlines([self.restructed_torque_value, -self.restructed_torque_value],
                  xmin=1, xmax=12, color="red", linestyle="dashed")
        ax.hlines(0, xmin=1, xmax=12, color="darkgray")

        return ax


class OapConstants:
    def __init__(self,
                 ideal_radius_of_curvature, ideal_off_axis_distance, ideal_clocking_angle_rad,
                 delta_radius_of_curvature, delta_off_axis_distance, delta_clocking_angle_rad):
        """


        Parameters
        ----------
        ideal_radius_of_curvature : float
            ideal
        ideal_off_axis_distance : float
            ideal
        ideal_clocking_angle_rad : float
            ideal
        delta_radius_of_curvature : float
            ideal+delta = init parameter for x0 in sp.optimize.minimize()
        delta_off_axis_distance : float
            ideal+delta = init parameter for x0 in sp.optimize.minimize()
        delta_clocking_angle_rad : float
            ideal+delta = init parameter for x0 in sp.optimize.minimize()

        Returns
        -------
        None.

        """

        self.ideal_radius_of_curvature = ideal_radius_of_curvature
        self.ideal_off_axis_distance = ideal_off_axis_distance
        self.ideal_clocking_angle_rad = ideal_clocking_angle_rad

        self.delta_radius_of_curvature = delta_radius_of_curvature
        self.delta_off_axis_distance = delta_off_axis_distance
        self.delta_clocking_angle_rad = delta_clocking_angle_rad

        self.minimize_init_list = self.__make_minimize_init_list()

    def h(self):
        mkhelp(self)

    def __make_minimize_init_list(self):
        minimize_init_list = [self.ideal_radius_of_curvature + self.delta_radius_of_curvature,
                              self.ideal_off_axis_distance + self.delta_off_axis_distance,
                              self.ideal_clocking_angle_rad + self.delta_clocking_angle_rad]

        return minimize_init_list


class OapMinimize:
    def __init__(self, constants, oap_constants, inputed_surface):
        self.consts = constants
        self.oap_consts = oap_constants
        self.inputed_surface = inputed_surface

        ideal_oap = OapSurface(constants=self.consts,
                               radius_of_curvature=self.oap_consts.ideal_radius_of_curvature,
                               off_axis_distance=self.oap_consts.ideal_off_axis_distance,
                               clocking_angle_rad=self.oap_consts.ideal_clocking_angle_rad)

        self.__ideal_oap_surface = ideal_oap.surface

        self.minimize_args_list = self.__make_minimize_args_list()

        print("OapMinimize start")
        start_time = time.time()

        minimize_result = sp.optimize.minimize(fun=self.__minimize_input_function,
                                               x0=self.oap_consts.minimize_init_list,
                                               args=(self.minimize_args_list),
                                               method="Powell")

        print("OapMinimize finished")
        end_time = time.time()

        self.result_all = minimize_result

        self.result_time = end_time - start_time
        self.result_parameters = minimize_result["x"]
        return

    def h(self):
        mkhelp(self)

    def __make_minimize_args_list(self):
        args_list = [self.consts,
                     self.inputed_surface,
                     self.__ideal_oap_surface]

        return args_list

    def __minimize_input_function(self, X, args_list):
        test_radius_of_curvature = X[0]
        test_off_axis_distance = X[1]
        test_clocking_angle_rad = X[2]

        arg_consts = args_list[0]
        arg_inputed_surface = args_list[1]
        arg_ideal_oap_surface = args_list[2]

        inputed_surface_physical_height = arg_inputed_surface + arg_ideal_oap_surface

        test_oap_obj = OapSurface(constants=arg_consts,
                                  radius_of_curvature=test_radius_of_curvature,
                                  off_axis_distance=test_off_axis_distance,
                                  clocking_angle_rad=test_clocking_angle_rad)

        test_oap_surface_physical_height = test_oap_obj.surface

        difference_of_surface = inputed_surface_physical_height - test_oap_surface_physical_height

        difference_surface_obj = Surface(constants=arg_consts,
                                         surface=difference_of_surface)

        volume = difference_surface_obj.volume

        del test_oap_obj
        del difference_surface_obj
        return volume


class CirclePathMeasurementReading:
    def __init__(self,
                 Constants,
                 original_csv_fpath: str,
                 deformed_csv_fpath: str
                 ) -> None:

        self.consts = Constants
        self.df_raw_original = pd.read_csv(original_csv_fpath)

        if deformed_csv_fpath == "":
            self.df_diff = pd.DataFrame({"degree": self.df_raw_original["angle"].values,
                                         "radian": np.deg2rad(self.df_raw_original["angle"].values),
                                         "height": self.df_raw_original["height"].values})

        else:
            self.df_raw_deformed = pd.read_csv(deformed_csv_fpath)
            height_diff = self.df_raw_deformed["height"].values - self.df_raw_original["height"].values

            self.df_diff = pd.DataFrame({"degree": self.df_raw_original["angle"].values,
                                         "radian": np.deg2rad(self.df_raw_original["angle"].values),
                                         "height": height_diff})

    def h(self):
        mkhelp(self)


class CirclePathZernikeFitting:
    def __init__(self,
                 Constants,
                 circle_path_radius: float,
                 df_diff: pd.DataFrame,
                 ignore_zernike_number_list: list[int]
                 ) -> None:
        """CirclePathZernikeFitting
        zernike多項式を除去

        Parameters
        ----------
        Constants : class

        circle_path_radius : float
            円環パスの半径 [m]
        df_diff : pd.DataFrame
            CirclePathMeasurementReading の df_diff
        ignore_zernike_number_list : list[int]
            zernikeでfittingした後に除去する項（中身は1~11まで）
        """

        self.consts = Constants
        self.circle_path_radius = circle_path_radius
        self.df_diff = df_diff
        self.ignore_zernike_number_list = ignore_zernike_number_list

        zernike_fitting_result = self.__zernike_fitting()
        self.optimize_result = zernike_fitting_result["optimize_result"]
        self.zernike_polynomial_array = zernike_fitting_result["zernike_polynomial_array"]

        zernike_removing_result = self.__zernike_removing(zernike_polynomial_array=self.zernike_polynomial_array)
        self.removing_zernike = zernike_removing_result["removing_zernike"]
        self.height_removed = zernike_removing_result["residual"]

    def h(self):
        mkhelp(self)

    def __zernike_fitting(self):

        def zernike_r_const_polynomial_calculation(coef_r_const: list[float],
                                                   radius: float,
                                                   theta: ndarray) -> ndarray:
            """zernike_r_const_polynomial_calculation
            radiusが固定値（＝円環パス）の場合のzernike 多項式の計算
            そのままのzernikeでフィッティングすると、
            thetaによる計算結果の部分が全く同じ値になる項（1, 4, 11など）が出て各係数が独立ではなくなってしまう

            Parameters
            ----------
            coef_r_const : list[float]
                radius固定値の時のzernike係数ベクトル
            radius : float
                円環パスの半径（固定値なのでarrayではない）
            theta : ndarray
                角度 [rad]

            Returns
            -------
            ndarray
                zernike 多項式の計算結果
            """

            zernike1_ = zernike_term_calculation(1, radius, theta)
            zernike2_ = zernike_term_calculation(2, radius, theta)
            zernike3_ = zernike_term_calculation(3, radius, theta)
            zernike4_ = zernike_term_calculation(4, radius, theta)
            zernike5_ = zernike_term_calculation(5, radius, theta)
            zernike6_ = zernike_term_calculation(6, radius, theta)
            zernike7_ = zernike_term_calculation(7, radius, theta)
            zernike8_ = zernike_term_calculation(8, radius, theta)
            zernike9_ = zernike_term_calculation(9, radius, theta)
            zernike10_ = zernike_term_calculation(10, radius, theta)
            zernike11_ = zernike_term_calculation(11, radius, theta)

            zernike1_4_11 = coef_r_const[0] * (zernike1_ + zernike4_ + zernike11_)
            zernike2_8 = coef_r_const[1] * (zernike2_ + zernike8_)
            zernike3_7 = coef_r_const[2] * (zernike3_ + zernike7_)
            zernike5 = coef_r_const[3] * zernike5_
            zernike6 = coef_r_const[4] * zernike6_
            zernike9 = coef_r_const[5] * zernike9_
            zernike10 = coef_r_const[6] * zernike10_

            zernike_r_const_polynomial = zernike1_4_11 + zernike2_8 + zernike3_7 + zernike5 + zernike6 + zernike9 + zernike10

            return zernike_r_const_polynomial

        def minimize_funciton_sq(x, params_):
            radius_, theta_array_, height_array_ = params_
            zernike_array_ = zernike_r_const_polynomial_calculation(coef_r_const=x,
                                                                    radius=radius_,
                                                                    theta=theta_array_)

            residual = height_array_ - zernike_array_
            return residual

        def make_zernike_polynomial_from_r_const_zernike(zernike_r_const_polynomial_list: list[float]
                                                         ) -> ndarray:
            """make_zernike_polynomial_from_r_const_zernike
            radius固定の場合のzernike係数ベクトルを通常のzernike係数ベクトルとして並べ直す

            Parameters
            ----------
            zernike_r_const_polynomial_list : list[float]
                radius固定の場合のzernike係数

            Returns
            -------
            ndarray
                通常のzernike係数ベクトル len(ndarray) = 11
            """

            zernike_polynomial_array_ = np.empty(11)
            zernike_polynomial_array_[0] = zernike_r_const_polynomial_list[0]
            zernike_polynomial_array_[1] = zernike_r_const_polynomial_list[1]
            zernike_polynomial_array_[2] = zernike_r_const_polynomial_list[2]
            zernike_polynomial_array_[3] = zernike_r_const_polynomial_list[0]
            zernike_polynomial_array_[4] = zernike_r_const_polynomial_list[3]
            zernike_polynomial_array_[5] = zernike_r_const_polynomial_list[4]
            zernike_polynomial_array_[6] = zernike_r_const_polynomial_list[2]
            zernike_polynomial_array_[7] = zernike_r_const_polynomial_list[1]
            zernike_polynomial_array_[8] = zernike_r_const_polynomial_list[5]
            zernike_polynomial_array_[9] = zernike_r_const_polynomial_list[6]
            zernike_polynomial_array_[10] = zernike_r_const_polynomial_list[0]

            return zernike_polynomial_array_

        radius = self.circle_path_radius
        theta_array = self.df_diff["radian"].values
        height_array = self.df_diff["height"].values

        params = [radius, theta_array, height_array]

        optimize_result_sq = optimize.least_squares(fun=minimize_funciton_sq,
                                                    x0=(np.ones(7) * 1e-9),
                                                    args=(params, ))

        zernike_polynomial_array = make_zernike_polynomial_from_r_const_zernike(optimize_result_sq["x"])

        result_dict = {"optimize_result": optimize_result_sq,
                       "zernike_polynomial_array": zernike_polynomial_array}

        return result_dict

    def __zernike_removing(self,
                           zernike_polynomial_array: ndarray
                           ) -> dict:
        """__zernike_removing
        zernike係数ベクトルを用いてzernike成分を除去

        Parameters
        ----------
        zernike_polynomial_array : ndarray
            通常のzernike係数ベクトル

        Returns
        -------
        dict
            removing_zernike: 除去するzernike成分の計算結果のarray
            residual: height - removing_zernike
        """
        def make_remaining_array(ignore_list: list[int],
                                 polynomial_array: ndarray) -> ndarray:
            tf_array = np.zeros(11)
            for num in ignore_list:
                tf_array[num - 1] = 1
            remaining_polynomial_array = tf_array * polynomial_array
            return remaining_polynomial_array

        radian_array = self.df_diff["radian"].values
        height_array = self.df_diff["height"].values
        ignore_zernike_number_list = self.ignore_zernike_number_list

        remaining_zernike_polynomial_array = make_remaining_array(ignore_list=ignore_zernike_number_list,
                                                                  polynomial_array=zernike_polynomial_array)

        removing_zernike_array = zernike_polynomial_calculation(coef=remaining_zernike_polynomial_array,
                                                                radius=self.circle_path_radius,
                                                                theta=radian_array)

        residual_height_array = height_array - removing_zernike_array

        result_dict = {"removing_zernike": removing_zernike_array,
                       "residual": residual_height_array}

        return result_dict
