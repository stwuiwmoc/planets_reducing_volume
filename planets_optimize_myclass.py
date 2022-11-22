# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 14:59:35 2021

@author: swimc
"""
import time

import cv2
import matplotlib
import mpl_toolkits.axes_grid1
import numpy as np
import pandas as pd
import PIL
import proper as pr
import scipy as sp
from matplotlib import cm
from matplotlib.colors import Normalize
from numpy import ndarray
from scipy import interpolate, optimize


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


def get_latest_commit_datetime() -> list[str]:
    """現在のツリーで直近のコミット日時を文字列として取得する

    Returns
    -------
    list[str]
        ["Latest commit datetime", コミットした日付, コミットした時刻]
    """

    import subprocess

    git_command = ["git", "log", "-1", "--format='%cI'"]

    proc = subprocess.run(
        git_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)

    latest_commit_date_text = proc.stdout.decode("utf-8")

    # タイムゾーンを表す "+0900" を削除
    latest_commit_datetime_text_without_timezone = latest_commit_date_text[1:-8]

    # 日付のみを含む文字列
    latest_commit_date_text = latest_commit_datetime_text_without_timezone[:10]

    # 時刻のみを含む文字列
    latest_commit_time_text = latest_commit_datetime_text_without_timezone[11:]

    return "Latest commit datetime", latest_commit_date_text, latest_commit_time_text


def have_some_change_in_git_status() -> bool:
    """git的な意味でのファイルの変更の有無をboolianで取得する

    Returns
    -------
    bool
        ファイルの変更箇所ありならTrue, 無しならFalse
    """

    import subprocess

    git_command = ["git", "status", "--short"]

    proc = subprocess.run(
        git_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)

    proc_text = proc.stdout.decode("utf-8")

    if len(proc_text) == 0:
        # git status --short で出力される文字列が無い
        # つまり「ファイルの変更箇所なし」を意味する
        have_some_change = False

    else:
        # git status --short で出力された文字列がある
        # つまり「ファイルに何かしらの変更箇所あり」を意味する
        have_some_change = True

    return have_some_change


def plot_parameter_table(
        fig: matplotlib.figure.Figure,
        position: matplotlib.gridspec.GridSpec,
        parameter_table: list,
        fontsize: int) -> matplotlib.axes._subplots.Axes:
    """パラメータ表示用のtableをax内に作成

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figureオブジェクト
    position: matplotlib.gridspec.GridSpec
        fig内での配置
    parameter_table : list
        パラメータリスト。縦横の要素数がそれぞれ一致している必要がある。（正方行列のように縦横で同じ要素数になる必要はない）
    fontsize : int
        テーブル内の文字サイズ

    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot
        Axesオブジェクト
    """

    ax = fig.add_subplot(position)

    table = ax.table(
        cellText=parameter_table,
        loc="center")

    # fontsizeの調整
    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)

    # 軸ラベルの消去
    ax.axis("off")

    # 縦幅をaxオブジェクトの縦幅に合わせて調整
    for pos, cell in table.get_celld().items():
        cell.set_height(1 / len(parameter_table))

    return ax


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


def zernike_term_calculation(
        term_number: int,
        pupil_radius: float,
        normalization_radius: ndarray,
        theta: ndarray) -> ndarray:
    """zernike_term_calculation
    zernikeの各項を計算する（係数は入ってない）

    Parameters
    ----------
    term_number : int
        [無次元] 計算したいzernikeの項の番号（1～11）
    pupil_radius : float
        [m] 瞳の半径
    radius : ndarray
        [m] 半径の1次元array
    theta : ndarray
        [rad] 角度の1次元array

    Returns
    -------
    ndarray
        zernikeのどれか一つの項の計算結果（係数は入ってない）
    """

    # zernike多項式は単位円上での話なので、瞳半径で1になるように正規化する
    normalization_radius = normalization_radius / pupil_radius

    if term_number == 1:
        return 1
    if term_number == 2:
        return 2 * normalization_radius * np.cos(theta)
    if term_number == 3:
        return 2 * normalization_radius * np.sin(theta)
    if term_number == 4:
        return np.sqrt(3) * (2.0 * normalization_radius ** 2 - 1.0)
    if term_number == 5:
        return np.sqrt(6) * (normalization_radius ** 2) * np.sin(2 * theta)
    if term_number == 6:
        return np.sqrt(6) * (normalization_radius ** 2) * np.cos(2 * theta)
    if term_number == 7:
        return np.sqrt(8) * (3.0 * normalization_radius ** 3 - 2.0 * normalization_radius) * np.sin(theta)
    if term_number == 8:
        return np.sqrt(8) * (3.0 * normalization_radius ** 3 - 2.0 * normalization_radius) * np.cos(theta)
    if term_number == 9:
        return np.sqrt(8) * (normalization_radius ** 3) * np.sin(3 * theta)
    if term_number == 10:
        return np.sqrt(8) * (normalization_radius ** 3) * np.cos(3 * theta)
    if term_number == 11:
        return np.sqrt(5) * (6.0 * normalization_radius ** 4 - 6.0 * normalization_radius ** 2 + 1.0)
    else:
        print("term_number must be 1~11")
        return


def zernike_polynomial_calculation(
        coef: list[float],
        pupil_radius: float,
        radius: ndarray,
        theta: ndarray) -> ndarray:
    """
    zernike多項式のarrayを計算する（係数込み）

    Parameters
    ----------
    coef : list[float]
        zernike係数ベクトル
    pupil_radius : float
        [m] 瞳の半径
    radius : ndarray
        [m] 半径の1次元array
    theta : ndarray
        [rad] 角度の1次元array

    Returns
    -------
    ndarray
        zernike多項式のarray
    """

    zernike1 = coef[0] * zernike_term_calculation(1, pupil_radius, radius, theta)
    zernike2 = coef[1] * zernike_term_calculation(2, pupil_radius, radius, theta)
    zernike3 = coef[2] * zernike_term_calculation(3, pupil_radius, radius, theta)
    zernike4 = coef[3] * zernike_term_calculation(4, pupil_radius, radius, theta)
    zernike5 = coef[4] * zernike_term_calculation(5, pupil_radius, radius, theta)
    zernike6 = coef[5] * zernike_term_calculation(6, pupil_radius, radius, theta)
    zernike7 = coef[6] * zernike_term_calculation(7, pupil_radius, radius, theta)
    zernike8 = coef[7] * zernike_term_calculation(8, pupil_radius, radius, theta)
    zernike9 = coef[8] * zernike_term_calculation(9, pupil_radius, radius, theta)
    zernike10 = coef[9] * zernike_term_calculation(10, pupil_radius, radius, theta)
    zernike11 = coef[10] * zernike_term_calculation(11, pupil_radius, radius, theta)

    zernike_polynomial = zernike1 + zernike2 + zernike3 + zernike4 + zernike5 + zernike6 + zernike7 + zernike8 + zernike9 + zernike10 + zernike11
    return zernike_polynomial


class Constants:

    def __init__(
            self, physical_radius, ignore_radius,
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
        self.xx, self.yy = make_meshgrid(
            -physical_radius, physical_radius,
            -physical_radius, physical_radius,
            pixel_number)
        self.tf = np.where(self.xx ** 2 + self.yy ** 2 <= self.varid_radius**2, True, False)
        self.mask = np.where(self.tf, 1, np.nan)
        self.zernike_max_degree = zernike_max_degree
        operation_matrix_fpath = "mkfolder/make_opration_matrix/WT06_zer" + str(self.zernike_max_degree) + "_opration_matrix[m].csv"
        self.operation_matrix = np.genfromtxt(operation_matrix_fpath, delimiter=",").T

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
        self.zernike_value_array = self._zernike_value_array_calculation(self.surface)

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

        lower_ignored_surface = np.where(
            self.surface >= offset_height_value,
            self.surface - offset_height_value,
            0)

        # 1pixelあたりの単位面積を計算
        physical_diameter = 2 * self.consts.physical_radius
        unit_pixel_area = (physical_diameter / self.consts.pixel_number)**2

        # 1 pixelあたりの体積を計算
        unit_pixel_volume = unit_pixel_area * lower_ignored_surface
        volume_in_m3 = unit_pixel_volume.sum()
        return (volume_in_m3, offset_height_value)

    def _zernike_value_array_calculation(self, surface):
        surface_without_nan = np.where(self.consts.tf, surface, 0)

        zernike_value_array = pr.prop_fit_zernikes(
            wavefront0=surface_without_nan,
            pupil0=self.consts.tf,
            pupilradius0=self.consts.pixel_number // 2,
            nzer=self.consts.zernike_max_degree,
            xc=self.consts.pixel_number // 2,
            yc=self.consts.pixel_number // 2)
        return zernike_value_array

    def _make_masked_zernike_surface(self, full_zernike_value_array: ndarray) -> ndarray:
        """_make_masked_zernike_surface
            決まった半径でnp.nanでmaskされたzernike平面を作る

        Parameters
        ----------
        full_zernike_value_array : ndarray
            zernike係数ベクトル len() = Constants.zernike_max_degree

        Returns
        -------
        ndarray
            zernike係数ベクトルによって作られる2次元平面
        """

        optical_wavelength = 500e-9

        full_zernike_number_list = [i + 1 for i in range(self.consts.zernike_max_degree)]

        wavestruct = pr.prop_begin(
            beam_diameter=2 * self.consts.varid_radius,
            lamda=optical_wavelength,
            grid_n=self.consts.pixel_number,
            beam_diam_fraction=1)

        wfe = pr.prop_zernikes(
            wavestruct,
            full_zernike_number_list,
            full_zernike_value_array)

        masked_wfe = self.consts.mask * wfe
        return masked_wfe

    def calc_circle_path_height(
            self,
            radius: float,
            angle_division_number: int) -> list[np.ndarray, np.ndarray]:
        """与えられた半径に対応した円環パスに沿った高さを計算する

        Parameters
        ----------
        radius : float
            [m] 円環パスの半径
        angle_division_number : int
            [無次元] 0deg～359degを何分割するか

        Returns
        -------
        list[np.ndarray, np.ndarray]
            list[0] : [deg] 円環パスの角度（CCW, x軸正方向を0degとする）\n
            list[1] : [m] 半径に対応した円環パスの高さ
        """
        varid_radius_pixel_number = int(self.consts.varid_radius / self.consts.physical_radius * self.consts.pixel_number / 2)
        measurement_radius_idx = int(radius * 1e3)

        image = np.where(self.consts.tf, self.surface, 0)
        flags = cv2.INTER_CUBIC + cv2.WARP_FILL_OUTLIERS + cv2.WARP_POLAR_LINEAR

        linear_polar_image = cv2.warpPolar(
            src=image,
            dsize=(int(self.consts.varid_radius * 1e3), angle_division_number),
            center=(self.consts.pixel_number / 2, self.consts.pixel_number / 2),
            maxRadius=varid_radius_pixel_number,
            flags=flags)

        # 与えられた半径に対応する高さだけを取り出し
        circle_path_line = linear_polar_image[:, measurement_radius_idx]

        # 0degから一周するまでの角度の1次元arrayを作成
        angle_array = np.linspace(0, 359, angle_division_number)

        return angle_array, circle_path_line

    def make_image_plot(
            self, figure,
            position=111,
            cbar_surface=None,
            cbar_min_percent=0, cbar_max_percent=100,
            pv_digits=2, rms_digits=2):

        cmap = cm.jet
        fontsize = 15
        title = "pv = " + str(round(self.pv * 1e6, pv_digits)) + " [um]" + "\n" + "RMS = " + str(round(self.rms * 1e6, rms_digits)) + " [um]"

        if cbar_surface is None:
            cbar_surface = self.surface
        else:
            pass

        cbar_pv = np.nanmax(cbar_surface) - np.nanmin(cbar_surface)
        cbar_min = np.nanmin(cbar_surface) + cbar_pv * cbar_min_percent / 100
        cbar_max = np.nanmin(cbar_surface) + cbar_pv * cbar_max_percent / 100

        extent = [
            -self.consts.physical_radius, self.consts.physical_radius,
            -self.consts.physical_radius, self.consts.physical_radius]

        ax = figure.add_subplot(position)
        ax.imshow(
            self.surface, interpolation="nearest", cmap=cmap,
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

    def make_3d_plot(
            self, figure, position=111,
            cbar_surface=None,
            cbar_min_percent=0, cbar_max_percent=100,
            zaxis_bottom_percent=30):

        fontsize = 10

        xx_ = self.consts.xx
        yy_ = self.consts.yy

        zz_bottom_value = np.nanmin(self.surface) + self.pv * zaxis_bottom_percent / 100

        zz_ = np.where(
            self.consts.tf,
            self.surface,
            zz_bottom_value)

        if cbar_surface is None:
            cbar_surface = self.surface
        else:
            pass

        cbar_pv = np.nanmax(cbar_surface) - np.nanmin(cbar_surface)
        cbar_min = np.nanmin(cbar_surface) + cbar_pv * cbar_min_percent / 100
        cbar_max = np.nanmin(cbar_surface) + cbar_pv * cbar_max_percent / 100

        ax = figure.add_subplot(position, projection="3d")
        ax.plot_surface(xx_, yy_, zz_, cmap=cm.jet, vmin=cbar_min, vmax=cbar_max)

        ax.set_zlim(np.nanmin(cbar_surface), np.nanmax(cbar_surface))
        ax.set_xlabel("x [m]", fontsize=fontsize)
        ax.set_ylabel("y [m]", fontsize=fontsize)
        ax.set_zticklabels([])

        ax.view_init(elev=30, azim=-60)

        norm = Normalize(vmin=cbar_min, vmax=cbar_max)
        mappable = cm.ScalarMappable(norm=norm, cmap=cm.jet)
        cbar = figure.colorbar(mappable, ax=ax, shrink=0.8)
        cbar_title = "[m]"
        cbar.set_label(cbar_title, fontsize=fontsize)

        return ax

    def make_circle_path_plot(
            self,
            figure: matplotlib.figure.Figure,
            position: matplotlib.gridspec.GridSpec = 111,
            radius: float = 0.870,
            angle_division_number: int = 360,
            height_magn: float = 1e9,
            height_unit_str: str = "[nm]",
            line_label: str = "") -> matplotlib.axes._subplots.Axes:
        """与えられた半径に対する円環パスの高さをプロット

        Parameters
        ----------
        figure : matplotlib.figure.Figure
            Figureインスタンス
        position : matplotlib.gridspec.GridSpec, optional
            GridSpecインスタンス, by default 111
        radius : float, optional
            [m] 円環パスの半径, by default 0.870
        angle_division_number : int, optional
            [無次元] 円環パスの角度0deg～359degを何分割するか, by default 360
        height_magn : float, optional
            [無次元] 高さ方向の倍率、デフォルトでは [m] 単位から [nm] 単位に換算してから出力される, by default 1e9
        height_unit_str : str, optional
            高さの軸ラベルで表記する単位, by default "[nm]"
        line_label : str, optional
            凡例のテキスト, by default ""

        Returns
        -------
        matplotlib.axes._subplots.Axes
            axインスタンス
        """

        fontsize = 15

        angle_array, height_array = self.calc_circle_path_height(
            radius=radius,
            angle_division_number=angle_division_number
        )

        circle_path_line = height_magn * height_array

        ax = figure.add_subplot(position)
        ax.plot(angle_array, circle_path_line, label=line_label)
        ax.grid()

        ax.set_xlabel("degree", fontsize=fontsize)
        ax.set_ylabel(
            "heignt " + height_unit_str + "\nat R=" + str(radius) + "m",
            fontsize=fontsize)

        return ax

    def make_torque_fulcrum_plot(
            self,
            ax,
            torque_value_array: np.ndarray,
            init_theta=150):

        def polar_coordinate_conversion(radius_: float, theta_, init_theta_):
            x_coordinate = radius_ * np.cos(np.deg2rad(theta_ + init_theta_))
            y_coordinate = radius_ * np.sin(np.deg2rad(theta_ + init_theta_))
            return x_coordinate, y_coordinate

        def decide_points_symbol(torque_value_array_, polar_coordinate_):
            plus_x, plus_y = np.where(
                torque_value_array_ > 0,
                polar_coordinate_,
                np.nan)
            minus_x, minus_y = np.where(torque_value_array_ < 0,
                                        polar_coordinate_,
                                        np.nan)

            return (plus_x, plus_y), (minus_x, minus_y)

        df = pd.read_csv("raw_data/torque_fulcrum_coordinate.csv", index_col="idx")

        polar_coordinate = polar_coordinate_conversion(
            radius_=df["radius"],
            theta_=df["theta"],
            init_theta_=init_theta)

        plus, minus = decide_points_symbol(
            torque_value_array,
            polar_coordinate)

        ax.scatter(plus[0], plus[1], s=120, c="black")
        ax.scatter(plus[0], plus[1], s=50, c="red")

        ax.scatter(minus[0], minus[1], s=120, c="white")
        ax.scatter(minus[0], minus[1], s=50, c="blue")

        return ax

    def make_zernike_value_plot(
            self,
            figure,
            position,
            label_=None):

        xaxis = np.arange(self.consts.zernike_max_degree) + 1
        xaxis_min = 0.5
        xaxis_max = self.consts.zernike_max_degree + 0.5

        ax = figure.add_subplot(position)

        ax.plot(xaxis, self.zernike_value_array, marker="s", label=label_)
        ax.hlines([0], xaxis_min, xaxis_max, color="gray")
        ax.grid()

        ax.set_xticks(xaxis)
        ax.set_xlim(xaxis_min, xaxis_max)

        ax.set_xlabel("zernike terms")
        ax.set_ylabel("zernike values [m]")

        if label_ is not None:
            ax.legend()

        return ax


class ZernikeRemovedSurface(Surface):
    def __init__(
            self, constants,
            inputed_surface: ndarray,
            removing_zernike_number_list: list[int]):
        """ZernikeRemovedSurface
        input_surfaceからremoving_zernike_number_listで指定したzernike項を除去してsurfaceとして返す

        Parameters
        ----------
        constants : [type]
            Constantsクラス
        inputed_surface : ndarray
            zernikeを除去したいsurface
        removing_zernike_number_list : list[int]
            除去したいzernikeの項数
        """
        self.consts = constants
        self.removing_zernike_number_list = removing_zernike_number_list

        inputed_zernike_value_array = super()._zernike_value_array_calculation(inputed_surface)

        removing_zernike_value_array = self.__make_full_length_zernike(
            original_value_array=inputed_zernike_value_array,
            removing_number_list=self.removing_zernike_number_list)

        removing_surface = super()._make_masked_zernike_surface(removing_zernike_value_array)
        self.zernike_value_array = inputed_zernike_value_array - removing_zernike_value_array
        self.surface = inputed_surface - removing_surface

        self.pv = super()._pv_calculation()
        self.rms = super()._rms_calculation()
        self.volume = super()._volume_calculation()[0]
        self.offset_height_value = super()._volume_calculation()[1]

    def h(self):
        mkhelp(self)

    def __make_full_length_zernike(self, original_value_array, removing_number_list):
        removing_value_array = np.zeros(self.consts.zernike_max_degree)
        for i in range(len(removing_number_list)):
            idx = removing_number_list[i] - 1
            removing_value_array[idx] = original_value_array[idx]
        return removing_value_array


class ZernikeToSurface(Surface):

    def __init__(
            self, constants,
            zernike_value_array: ndarray):
        """ZernikeToSurface
        zernike係数ベクトルから平面を計算する

        Parameters
        ----------
        constants : [type]
            Constant クラス
        zernike_value_array : ndarray
            zernike係数ベクトル len() = Constants.zernike_max_degree
        """

        self.consts = constants
        self.zernike_value_array = zernike_value_array

        self.surface = super()._make_masked_zernike_surface(self.zernike_value_array)
        self.pv = super()._pv_calculation()
        self.rms = super()._rms_calculation()
        self.volume = super()._volume_calculation()[0]
        self.offset_height_value = super()._volume_calculation()[1]

    def h(self):
        mkhelp(self)


class KagiStitchToSurface(Surface):
    def __init__(
            self,
            constants,
            txt_fpath: str):
        """__init__
        鍵谷先生のステッチ出力からSurfaceを作成

        Parameters
        ----------
        constants :
            Constantsクラス
        txt_fpath : str
            ステッチ出力のファイルパス
        """

        self.consts = constants
        self.stitched_txt_filepath = txt_fpath
        self.surface = self.consts.mask * self.__read_txt_and_interpolation()
        self.pv = super()._pv_calculation()
        self.rms = super()._rms_calculation()
        self.volume = super()._volume_calculation()[0]
        self.offset_height_value = super()._volume_calculation()[1]
        self.zernike_value_array = super()._zernike_value_array_calculation(self.surface)

    def h(self):
        mkhelp(self)

    def __read_txt_and_interpolation(self):
        """__read_txt_and_interpolation
        txtを読み込んで、メッシュに補間し、単位をmeterに換算
        ステッチファイル出力ははxy座標は[mm]、z座標は[nm]
        """
        def stitched_data_interpolation(
                x_old_array_mm: ndarray,
                y_old_array_mm: ndarray,
                z_old_array_nm: ndarray,
                x_new_mesh_mm: ndarray,
                y_new_mesh_mm: ndarray) -> ndarray:
            """stitched_data_interpolation
            データをメッシュに補間

            Parameters
            ----------
            x_old_array_mm : ndarray
                1次元array [mm]
            y_old_array_mm : ndarray
                1次元array [mm]
            z_old_array_nm : ndarray
                1次元array [nm] ※単位注意
            x_new_mesh_mm : ndarray
                2次元array [mm]
            y_new_mesh_mm : ndarray
                2次元array [mm]

            Returns
            -------
            ndarray
                z座標値の2次元array [nm]
            """

            xy_old_mm = np.stack([x_old_array_mm, y_old_array_mm], axis=1)
            z_new_mesh_nm = interpolate.griddata(
                points=xy_old_mm,
                values=z_old_array_nm,
                xi=(x_new_mesh_mm, y_new_mesh_mm),
                method="cubic",
                fill_value=0)

            return z_new_mesh_nm

        filepath = self.stitched_txt_filepath

        raw = np.loadtxt(filepath)
        x_raw_array = raw[:, 1]
        y_raw_array = raw[:, 2]
        z_raw_array = raw[:, 3]

        z_mesh_nm = stitched_data_interpolation(
            x_old_array_mm=x_raw_array,
            y_old_array_mm=y_raw_array,
            z_old_array_nm=z_raw_array,
            x_new_mesh_mm=self.consts.xx * 1e3,
            y_new_mesh_mm=self.consts.yy * 1e3)

        z_mesh_meter = z_mesh_nm * 1e-9
        return z_mesh_meter


class ExelisCsvToSurface(Surface):
    def __init__(
            self,
            constants,
            csv_filepath: str = "raw_data/digitFig01.csv"):

        self.consts = constants
        self.filepath = csv_filepath
        reshaped_surface = self.__raw_data_reshape()
        self.surface = 1e-6 * self.consts.mask * self.__data_resize(reshaped_surface)

        self.pv = super()._pv_calculation()
        self.rms = super()._rms_calculation()
        self.volume = super()._volume_calculation()[0]
        self.offset_height_value = super()._volume_calculation()[1]
        self.zernike_value_array = super()._zernike_value_array_calculation(self.surface)

    def h(self):
        mkhelp(self)

    def __raw_data_reshape(self):
        filepath = self.filepath
        raw = np.genfromtxt(
            fname=filepath,
            delimiter=",",
            encoding="utf-8_sig")

        reshaped = raw.reshape((1024, 1026))

        # 縦横の5.52182だけが並ぶ行と列を削除して1023*1023に成形
        deleted_temp = np.delete(reshaped, 0, 0)
        deleted = np.delete(deleted_temp, [0, 1024, 1025], 1)

        return deleted

    def __data_resize(self, surface: ndarray) -> ndarray:
        """__data_resize
        Exelisの1023 ** 2をpixel_number**2にリサイズ

        Parameters
        ----------
        surface : ndarray
            元の2次元array

        Returns
        -------
        ndarray
            size()=pixel_number**2 の2次元array
        """
        inner_nan_removed = np.where(
            np.isnan(surface),
            np.nanmean(surface),
            surface)

        image = PIL.Image.fromarray(inner_nan_removed)
        img_resize = image.resize(
            size=(
                self.consts.pixel_number,
                self.consts.pixel_number))

        return img_resize


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
        self.zernike_value_array = super()._zernike_value_array_calculation(self.surface)

    def h(self):
        mkhelp(self)

    def __smoothing_filter(self):
        surface_without_nan = np.where(
            self.consts.tf,
            self.inputed_surface,
            0)
        filtered_surface = sp.ndimage.filters.uniform_filter(
            surface_without_nan,
            size=self.filter_parameter)

        masked_filtered_surface = self.consts.mask * filtered_surface
        return masked_filtered_surface


class FemTxtToSurface(Surface):

    def __init__(
            self,
            constants,
            wt_version: int,
            wt_number: int):
        """__init__
        鍵谷先生のFEMモデル計算結果のtxtファイルからsurfaceを作る
        surfaceの値は (Fxx - F00) * (Constants.act_tuning_array のwt_numberに対応した値)

        Parameters
        ----------
        constants : Constants
            Constantsクラス
        wt_version : int
            鍵谷先生の出力ファイルのWTのバージョン
        wt_number : int
            計算したWTのact番号 1~36
        """

        self.consts = constants
        self.wt_number = wt_number

        self.fpath_original = "raw_data/Fxx/PM3.5_36ptAxWT" + str(wt_version).zfill(2) + "_F00.smesh.txt"
        self.fpath_deformed = "raw_data/Fxx/PM3.5_36ptAxWT" + str(wt_version).zfill(2) + "_F" + str(wt_number).zfill(2) + ".smesh.txt"
        df_original = self.__read_file(self.fpath_original)
        df_deformed = self.__read_file(self.fpath_deformed)

        self.act_tuning_array = self.__make_act_tuning_array()
        self.surface = self.__interpolation(df00=df_original, dfxx=df_deformed)

        self.pv = super()._pv_calculation()
        self.rms = super()._rms_calculation()
        self.volume = super()._volume_calculation()[0]
        self.offset_height_value = super()._volume_calculation()[1]
        self.zernike_value_array = super()._zernike_value_array_calculation(self.surface)

    def h(self):
        mkhelp(self)

    def __read_file(self, fpath: str) -> pd.DataFrame:
        """__read_file
        鍵谷先生のFEM出力データを読み出してdataframeとして出力

        Parameters
        ----------
        fpath : str
            読み込むfilepath

        Returns
        -------
        pd.DataFrame
            主に使うのはx, y, dz
        """

        df = pd.read_csv(
            fpath,
            sep="\\s+",
            skiprows=7,
            header=None)

        df.columns = ["x", "y", "z", "color", "dx", "dy", "dz"]
        return df

    def __make_act_tuning_array(self):
        act_tuning = np.array([
            37., 37., 20., 20., -38., 38.,
            37., 37., 20., 20., -38., 38.,
            37., 37., 20., 20., -38., 38.,
            37., 37., 20., 20., -38., 38.,
            37., 37., 20., 20., -38., 38.,
            37., 37., 20., 20., -38., 38.])

        return act_tuning

    def __interpolation(
            self,
            df00: pd.DataFrame,
            dfxx: pd.DataFrame) -> ndarray:
        """__interpolation
        FEM出力の点群をメッシュとして補間

        Parameters
        ----------
        df00 : pd.DataFrame
            変形前のFEM出力のDataFrame
        dfxx : pd.DataFrame
            変形後のFEM出力のDataFrame

        Returns
        -------
        ndarray
            変形前後でのz方向変位の2次元array (size = pixel_number * pixel_number)
        """

        x_1d = df00["x"]
        y_1d = df00["y"]
        dw_1d = dfxx["dz"] - df00["dz"]

        xy_old = np.stack([x_1d, y_1d], axis=1)

        dw_mesh = interpolate.griddata(
            points=xy_old,
            values=dw_1d,
            xi=(self.consts.xx, self.consts.yy),
            method="linear",
            fill_value=0)

        act_tuning_value = self.act_tuning_array[self.wt_number - 1]
        masked_surface = dw_mesh * self.consts.mask * act_tuning_value

        return masked_surface


class OapSurface(Surface):
    def __init__(self, constants, radius_of_curvature, off_axis_distance, clocking_angle_rad):
        self.consts = constants
        self.radius_of_curvature = radius_of_curvature
        self.off_axis_distance = off_axis_distance
        self.clocking_angle_rad = clocking_angle_rad

        oap = oap_calculation(
            clocking_angle_rad=self.clocking_angle_rad,
            radius_of_curvature=self.radius_of_curvature,
            off_axis_distance=self.off_axis_distance,
            x_mesh=self.consts.xx,
            y_mesh=self.consts.yy)

        self.surface = self.consts.mask * oap

        self.pv = super()._pv_calculation()
        self.rms = super()._rms_calculation()
        self.volume = super()._volume_calculation()[0]
        self.offset_height_value = super()._volume_calculation()[1]
        self.zernike_value_array = super()._zernike_value_array_calculation(self.surface)

    def h(self):
        mkhelp(self)


class ZernikeToTorque:
    def __init__(
            self,
            constants,
            target_zernike_value_array: ndarray,
            ignore_zernike_number_list: list[int],
            restructed_torque_value: float):

        """ZernikeToTorque
        与えられたzernikeベクトルに制約付き最小二乗fittingするようなtorqueベクトルを計算

        Parameters
        ----------
        constants : [type]
            Constants クラス
        target_zernike_value_array : ndarray[float]
            zernike多項式の値
        ignore_zernike_number_list : list[int]
            WH研磨体積削減で無視するzernike多項式の項番号
            無視しない場合は空のリストを渡す
        restructed_torque_value : float
            トルクの制限値
        """

        self.consts = constants
        self.target_zernike_value_array = target_zernike_value_array
        self.ignore_zernike_number_list = ignore_zernike_number_list
        self.restructed_torque_value = abs(restructed_torque_value)

        if len(self.ignore_zernike_number_list) == 0:
            make_torque_value_array_result = self.__make_torque_value_array(
                operation_matrix=self.consts.operation_matrix,
                zernike_value_array=self.target_zernike_value_array)

            self.torque_value_array = make_torque_value_array_result["torque"]
            self.optimize_result = make_torque_value_array_result["optimize_result"]

        else:
            remaining_operation_matrix = make_remaining_matrix(
                self.consts.operation_matrix,
                self.ignore_zernike_number_list)

            remaining_zernike_value_array = make_remaining_matrix(
                self.target_zernike_value_array,
                self.ignore_zernike_number_list)

            make_torque_value_array_result = self.__make_torque_value_array(
                operation_matrix=remaining_operation_matrix,
                zernike_value_array=remaining_zernike_value_array)

            self.torque_value_array = make_torque_value_array_result["torque"]
            self.optimize_result = make_torque_value_array_result["optimize_result"]

    def h(self):
        mkhelp(self)

    def __make_torque_value_array(
            self,
            operation_matrix: ndarray,
            zernike_value_array: ndarray) -> dict:
        """__make_torque_value_array
        線形モデルb = Axの行列形式に対して制約付き最小二乗フィッティング

        Parameters
        ----------
        operation_matrix : ndarray
            作用行列（線形モデルのAに該当）
        zernike_value_array : ndarray
            zernike係数ベクトル（線形モデルのbに該当）

        Returns
        -------
        dict
            "torque" : fittingに必要なtorque（線形モデルのxに該当）
            "optimize_result" : fittingの詳細結果
        """

        optimize_result = optimize.lsq_linear(
            A=operation_matrix,
            b=zernike_value_array,
            bounds=(-self.restructed_torque_value, self.restructed_torque_value))

        torque_value_array = optimize_result["x"]

        result_dict = {
            "torque": torque_value_array,
            "optimize_result": optimize_result}

        return result_dict


class TorqueToZernike:

    def __init__(
            self,
            constants,
            torque_value_array: ndarray):

        """
        与えられたトルクに作用行列をかけてzernikeに変換
        len(self.zernike_value_array) = Constants.zernike_max_degree

        Parameters
        ----------
        constants : [type]
            Constant クラス
        torque_value_array : ndarray
            トルクベクトル
        """

        self.consts = constants
        self.torque_value_array = torque_value_array

        self.zernike_value_array = np.dot(
            self.consts.operation_matrix,
            self.torque_value_array)

    def h(self):
        mkhelp(self)

    def make_torque_plot(self, figure, position=111):
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

        ax.hlines(0, xmin=1, xmax=12, color="darkgray")

        return ax


class OapConstants:
    def __init__(
            self,
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
        minimize_init_list = [
            self.ideal_radius_of_curvature + self.delta_radius_of_curvature,
            self.ideal_off_axis_distance + self.delta_off_axis_distance,
            self.ideal_clocking_angle_rad + self.delta_clocking_angle_rad]

        return minimize_init_list


class OapMinimize:
    def __init__(self, constants, oap_constants, inputed_surface):
        self.consts = constants
        self.oap_consts = oap_constants
        self.inputed_surface = inputed_surface

        ideal_oap = OapSurface(
            constants=self.consts,
            radius_of_curvature=self.oap_consts.ideal_radius_of_curvature,
            off_axis_distance=self.oap_consts.ideal_off_axis_distance,
            clocking_angle_rad=self.oap_consts.ideal_clocking_angle_rad)

        self.__ideal_oap_surface = ideal_oap.surface

        self.minimize_args_list = self.__make_minimize_args_list()

        print("OapMinimize start")
        start_time = time.time()

        minimize_result = sp.optimize.minimize(
            fun=self.__minimize_input_function,
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
        args_list = [
            self.consts,
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

        test_oap_obj = OapSurface(
            constants=arg_consts,
            radius_of_curvature=test_radius_of_curvature,
            off_axis_distance=test_off_axis_distance,
            clocking_angle_rad=test_clocking_angle_rad)

        test_oap_surface_physical_height = test_oap_obj.surface

        difference_of_surface = inputed_surface_physical_height - test_oap_surface_physical_height

        difference_surface_obj = Surface(
            constants=arg_consts,
            surface=difference_of_surface)

        volume = difference_surface_obj.volume

        del test_oap_obj
        del difference_surface_obj
        return volume


class CirclePathMeasurementCsvReading:
    def __init__(
            self,
            Constants,
            original_csv_fpath: str,
            deformed_csv_fpath: str) -> None:

        self.consts = Constants
        self.df_raw_original = pd.read_csv(original_csv_fpath)

        if deformed_csv_fpath == "":
            self.df_diff = pd.DataFrame({
                "degree": self.df_raw_original["angle"].values,
                "radian": np.deg2rad(self.df_raw_original["angle"].values),
                "height": self.df_raw_original["height"].values})

        else:
            self.df_raw_deformed = pd.read_csv(deformed_csv_fpath)
            height_diff = self.df_raw_deformed["height"].values - self.df_raw_original["height"].values

            self.df_diff = pd.DataFrame({
                "degree": self.df_raw_original["angle"].values,
                "radian": np.deg2rad(self.df_raw_original["angle"].values),
                "height": height_diff})

    def h(self):
        mkhelp(self)


class CirclePathMeasurementTxtReading:

    def __init__(
            self,
            Constants: Constants,
            original_txt_fpath: str,
            deformed_txt_fpath: str) -> None:
        """鍵谷先生の円環パス逐次積分結果を読み込んで、
        CirclePathMeasurementCsvReadingで使っているものと同じ形式のdataframeを作る

        Parameters
        ----------
        Constants : Constants
            自作インスタンス
        original_txt_fpath : str
            変形前の円環パス逐次積分結果のファイルパス（.txt）
        deformed_txt_fpath : str
            変形後の円環パス逐次積分結果のファイルパス（.txt）
        """

        self.consts = Constants

        # データの読み出し
        self.df_raw_original = self.__read_txt(
            txt_filepath=original_txt_fpath)

        self.df_raw_deformed = self.__read_txt(
            txt_filepath=deformed_txt_fpath)

        # 円環パスの半径を計算
        radius_original = np.mean(np.sqrt(
            self.df_raw_original["x"]**2 + self.df_raw_original["y"]**2))

        radius_deformed = np.mean(np.sqrt(
            self.df_raw_deformed["x"]**2 + self.df_raw_deformed["y"]**2))

        self.circle_path_radius = (radius_original + radius_deformed) / 2

        # 高さの差分を出す
        # データ点数が合わない場合、少ない方に合わせてから高さ差分をとりたい
        if len(self.df_raw_original) <= len(self.df_raw_deformed):
            idx_max = len(self.df_raw_original)
        else:
            idx_max = len(self.df_raw_deformed)

        height_diff = self.df_raw_deformed["height"].values[:idx_max] - self.df_raw_original["height"].values[:idx_max]

        self.df_diff = pd.DataFrame({
            "degree": self.df_raw_original["degree"].values[:idx_max],
            "height": height_diff,
        })

    def h(self):
        mkhelp(self)

    def __read_txt(self, txt_filepath: str) -> pd.DataFrame:
        """鍵谷先生の円環パス逐次積分結果を読み込んで、
        CirclePathMeasurementCsvReadingで使っているものと同じ形式のdataframeを作る

        Parameters
        ----------
        txt_filepath : str
            鍵谷先生の作る円環パス逐次積分結果（.txt）

        Returns
        -------
        pd.DataFrame
            angle : [deg] 角度
            height : [m] 高さ
        """

        df_raw = pd.read_csv(
            txt_filepath, skiprows=1, names=["x", "y", "z", "beam"],
            delimiter=" ", index_col=None, skipfooter=2, engine="python")

        x_array = df_raw["x"].values
        y_array = df_raw["y"].values

        angle_raw_array = np.rad2deg(np.arctan2(y_array, x_array))

        # 角度を-180deg to +180degから0deg to +360deg へ変換
        angle_array = np.where(
            angle_raw_array >= 0,
            angle_raw_array.copy(),
            angle_raw_array.copy() + 360
        )

        df = pd.DataFrame({
            "x": df_raw["x"].values * 1e-3,
            "y": df_raw["y"].values * 1e-3,
            "degree": angle_array,
            "height": df_raw["z"].values * 1e-9
        })

        return df


class CirclePathZernikeFitting:
    def __init__(
            self,
            Constants: Constants,
            circle_path_radius: float,
            degree_array: ndarray,
            unprocessed_height_array: ndarray,
            ignore_zernike_number_list: list[int]) -> None:
        """CirclePathZernikeFitting
        zernike多項式を除去

        Parameters
        ----------
        Constants : Constants
            自作インスタンス
        circle_path_radius : float
            [m] 円環パスの半径
        degree_array : ndarray
            [deg] 角度の1次元array
        unprocessed_height_array : ndarray
            [m] WH駆動による高さの変形量の1次元array
        ignore_zernike_number_list : list[int]
            zernikeでfittingした後に除去する項（中身は1~11まで）
        """

        self.consts = Constants
        self.circle_path_radius = circle_path_radius
        self.degree_array = degree_array
        self.unprocessed_height_array = unprocessed_height_array
        self.ignore_zernike_number_list = ignore_zernike_number_list

        self.radian_array = np.deg2rad(self.degree_array)

        zernike_fitting_result = self.__zernike_fitting()
        self.optimize_result = zernike_fitting_result["optimize_result"]
        self.r_const_zernike_polynomial_array = zernike_fitting_result["r_const_zernike_polynomial_array"]
        self.r_const_zernike_number_meaning_list = ["1+4+11", "2+8", "3+7", "5", "6", "9", "10"]

        zernike_removing_result = self.__zernike_removing(
            r_const_zernike_polynomial_array=self.r_const_zernike_polynomial_array)

        self.removing_zernike_height_array = zernike_removing_result["removing_zernike_height"]
        self.zernike_removed_height_array = zernike_removing_result["residual"]

    def h(self):
        mkhelp(self)

    def __r_const_zernike_polynomial_calculation(
            self,
            r_const_coef: list[float],
            pupil_radius: float,
            radius: float,
            theta: ndarray) -> ndarray:
        """
        radiusが固定値（＝円環パス）の場合のzernike 多項式の計算
        そのままのzernikeでフィッティングすると、
        thetaによる計算結果の部分が全く同じ値になる項（1, 4, 11など）が出て各係数が独立ではなくなってしまう

        Parameters
        ----------
        r_const_coef : list[float]
            radius固定値の時のzernike係数ベクトル（z=11まででlen=7）
        pupil_radius : float
            [m] 瞳の半径
        radius : float
            [m] 円環パスの半径（固定値なのでarrayではない）
        theta : ndarray
            [rad] 角度

        Returns
        -------
        ndarray
            zernike 多項式の計算結果
        """

        zernike1_ = zernike_term_calculation(1, pupil_radius, radius, theta)
        zernike2_ = zernike_term_calculation(2, pupil_radius, radius, theta)
        zernike3_ = zernike_term_calculation(3, pupil_radius, radius, theta)
        zernike4_ = zernike_term_calculation(4, pupil_radius, radius, theta)
        zernike5_ = zernike_term_calculation(5, pupil_radius, radius, theta)
        zernike6_ = zernike_term_calculation(6, pupil_radius, radius, theta)
        zernike7_ = zernike_term_calculation(7, pupil_radius, radius, theta)
        zernike8_ = zernike_term_calculation(8, pupil_radius, radius, theta)
        zernike9_ = zernike_term_calculation(9, pupil_radius, radius, theta)
        zernike10_ = zernike_term_calculation(10, pupil_radius, radius, theta)
        zernike11_ = zernike_term_calculation(11, pupil_radius, radius, theta)

        zernike1_4_11 = r_const_coef[0] * (zernike1_ + zernike4_ + zernike11_)
        zernike2_8 = r_const_coef[1] * (zernike2_ + zernike8_)
        zernike3_7 = r_const_coef[2] * (zernike3_ + zernike7_)
        zernike5 = r_const_coef[3] * zernike5_
        zernike6 = r_const_coef[4] * zernike6_
        zernike9 = r_const_coef[5] * zernike9_
        zernike10 = r_const_coef[6] * zernike10_

        r_const_zernike_polynomial = zernike1_4_11 + zernike2_8 + zernike3_7 + zernike5 + zernike6 + zernike9 + zernike10

        return r_const_zernike_polynomial

    def __zernike_fitting(self):

        def minimize_funciton_sq(x, params_):
            pupil_radius_, radius_, theta_array_, height_array_ = params_
            zernike_array_ = self.__r_const_zernike_polynomial_calculation(
                r_const_coef=x,
                pupil_radius=pupil_radius_,
                radius=radius_,
                theta=theta_array_)

            residual = height_array_ - zernike_array_
            return residual

        pupil_radius = self.consts.varid_radius
        radius = self.circle_path_radius
        theta_array = self.radian_array
        height_array = self.unprocessed_height_array

        params = [pupil_radius, radius, theta_array, height_array]

        optimize_result_sq = optimize.least_squares(fun=minimize_funciton_sq,
                                                    x0=(np.ones(7) * 1e-8),
                                                    args=(params, ))

        result_dict = {
            "optimize_result": optimize_result_sq,
            "r_const_zernike_polynomial_array": optimize_result_sq["x"]}

        return result_dict

    def __zernike_removing(
            self,
            r_const_zernike_polynomial_array: ndarray) -> dict:
        """
        zernike係数ベクトルを用いてzernike成分を除去

        Parameters
        ----------
        r_const_zernike_polynomial_array : ndarray
            半径一定の場合のzernike係数ベクトル

        Returns
        -------
        dict
            removing_zernike_height: 除去するzernike成分の計算結果のarray
            residual: height - removing_zernike
        """

        def make_removing_r_const_zernike_polynomial_array(
                r_const_zernike_polynomial_array_: ndarray,
                ignore_zernike_number_list: list[int]) -> list[int]:

            removing_r_const_zernike_polynomial_tf = np.zeros(7)

            # 半径一定だとzernike 1, 4, 11項を区別できない
            if 1 in ignore_zernike_number_list or 4 in ignore_zernike_number_list or 11 in ignore_zernike_number_list:
                # もし ignore_zernike_number_list に 1, 4, 11 のどれかが含まれる場合
                # 1, 4, 11 は区別できないので、除くのであれば 1, 4, 11 全てをまとめて除く必要がある

                if 1 in ignore_zernike_number_list and 4 in ignore_zernike_number_list and 11 in ignore_zernike_number_list:
                    removing_r_const_zernike_polynomial_tf[0] = 1

                else:
                    print("Error!")
                    print("When radius is const value,")
                    print("if you input 1 or 4 or 11 for ignore_zernike_number_list,")
                    print("you must input all of [1, 4, 11]")

                    # Errorに気づくようにnp.nanを入れる
                    removing_r_const_zernike_polynomial_tf[0] = np.nan

            # 半径一定だとzernike 2, 8 項を区別できない
            if 2 in ignore_zernike_number_list or 8 in ignore_zernike_number_list:
                # もし ignore_zernike_number_list に 2, 8 のどれかが含まれる場合
                # 2, 8 は区別できないので、除くのであれば 2, 8 両方をまとめて除く必要がある

                if 2 in ignore_zernike_number_list and 8 in ignore_zernike_number_list:
                    removing_r_const_zernike_polynomial_tf[1] = 1

                else:
                    print("Error!")
                    print("When radius is const value,")
                    print("if you input 2 or 8 for ignore_zernike_number_list,")
                    print("you must input both of [2, 8]")

                    # Errorに気づくようにnp.nanを入れる
                    removing_r_const_zernike_polynomial_tf[1] = np.nan

            # 半径一定だとzernike 2, 8 項を区別できない
            if 3 in ignore_zernike_number_list or 7 in ignore_zernike_number_list:
                # もし ignore_zernike_number_list に 3, 7 のどれかが含まれる場合
                # 3, 7 は区別できないので、除くのであれば 3, 7 両方をまとめて除く必要がある

                if 3 in ignore_zernike_number_list and 7 in ignore_zernike_number_list:
                    removing_r_const_zernike_polynomial_tf[2] = 1

                else:
                    print("Error!")
                    print("When radius is const value,")
                    print("if you input 3 or 7 for ignore_zernike_number_list,")
                    print("you must input both of [3, 7]")

                    # Errorに気づくようにnp.nanを入れる
                    removing_r_const_zernike_polynomial_tf[2] = np.nan

            if 5 in ignore_zernike_number_list:
                removing_r_const_zernike_polynomial_tf[3] = 1

            if 6 in ignore_zernike_number_list:
                removing_r_const_zernike_polynomial_tf[4] = 1

            if 9 in ignore_zernike_number_list:
                removing_r_const_zernike_polynomial_tf[5] = 1

            if 10 in ignore_zernike_number_list:
                removing_r_const_zernike_polynomial_tf[6] = 1

            removing_r_const_zernike_polynomial_array_ = r_const_zernike_polynomial_array_ * removing_r_const_zernike_polynomial_tf

            return removing_r_const_zernike_polynomial_array_

        radian_array = self.radian_array
        height_array = self.unprocessed_height_array

        removing_r_const_zernike_polynomial_array = make_removing_r_const_zernike_polynomial_array(
            r_const_zernike_polynomial_array_=r_const_zernike_polynomial_array,
            ignore_zernike_number_list=self.ignore_zernike_number_list
        )

        removing_zernike_height_array = self.__r_const_zernike_polynomial_calculation(
            r_const_coef=removing_r_const_zernike_polynomial_array,
            pupil_radius=self.consts.varid_radius,
            radius=self.circle_path_radius,
            theta=radian_array)

        residual_height_array = height_array - removing_zernike_height_array

        result_dict = {
            "removing_zernike_height": removing_zernike_height_array,
            "residual": residual_height_array}

        return result_dict
