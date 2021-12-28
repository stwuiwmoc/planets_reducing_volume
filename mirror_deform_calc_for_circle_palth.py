# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 14:59:35 2021

@author: swimc
"""

import numpy as np
import proper as pr
import matplotlib.pyplot as plt
import scipy as sp
import cv2
import PIL

from matplotlib import cm
from matplotlib.colors import Normalize
import mpl_toolkits.axes_grid1


def mkfolder(suffix = ""):
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
        if attr.startswith("_"): continue
        print(attr)
    for method in inspect.getmembers(instance, inspect.ismethod):
        if method[0].startswith("_"): continue
        print(method[0]+"()")


def make_meshgrid(x_min, x_max, y_min, y_max, pixel_number):
    x_array = np.linspace(x_min, x_max, pixel_number)
    y_array = np.linspace(y_min, y_max, pixel_number)
    return np.meshgrid(x_array, y_array)

def pv_calculation(array2d):
    peak = np.nanmax(array2d)
    valley = np.nanmin(array2d)
    pv = peak - valley
    return pv

def rms_calculation(array2d):
    sigma = np.nansum(array2d**2)
    data_count = np.sum(~np.isnan(array2d))
    rms = np.sqrt( sigma / data_count )
    return rms

def make_remaining_matrix(matrix, ignore_zernike_number_list):
    idx_array = np.array(ignore_zernike_number_list) - 1
    remaining_matrix = np.delete(arr=matrix, obj=idx_array, axis=0)
    return remaining_matrix

def make_full_torque_value_array(torque_number_list, torque_value_aray):
    full_value_array = np.zeros(36)
    idx_array = np.array(torque_number_list) -1
    for i in range(len(idx_array)):
        idx = idx_array[i]
        full_value_array[idx] = torque_value_aray[i]
    return full_value_array


class Constants:
    
    
    def __init__(self, physical_radius, ignore_radius, 
                 pixel_number, zernike_max_degree):
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

        Returns
        -------
        None.

        """
        self.physical_radius = physical_radius
        self.ignore_radius = ignore_radius
        self.pixel_number = pixel_number

        self.varid_radius = physical_radius - ignore_radius
        self.xx, self.yy = make_meshgrid(-physical_radius, physical_radius, 
                                         -physical_radius, physical_radius, 
                                         pixel_number)
        self.tf = np.where(self.xx**2+self.yy**2<=self.varid_radius**2, True, False)
        self.mask = np.where(self.tf==True, 1, np.nan)
        self.zernike_max_degree = zernike_max_degree
        self.operation_matrix = np.genfromtxt("raw_data/WT06_zer10_operation_matrix[m].csv", delimiter=",").T
        
    def h(self):
        mkhelp(self)

class ZernikeToSurface:
    
    
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

        Returns
        -------
        None.

        """
        
        self.consts = constants
        self.zernike_number_list = zernike_number_list
        self.zernike_value_array = zernike_value_array

        self.surface = self.__make_masked_zernike_surface()
        self.pv=pv_calculation(self.surface)
        self.rms = rms_calculation(self.surface)
        
    def h(self):
        mkhelp(self)
    
    def __make_masked_zernike_surface(self):
        optical_wavelength = 500e-9
        
        wavestruct = pr.prop_begin(beam_diameter = 2*self.consts.varid_radius,
                                   lamda = optical_wavelength,
                                   grid_n = self.consts.pixel_number,
                                   beam_diam_fraction=1)
        
        wfe = pr.prop_zernikes(wavestruct, 
                               self.zernike_number_list, 
                               self.zernike_value_array)
        
        masked_wfe = self.consts.mask * wfe
        return masked_wfe
    
    def make_image_plot(self, figure=plt.figure(), position=111, 
                        color_scale=False, cbar_min_percent=0, cbar_max_percent=100, 
                        pv_digits=2, rms_digits=2):
        
        cmap = cm.jet
        fontsize = 15
        title = "pv = " + str(round(self.pv*1e6, pv_digits)) + " [um]" + "\n" \
              + "RMS = " + str(round(self.rms*1e6, rms_digits)) + " [um]"
        
        cbar_min = np.nanmin(self.surface) + self.pv * cbar_min_percent/100
        cbar_max = np.nanmin(self.surface) + self.pv * cbar_max_percent/100
        
        extent = [-self.consts.physical_radius, self.consts.physical_radius,
                  -self.consts.physical_radius, self.consts.physical_radius]
        
        ax = figure.add_subplot(position)
        ax.imshow(self.surface, interpolation="nearest", cmap=cmap, 
                  vmin=cbar_min, vmax=cbar_max, origin="lower", extent=extent)
        ax.set_title(title, fontsize=fontsize)
        
        divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
        cax = divider.append_axes("right", "5%", pad="3%")
        
        norm = Normalize(vmin=cbar_min, vmax=cbar_max)
        cbar_title = "[mm]"
        mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
        
        cbar = figure.colorbar(mappable, ax=ax, cax=cax)
        cbar.set_label(cbar_title, fontsize=fontsize)
        return ax
        
    def make_circle_path_plot(self, figure=plt.figure(), position=111, 
                              radius=0.870, height_magn=1e9, height_unit_str="[nm]"):
        fontsize = 15
        varid_radius_pixel_number = int(self.consts.varid_radius/self.consts.physical_radius*self.consts.pixel_number/2)
        measurement_radius_idx = int(radius*1e3)
        
        image = np.where(self.consts.tf==True, self.surface, 0)
        flags = cv2.INTER_CUBIC + cv2.WARP_FILL_OUTLIERS + cv2.WARP_POLAR_LINEAR
        
        linear_polar_image = cv2.warpPolar(src=image,
                                           dsize=(int(self.consts.varid_radius*1e3),360), 
                                           center=(self.consts.pixel_number/2, self.consts.pixel_number/2),
                                           maxRadius=varid_radius_pixel_number, 
                                           flags=flags)
        
        circle_path_line = height_magn*linear_polar_image[:, measurement_radius_idx]
        height_pv = np.nanmax(circle_path_line) - np.nanmin(circle_path_line)
        height_pv_str = str(round(height_pv, 2))
        
        self.make_image_plot()
        
        ax = figure.add_subplot(position)
        ax.plot(circle_path_line)
        ax.grid()
        ax.set_title("circle_path ( pv = " + height_pv_str + " " + height_unit_str + " )", fontsize=fontsize)
        ax.set_xlabel("degree", fontsize=fontsize)
        ax.set_ylabel("heignt " + height_unit_str + "\nat R=" + str(measurement_radius_idx), fontsize=fontsize)
        return ax
    
class StitchedCsvToSurface(ZernikeToSurface):
    def __init__(self, constants, original_stitched_csv_fpath, deformed_stitched_csv_fpath):
        """
        class : 2d-surface from measured (stitched) 2d-csv data

        Parameters
        ----------
        constants : TYPE
            instance by class : Constants
        original_stitched_csv_fpath : str
            filepath of measured (stitched) and not deformed 2d-csv data 
        deformed_stitched_csv_fpath : TYPE
            filepath of measured (stitched) and deformed 2d-csv data

        Returns
        -------
        None.

        """
        self.consts = constants
    
        self.original_masked_surface = self.__read_csv_to_masked_surface(original_stitched_csv_fpath)
        self.deformed_masked_surface = self.__read_csv_to_masked_surface(deformed_stitched_csv_fpath)
        self.surface = self.deformed_masked_surface - self.original_masked_surface
        
        self.pv=pv_calculation(self.surface)
        self.rms=rms_calculation(self.surface)
    
    def __read_csv_to_masked_surface(self,filepath):
        raw = np.loadtxt(filepath)
        raw_zero_fill = np.where(np.isnan(raw), 0, raw)
        image = PIL.Image.fromarray(raw_zero_fill)
        img_resize = image.resize(size=(self.consts.pixel_number, self.consts.pixel_number))
        masked_surface = self.consts.mask * np.array(img_resize)
        return masked_surface
        
class ZernikeToTorque:
    def __init__(self, constants, target_zernike_number_list, target_zernike_value_array, ignore_zernike_number_list):
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
        self.remaining_zernike_number_list = make_remaining_matrix(1+np.arange(self.consts.zernike_max_degree), self.ignore_zernike_number_list)
        
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
        inverse_operation_matrix = sp.linalg.pinv(1e9*self.remaining_operation_matrix)*1e9
        
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
        self.remaining_zernike_number_list = make_remaining_matrix(1+np.arange(self.consts.zernike_max_degree), 
                                                                   self.ignore_zernike_number_list)

    def h(self):
        mkhelp(self)

    def __make_restructed_torque_value_array(self):
        only_max_restructed_torque_value_array = np.where(self.torque_value_array<self.restructed_torque_value,
                                                         self.torque_value_array,
                                                         self.restructed_torque_value)
        
        restructed_torque_value_array = np.where(only_max_restructed_torque_value_array>-self.restructed_torque_value,
                                                 only_max_restructed_torque_value_array,
                                                 -self.restructed_torque_value)
        
        return restructed_torque_value_array
    
    def __make_reproducted_zernike_value_array(self, using_torque_value_array):
        remaining_reproducted_zernike_value_array = np.dot(self.remaining_operation_matrix, using_torque_value_array)
        return remaining_reproducted_zernike_value_array
    
    def make_torque_plot(self, figure=plt.figure(), position=111):
        fontsize = 15
        ax_title = "WH (black : 1-12, green : 13-24, violet : 25-36)"
        
        x = np.arange(1, 13)
        x_str = []
        
        for i in range(12):
            text = ""
            for j in range(3):
                num = str(12*j + i + 1).zfill(2)
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
                  xmin=1, xmax=12, color ="red", linestyle="dashed")
        ax.hlines(0, xmin=1, xmax=12, color ="darkgray")
        
        return ax

    
if __name__ == "__main__":
    
    CONSTS = Constants(physical_radius=925e-3, 
                       ignore_radius=25e-3,
                       pixel_number=256,
                       zernike_max_degree = 10)
    
    """
    # for parameter study
    original_torque_value_array = make_full_torque_value_array([1,7,13,19,25,31],
                                                               [5,-5,5,-5,5,-5])
    original_torque_value_array = make_full_torque_value_array([2,8,14,20,26,32],
                                                               [5,-5,5,-5,5,-5])
    original_torque_value_array = make_full_torque_value_array([4,10,16,22,28,34],
                                                               [5,-5,5,-5,5,-5])
    original_torque_value_array = make_full_torque_value_array([5,11,17,23,29,35],
                                                               [5,-5,5,-5,5,-5])    
    original_torque_value_array = make_full_torque_value_array([6,12,18,24,30,36],
                                                               [5,-5,5,-5,5,-5])
    
    original_torque_value_array = make_full_torque_value_array([1,7,13,19,25,31, 3,9,15,21,27,33, 4,10,16,22,28,34, 5,11,17,23,29,35, 6,12,18,24,30,36],
                                                               [-5,5,-5,5,-5,5, -5,5,-5,5,-5,5, -5,5,-5,5,-5,5, -5,5,-5,5,-5,5, 5,-5,5,-5,5,-5])
    """


    original_torque_value_array = make_full_torque_value_array([3,9,15,21,27,33],
                                                               [5,-5,5,-5,5,-5])
    wh_deformed_zernike = TorqueToZernike(constants=CONSTS,
                                          torque_value_array=original_torque_value_array,
                                          restructed_torque_value=5,
                                          ignore_zernike_number_list=[1,2,3,4,5,6])
    
    wh_deformed_surface = ZernikeToSurface(constants=CONSTS,
                                           zernike_number_list=wh_deformed_zernike.remaining_zernike_number_list,
                                           zernike_value_array=wh_deformed_zernike.remaining_reproducted_zernike_value_array)
    
    fig = plt.figure(figsize=(15,7))
    gs = fig.add_gridspec(2,3)
    wh_deformed_zernike.make_torque_plot(figure=fig, position=gs[0,2])
    wh_deformed_surface.make_image_plot(figure=fig, position=gs[0:2,0:2])
    wh_deformed_surface.make_circle_path_plot(figure=fig, position=gs[1,2])
    fig.tight_layout()
    
    
    diff = StitchedCsvToSurface(constants=CONSTS,
                                original_stitched_csv_fpath="mkfolder/stitch2mesh/zer10_1215xm1301214ym870-510cir.v4.22.hei_dense.csv", 
                                deformed_stitched_csv_fpath="mkfolder/stitch2mesh/zer10_1215xm1301216ym870-510cir.v4.21.hei_dense.csv")
    
    diff.make_image_plot()
    diff.make_circle_path_plot()