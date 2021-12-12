# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 14:59:35 2021

@author: swimc
"""

import numpy as np
import proper as pr
import matplotlib.pyplot as plt
import scipy as sp

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
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

class Constants:
    def __init__(self, physical_radius, ignore_radius, pixel_number, zernike_max_degree):
        self.physical_radius = physical_radius
        self.ignore_radius = ignore_radius
        self.pixel_number = pixel_number

        self.varid_radius = physical_radius - ignore_radius
        self.xx, self.yy = make_meshgrid(-physical_radius, physical_radius, -physical_radius, physical_radius, pixel_number)
        self.tf = np.where(self.xx**2+self.yy**2<=self.varid_radius**2, True, False)
        self.mask = np.where(self.tf==True, 1, np.nan)
        self.zernike_max_degree = zernike_max_degree

class ZernikeToSurface:
    def __init__(self, constants, zernike_number_list, zernike_value_array):
        self.__constants = constants
        self.zernike_number_list = zernike_number_list
        self.zernike_value_array = zernike_value_array

        self.surface = self.__make_masked_zernike_surface()
        self.pv = pv_calculation(self.surface)
        self.rms = rms_calculation(self.surface)
    
    def __make_masked_zernike_surface(self):
        optical_wavelength = 500e-9
        
        wavestruct = pr.prop_begin(beam_diameter = 2*self.__constants.varid_radius,
                                   lamda = optical_wavelength,
                                   grid_n = self.__constants.pixel_number,
                                   beam_diam_fraction=1)
        
        wfe = pr.prop_zernikes(wavestruct, 
                               self.zernike_number_list, 
                               self.zernike_value_array)
        
        masked_wfe = self.__constants.mask * wfe
        return masked_wfe
    
    def make_image_plot(self, fig=False, position=False, color_scale=False, cbar_min_percent=0, cbar_max_percent=100, pv_digits=2, rms_digits=2):
        figure = plt.figure()
        position = 111
        
        cmap = cm.jet
        fontsize = 15
        title = "pv = " + str(round(self.pv*1e6, 2)) + " [um]" + "\n" \
              + "RMS = " + str(round(self.rms*1e6, 2)) + " [um]"
        
        cbar_min = np.nanmin(self.surface) + self.pv * cbar_min_percent/100
        cbar_max = np.nanmin(self.surface) + self.pv * cbar_max_percent/100
        
        extent = [-self.__constants.physical_radius, self.__constants.physical_radius,
                  -self.__constants.physical_radius, self.__constants.physical_radius]
        
        ax = figure.add_subplot(position)
        ax.imshow(self.surface, interpolation="nearest", cmap=cmap, vmin=cbar_min, vmax=cbar_max, origin="lower", extent=extent)
        ax.set_title(title, fontsize=fontsize)
        
        divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
        cax = divider.append_axes("right", "5%", pad="3%")
        
        norm = Normalize(vmin=cbar_min, vmax=cbar_max)
        cbar_title = r"$[\mu$m]"
        mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
        
        cbar = figure.colorbar(mappable, ax=ax, cax=cax)
        cbar.set_label(cbar_title, fontsize=fontsize)
        return ax
        

class WhReproduction:
    def __init__(self, constants, target_zernike_number_list, target_zernike_value_array, restructed_torque_value, ignore_zernike_number_list):
        self.__constants = constants
        self.target_zernike_number_list = target_zernike_number_list
        self.target_zernike_value_array = target_zernike_value_array
        self.restructed_torque_value = abs(restructed_torque_value)
        self.ignore_zernike_number_list = ignore_zernike_number_list
     
        self.operation_matrix = np.genfromtxt("raw_data/WT06_zer10_operation_matrix[m].csv", delimiter=",").T
        self.full_zernike_value_array = self.__make_full_zernike_value_array()
        self.remaining_operation_matrix = self.__make_remaining_matrix(self.operation_matrix)
        self.remaining_zernike_value_array = self.__make_remaining_matrix(self.full_zernike_value_array)
        self.remaining_zernike_number_list = self.__make_remaining_matrix(1+np.arange(self.__constants.zernike_max_degree))
        
        self.torque_value_array, self.restructed_torque_value_array = self.__make_torque_value_array()
        self.remaining_reproducted_zernike_value_array = self.__make_reproducted_zernike_value_array()
        

    def __make_full_zernike_value_array(self):
        target_zernike_number_idx_array = np.array(self.target_zernike_number_list) - 1
        
        full_zernike_value_array = np.zeros(self.__constants.zernike_max_degree)
        for i in range(len(self.target_zernike_value_array)):
            target_zernike_number_idx = target_zernike_number_idx_array[i]
            full_zernike_value_array[target_zernike_number_idx] = self.target_zernike_value_array[i]
        return full_zernike_value_array
    
    def __make_remaining_matrix(self, matrix):
        idx_array = np.array(self.ignore_zernike_number_list) - 1
        remaining_matrix = np.delete(arr=matrix, obj=idx_array, axis=0)
        return remaining_matrix
    
    def __make_torque_value_array(self):
        inverse_operation_matrix = sp.linalg.pinv2(1e9*self.remaining_operation_matrix)*1e9
        
        torque_value_array = np.dot(inverse_operation_matrix, self.remaining_zernike_value_array)
        
        
        only_max_restructed_torque_value_array = np.where(torque_value_array<self.restructed_torque_value,
                                                         torque_value_array,
                                                         self.restructed_torque_value)
            
        restructed_torque_value_array = np.where(only_max_restructed_torque_value_array>-self.restructed_torque_value,
                                                 only_max_restructed_torque_value_array,
                                                 -self.restructed_torque_value)
        
        return torque_value_array, restructed_torque_value_array
    
    def __make_reproducted_zernike_value_array(self):
        remaining_reproducted_zernike_value_array = np.dot(self.remaining_operation_matrix, self.restructed_torque_value_array)
        return remaining_reproducted_zernike_value_array
    
    def make_torque_plot(self,fig=False, position=False, title=False):
        figure = plt.figure()
        position = 111
        
        fontsize = 15
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
        
        ax.set_title(title, fontsize=fontsize)
        ax.grid()
        ax.set_xlabel("Motor Number", fontsize=fontsize)
        ax.set_xticks(x)
        ax.set_xticklabels(x_str)
        ax.set_ylabel("Motor drive amount [mm]", fontsize=fontsize)
        
        ax.hlines([self.restructed_torque_value, -self.restructed_torque_value],
                  xmin=1, xmax=12, color ="red", linestyle="dashed")
        ax.hlines(0, xmin=1, xmax=12, color ="darkgray")
        
        return ax
        

if __name__ == "__main__":
    
    consts = Constants(physical_radius=925e-3, 
                       ignore_radius=25e-3,
                       pixel_number=256,
                       zernike_max_degree = 10)
    
    zernike = ZernikeToSurface(constants = consts, 
                               zernike_number_list = [2],
                               zernike_value_array = np.array([2e-6])) 
    
    idealized_zernike = ZernikeToSurface(constants = consts, 
                                         zernike_number_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                                         zernike_value_array = np.array([ 8.11304455e-08, -1.01586202e-07, -3.38411547e-07,  3.02566783e-07, 2.10233957e-07, -2.01693302e-07, -6.40135092e-08,  1.15529214e-08,3.01199936e-07, -1.78044987e-08]))
                
    
    reproduction_same_for_cb = WhReproduction(constants = consts,
                                              target_zernike_number_list = idealized_zernike.zernike_number_list,
                                              target_zernike_value_array = idealized_zernike.zernike_value_array,
                                              restructed_torque_value=5,
                                              ignore_zernike_number_list=[1,2,3,4])
    
    reproducted_zernike = ZernikeToSurface(constants = consts,
                                           zernike_number_list = reproduction_same_for_cb.remaining_zernike_number_list,
                                           zernike_value_array = reproduction_same_for_cb.remaining_reproducted_zernike_value_array)