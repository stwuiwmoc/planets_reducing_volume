# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 14:59:35 2021

@author: swimc
"""

import numpy as np
import proper as pr
import matplotlib.pyplot as plt
import scipy as sp

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

class ZernikeSurface:
    def __init__(self, constants, zernike_number_list, zernike_value_array):
        self.__constants = constants
        self.zernike_number_list = zernike_number_list
        self.zernike_value_array = zernike_value_array

        self.surface = self.make_masked_zernike_surface()
        self.pv = pv_calculation(self.surface)
        self.rms = rms_calculation(self.surface)
    
    def make_masked_zernike_surface(self):
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

class WhReproductedSurface:
    def __init__(self, constants, target_surface, torque_max_value, ignore_zernike_number_list):
        self.__constants = constants
        self.target_surface = target_surface
        self.torque_max_value = torque_max_value
        self.ignore_zernike_number_list = ignore_zernike_number_list
     
        self.operation_matrix = np.genfromtxt("raw_data/WT06_zer10_operation_matrix[m].csv", delimiter=",").T
        self.fitted_zernike_value_array = self.make_fitted_zernike_value_array()
        self.__remaining_operation_matrix = self.make_remaining_matrix(self.operation_matrix)
        self.__remaining_fitted_zernike_value_array = self.make_remaining_matrix(self.fitted_zernike_value_array)
        self.__remaining_zernike_number_list = self.make_remaining_matrix(1+np.arange(self.__constants.zernike_max_degree))
        self.torque_value_array, self.restructed_torque_value_array = self.make_torque_value_array()
        
    
    def make_fitted_zernike_value_array(self):
        fitted_zernike_value_array = pr.prop_fit_zernikes(wavefront0=self.target_surface,
                                                          pupil0=self.__constants.tf,
                                                          pupilradius0=self.__constants.pixel_number/2,
                                                          nzer=self.__constants.zernike_max_degree,
                                                          xc=self.__constants.pixel_number/2,
                                                          yc=self.__constants.pixel_number/2)
        return fitted_zernike_value_array
    
    def make_remaining_matrix(self, matrix):
        idx_array = np.array(self.ignore_zernike_number_list) - 1
        remaining_matrix = np.delete(arr=matrix, obj=idx_array, axis=0)
        return remaining_matrix
    
    def make_torque_value_array(self):
        inverse_operation_matrix = sp.linalg.pinv2(1e9*self.__remaining_operation_matrix)*1e9
        print(inverse_operation_matrix)
        
        #torque_value_array = np.dot(inverse_operation_matrix, self.__remaining_fitted_zernike_value_array)
        torque_value_array = np.dot(inverse_operation_matrix, np.array([-1.01586202e-07, -3.38411547e-07,  3.02566783e-07,  2.10233957e-07,-2.01693302e-07, -6.40135092e-08,  1.15529214e-08,  3.01199936e-07,-1.78044987e-08]))
        
        print(torque_value_array)
        restructed_torque_value_array = np.where(torque_value_array<self.torque_max_value,
                                                 torque_value_array,
                                                 self.torque_max_value)
        return torque_value_array, restructed_torque_value_array
    
    def make_reproducted_surface(self):
        remaining_reproducted_zernike_value_array = np.dot(self.__remaining_operation_matrix, self.torque_value_array)
        

if __name__ == "__main__":
    
    consts = Constants(physical_radius=925e-3, 
                       ignore_radius=25e-3,
                       pixel_number=256,
                       zernike_max_degree = 10)
    
    zernike = ZernikeSurface(constants = consts, 
                             zernike_number_list = [2],
                             zernike_value_array = np.array([2e-6])) 
    
    reprod = WhReproductedSurface(constants = consts,
                                  target_surface=zernike.surface,
                                  torque_max_value=1e3,
                                  ignore_zernike_number_list=[1])
   