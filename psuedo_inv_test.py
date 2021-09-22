# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 19:00:00 2021

@author: swimc
"""


import numpy as np

def psuedo(array):
    a = array
    a_t = a.T
    a_dot = np.dot(a,a_t)
    a_inv = np.linalg.inv(a_dot)
    
    a_pinv = np.dot(a_t, a_inv)
    return a_pinv

if __name__ == '__main__':
    a = 10**3 * np.genfromtxt("reduce_polishing_code_pparc\WT03_zer10_opration_matrix[m].csv", delimiter=",")
    b = a * 10**6
    
    #a_pinv = psuedo(a)
    a_pinv = np.linalg.pinv(a)
    b_pinv = np.linalg.pinv(b) * 10**6
    
    test1 = np.mean( a - np.dot(a, np.dot(a_pinv,a)) )
    test2 = np.mean( a_pinv - np.dot(a_pinv, np.dot(a, a_pinv)) )
    test3 = np.mean( np.dot(a, a_pinv).conj().T - np.dot(a,a_pinv) )
    test4 = np.mean( np.dot(a_pinv, a).conj().T - np.dot(a_pinv, a) )