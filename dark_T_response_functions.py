# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 12:18:25 2015

@author: imchugh
"""

from scipy.optimize import curve_fit
import numpy as np

# No fixed parameters    
def TRF(data_d,Eo,rb):
    return rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))

# Eo fixed
def make_TRF(Eo):
    def TRF(data_d,rb):
        return rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))
    return TRF

# rb and Eo
def optimise_all(data_dict, params_dict):

    # Initialise error state variable
    error_state = 0              

    drivers_dict = {driver: data_dict[driver] for driver in ['TempC']}
    response_array = data_dict['NEE']

    try:
        params = curve_fit(TRF, drivers_dict, response_array, 
                           p0 = [params_dict['Eo_prior'], 
                                 params_dict['rb_prior']])[0]
    except RuntimeError:
        params = [np.nan, np.nan]
        error_state = 3

    return params, error_state        

# rb   
def optimise_rb(data_dict, params_dict):

    # Initialise error state variable
    error_state = 0              
    
    # Get drivers and response
    drivers_d = {driver: data_dict[driver] for driver in ['TempC']}
    response_array = data_dict['NEE']        
    
    try:
        params = curve_fit(make_TRF(params_dict['Eo_default']), 
                           drivers_d, response_array, 
                           p0 = [params_dict['rb_prior']])[0]
    except RuntimeError:
        params = [np.nan]

    # If negative rb returned, set to nan
    if params[0] < 0:
        error_state = 1
        params = [np.nan]
       
    return params, error_state
