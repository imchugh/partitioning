# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 10:54:05 2015

@author: imchugh
"""
from scipy.optimize import curve_fit
import numpy as np
import pdb

#------------------------------------------------------------------------------
# Data optimisation algorithms

# No fixed parameters
def LRF(data_d,Eo,rb,alpha,Aopt,k):
    
    Aopt_VPD = Aopt * np.exp(-k * (data_d['VPD'] - 1))
    index = np.where(data_d['VPD'] <= 1)[0]
    Aopt_VPD[index] = Aopt
    GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + 
           (alpha * data_d['PAR'] / Aopt_VPD))
    index = np.where(data_d['PAR'] < 10)[0]
    GPP[index] = 0
    Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))
    return GPP, Reco

# Eo fixed
def make_LRF_1(Eo):
    def LRF(data_d,rb,alpha,Aopt,k):
        Aopt_VPD = Aopt * np.exp(-k * (data_d['VPD'] - 1))
        index = np.where(data_d['VPD'] <= 1)[0]
        Aopt_VPD[index] = Aopt
        GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + 
               (alpha * data_d['PAR'] / Aopt_VPD))
        index = np.where(data_d['PAR'] < 10)[0]
        GPP[index] = 0
        Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))
        return GPP + Reco
    return LRF

# Eo and k fixed
def make_LRF_2(Eo, k):
    def LRF(data_d,rb,alpha,Aopt):
        Aopt_VPD = Aopt * np.exp(-k * (data_d['VPD'] - 1))
        index = np.where(data_d['VPD'] <= 1)[0]
        Aopt_VPD[index] = Aopt
        GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + 
               (alpha * data_d['PAR'] / Aopt_VPD))
        index = np.where(data_d['PAR'] < 10)[0]
        GPP[index] = 0
        Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))
        return GPP + Reco
    return LRF

# Eo, k and alpha fixed
def make_LRF_3(Eo, k, alpha):
    def LRF(data_d,rb,Aopt):
        Aopt_VPD = Aopt * np.exp(-k * (data_d['VPD'] - 1))
        index = np.where(data_d['VPD'] <= 1)[0]
        Aopt_VPD[index] = Aopt
        GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + 
               (alpha * data_d['PAR'] / Aopt_VPD))
        index = np.where(data_d['PAR'] < 10)[0]
        GPP[index] = 0
        Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))
        return GPP + Reco
    return LRF

# Eo and rb fixed
def make_LRF_4(Eo, rb):
    def LRF(data_d,alpha,Aopt,k):
#        return uni_LRF(data_d, Eo, rb, alpha, Aopt, k)
        Aopt_VPD = Aopt * np.exp(-k * (data_d['VPD'] - 1))
        index = np.where(data_d['VPD'] <= 1)[0]
        Aopt_VPD[index] = Aopt
        GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + 
               (alpha * data_d['PAR'] / Aopt_VPD))
        index = np.where(data_d['PAR'] < 10)[0]
        GPP[index] = 0
        Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))
        return GPP + Reco
    return LRF

# Eo, rb and k fixed
def make_LRF_5(Eo, rb, k):
    def LRF(data_d,alpha,Aopt):
        Aopt_VPD = Aopt * np.exp(-k * (data_d['VPD'] - 1))
        index = np.where(data_d['VPD'] <= 1)[0]
        Aopt_VPD[index] = Aopt
        GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + 
               (alpha * data_d['PAR'] / Aopt_VPD))
        index = np.where(data_d['PAR'] < 10)[0]
        GPP[index] = 0
        Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))
        return GPP + Reco
    return LRF

# Eo, rb, k and alpha fixed
def make_LRF_6(Eo, rb, k, alpha):
    def LRF(data_d,Aopt):
        Aopt_VPD = Aopt * np.exp(-k * (data_d['VPD'] - 1))
        index = np.where(data_d['VPD'] <= 1)[0]
        Aopt_VPD[index] = Aopt
        GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + 
               (alpha * data_d['PAR'] / Aopt_VPD))
        index = np.where(data_d['PAR'] < 10)[0]
        GPP[index] = 0
        Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))
        return GPP + Reco
    return LRF
    
#------------------------------------------------------------------------------
# Daytime rb, alpha, Aopt, k

def optimise_with_rb(data_dict, params_dict):
    """
    This script simultaneously finds optimal parameter values of: 
        i)  a rectangular hyperbolic light response function of the form:
            GPP = (alpha * PAR) / (1 - (PAR / 2000) + (alpha * PAR / Aopt))
            where alpha is the initial slope of the light response curve and Aopt is
            the magnitude of GPP at 2000umol photons m^-2 s^-1
        ii) the reference respiration parameter of the lloyd and Taylor function:
            Re = rb * e^(Eo * (1 / (10 + 46.02) - 1 / (T + 46.02)))
    Positional arguments to pass in are: 
        1) a dictionary of numpy arrays (data_dict, np.nan as missing data 
           value) with the following keys:
               - 'NEE' - carbon flux in umol m^-2 s^-1
               - 'TempC' - temperature in celsius
               - 'PAR' - photosynthetically active radiation in umol photons m^-2 s^-1 
               - 'VPD' - vapour pressure deficit in kPa 
        2) a dictionary containing required configurations (configs_dict), with the 
           following keys:
               'Eo_default' - used to fix value of Eo
               'alpha_default' - used if optimisation with free alpha fails
               'k_default' - used if optimisation with free k fails
               'rb_prior' - initial guess for rb
               'alpha_prior' - initial guess for alpha
               'Aopt_prior' - initial guess for Aopt
               'k_prior' - intial guess for k
    """

    # Initialise error state variable
    error_state = 0        

    # Get drivers and response
    drivers_d = {driver: data_dict[driver] for driver in ['PAR', 'TempC', 'VPD']}
    response_array = data_dict['NEE']      
    
    try:
        params = curve_fit(make_LRF_1(params_dict['Eo_default']), 
                           drivers_d, response_array, 
                           p0 = [params_dict['rb_prior'], 
                                 params_dict['alpha_prior'], 
                                 params_dict['Aopt_prior'], 
                                 params_dict['k_prior']])[0] 
    except RuntimeError:
        params = [np.nan, np.nan, np.nan, np.nan]
    rb_day, alpha, Aopt, k = params[0], params[1], params[2], params[3]
           
    # If nan returned or a negative VPD coefficient, rerun with VPD set to zero
    if ((k == np.nan) | (k < 0)):
        error_state = 1            
        k = params_dict['k_default']
        try:
            params = curve_fit(make_LRF_2(params_dict['Eo_default'], 
                                          params_dict['k_default']), 
                               drivers_d, response_array, 
                               p0 = [params_dict['rb_prior'], 
                                     params_dict['alpha_prior'], 
                                     params_dict['Aopt_prior']])[0]
        except RuntimeError:
            params = [np.nan, np.nan, np.nan]
        rb_day, alpha, Aopt = params[0], params[1], params[2]

    # If a nan, positive or otherwise out of range value of alpha was returned,
    # rerun with previous value of alpha if available, otherwise set to zero
    if ((alpha == np.nan) | (alpha[1] > 0) | (alpha[1] < -0.22)):
        error_state = 2   
        k = params_dict['k_default']
        alpha = params_dict['alpha_default']
        try:            
            params = curve_fit(make_LRF_3(params_dict['Eo_default'], 
                                          params_dict['k_default'], 
                                          params_dict['alpha_default']), 
                               drivers_d, response_array, 
                               p0 = [params_dict['rb_prior'], 
                                     params_dict['Aopt_prior']])[0]
        except RuntimeError:
            error_state = 3
            params = [np.nan, np.nan]
        rb_day, Aopt = params[0], params[1]
    
    # If Aopt or rb is of the wrong sign, reject all parameters
    if Aopt > 0 or rb_day < 0:
        if Aopt > 0 and rb_day < 0:
            error_state = 4
        elif Aopt > 0:
            error_state = 5
        else:
            error_state = 6
        params = [np.nan, np.nan, np.nan, np.nan]
        rb_day, alpha, Aopt, k = params[0], params[1], params[2], params[3]
            
    params = np.array([rb_day, alpha, Aopt, k])
    
    return params, error_state


def optimise_no_rb(data_dict, params_dict): 
    """
    This script finds optimal parameter values of a rectangular hyperbolic light 
    response function of the form:
        GPP = (alpha * PAR) / (1 - (PAR / 2000) + (alpha * PAR / Aopt))
        where alpha is the initial slope of the light response curve and Aopt is
        the magnitude of GPP at 2000umol photons m^-2 s^-1
    Positional arguments to pass in are: 
        1) a dictionary of numpy arrays (data_dict, np.nan as missing data 
           value) with the following keys:
               - 'NEE' - carbon flux in umol m^-2 s^-1
               - 'TempC' - temperature in celsius
               - 'PAR' - photosynthetically active radiation in umol photons m^-2 s^-1 
               - 'VPD' - vapour pressure deficit in kPa 
        2) a dictionary containing required configurations (configs_dict), with the 
           following keys:
               'Eo_default' - used to fix value of Eo
               'alpha_default' - used if optimisation with free alpha fails
               'k_default' - used if optimisation with free k fails
               'alpha_prior' - initial guess for alpha
               'Aopt_prior' - initial guess for Aopt
               'k_prior' - initial guess for k
    """
    
    # Initialise error state variable
    error_state = 0        

    # Get drivers and response
    drivers_d = {driver: data_dict[driver] for driver in ['PAR', 'TempC', 'VPD']}
    response_array = data_dict['NEE']      
    
    try:
        params = curve_fit(make_LRF_4(params_dict['Eo_default'], 
                                      params_dict['rb_default']), 
                           drivers_d, response_array, 
                           p0 = [params_dict['alpha_prior'], 
                                 params_dict['Aopt_prior'], 
                                 params_dict['k_prior']])[0]

    except RuntimeError:
        params = [np.nan, np.nan, np.nan]
    alpha, Aopt, k = params[0], params[1], params[2]
           
    # If nan returned or a negative VPD coefficient, rerun with VPD set to zero
    if ((k == np.nan) | (k < 0)):
        error_state = 1            
        k = params_dict['k_default']
        try:
            params = curve_fit(make_LRF_5(params_dict['Eo_default'], 
                                          params_dict['rb_default'], 
                                          params_dict['k_default']), 
                               drivers_d, response_array, 
                               p0 = [params_dict['alpha_prior'], 
                                     params_dict['Aopt_prior']])[0]
        except RuntimeError:
            params = [np.nan, np.nan]
        alpha, Aopt = params[0], params[1]

    # If a nan, positive or otherwise out of range value of alpha was returned,
    # rerun with previous value of alpha if available, otherwise set to zero
    if ((alpha == np.nan) | (alpha > 0) | (alpha < -0.22)):
        error_state = 2   
        k = params_dict['k_default']
        alpha = params_dict['alpha_default']
        try:            
            params = curve_fit(make_LRF_6(params_dict['Eo_default'],
                                          params_dict['rb_default'],
                                          params_dict['k_default'], 
                                          params_dict['alpha_default']), 
                               drivers_d, response_array, 
                               p0 = [params_dict['Aopt_prior']])[0]
        except RuntimeError:
            error_state = 3
            params = [np.nan]
        Aopt = params[0]
        
    params = np.array([alpha, Aopt, k])
    
    return params, error_state