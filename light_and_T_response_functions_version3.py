# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 10:54:05 2015

@author: imchugh
"""
from scipy.optimize import curve_fit
import numpy as np
import pdb

#------------------------------------------------------------------------------
# Data optimisation algorithm

def LRF(data_d, Eo, rb, alpha, beta, k):
    beta_VPD = beta * np.exp(-k * (data_d['VPD'] - 1))
    index = data_d['VPD'] <= 1
    beta_VPD[index] = beta
    GPP = (alpha * beta * data_d['PAR']) / (alpha * data_d['PAR'] + beta)    
    index = data_d['PAR'] < 10
    GPP[index] = 0
    Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))
    return GPP + Reco

#------------------------------------------------------------------------------
# Daytime rb, alpha, beta, k

def optimise_with_rb(data_dict, params_dict):
    """
    This script simultaneously finds optimal parameter values of: 
        i)  a rectangular hyperbolic light response function of the form:
            GPP = (alpha * PAR) / (1 - (PAR / 2000) + (alpha * PAR / beta))
            where alpha is the initial slope of the light response curve and beta is
            the magnitude of GPP at 2000umol photons m^-2 s^-1
        ii) the reference respiration parameter of the lloyd and Taylor function:
            Re = rb * e^(Eo * (1 / (10 + 46.02) - 1 / (T + 46.02)))
    Positional arguments to pass in are: 
        1) a dictionary of numpy arrays (data_dict) with the following keys:
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
               'beta_prior' - initial guess for beta
               'k_prior' - intial guess for k
    Note - no filtering is done here - any missing data values or nans will 
           propagate!!!               
    """

    # Initialise error state variable
    error_state = 0        

    # Get drivers and response
    drivers_d = {driver: data_dict[driver] for driver in ['PAR', 'TempC', 'VPD']}
    response_array = data_dict['NEE']      
    
    # Prepare the beta variations for iteration
    beta_dict = {'beta_1': params_dict['beta_prior'], 
                 'beta_2': 0.5 * params_dict['beta_prior'], 
                 'beta_3': 2 * params_dict['beta_prior']}
    beta_params_dict = {}
    beta_RMSE_dict = {}

    # Go!
    for key in beta_dict.keys():
        
        # Set beta
        params_dict['beta_prior'] = beta_dict[key]

        # Do the initial fit with largest number of parameters
        try:
            params = curve_fit(lambda x, b, c, d, e: 
                               LRF(x,
                                   params_dict['Eo_default'], 
                                   b, c, d, e),
                               drivers_d, 
                               response_array, 
                               p0 = [params_dict['rb_prior'], 
                                     params_dict['alpha_prior'], 
                                     params_dict['beta_prior'], 
                                     params_dict['k_prior']])[0] 
        except RuntimeError:
            params = [np.nan, np.nan, np.nan, np.nan]
        rb, alpha, beta, k = params[0], params[1], params[2], params[3]
               
        # If nan returned or a negative VPD coefficient, rerun with VPD set to zero
        if ((np.isnan(k)) | (k < 0)):
            error_state = 1            
            k = params_dict['k_default']
            try:
                params = curve_fit(lambda x, b, c, d:
                                   LRF(x,
                                       params_dict['Eo_default'],
                                       b, c, d,  
                                       params_dict['k_default']), 
                                   drivers_d, 
                                   response_array, 
                                   p0 = [params_dict['rb_prior'], 
                                         params_dict['alpha_prior'], 
                                         params_dict['beta_prior']])[0]
            except RuntimeError:
                params = [np.nan, np.nan, np.nan]
            rb, alpha, beta = params[0], params[1], params[2]
    
        # If a nan, positive or otherwise out of range value of alpha was returned,
        # rerun with previous value of alpha if available, otherwise set to zero
        if ((np.isnan(alpha)) | (alpha > 0) | (alpha < -0.22)):
            error_state = 2   
            k = params_dict['k_default']
            alpha = params_dict['alpha_default']
            try:            
                params = curve_fit(lambda x, b, d:
                                   LRF(x, 
                                       params_dict['Eo_default'],
                                       b,
                                       params_dict['alpha_default'],
                                       d, 
                                       params_dict['k_default']), 
                                   drivers_d, 
                                   response_array, 
                                   p0 = [params_dict['rb_prior'], 
                                         params_dict['beta_prior']])[0]
            except RuntimeError:
                error_state = 3
                params = [np.nan, np.nan]
            rb, beta = params[0], params[1]
    
        # If beta or rb is of the wrong sign, reject all parameters
        if beta > 0 or rb < 0:
            if beta > 0 and rb < 0:
                error_state = 4
            elif beta > 0:
                error_state = 5
            else:
                error_state = 6
            params = [np.nan, np.nan, np.nan, np.nan]
            rb, alpha, beta, k = params[0], params[1], params[2], params[3]
        
        # If valid beta was obtained, calculate RMSE
        if not np.isnan(beta):
            
            beta_params_dict[key] = params
            GPP_est, Reco_est = LRF(drivers_d, 
                                    params_dict['Eo_default'], 
                                    rb, alpha, beta, k)
            NEE_est = GPP_est + Reco_est
            beta_RMSE_dict[key] = np.sqrt(((response_array - NEE_est[0])**2).mean())
            
    # If the beta RMSE dictionary is not empty, find the lowest value, 
    # get its key and then get the corresponding parameters
    for key in beta_RMSE_dict:
        if beta_RMSE_dict[key] == min(beta_RMSE_dict.values()):
            rb, alpha, beta, k = (beta_params_dict[key][0], 
                                  beta_params_dict[key][1], 
                                  beta_params_dict[key][2], 
                                  beta_params_dict[key][3])
    
    params = np.array([rb, alpha, beta, k])
    
    return params, error_state


def optimise_no_rb(data_dict, params_dict): 
    """
    This script finds optimal parameter values of a rectangular hyperbolic light 
    response function of the form:
        GPP = (alpha * PAR) / (1 - (PAR / 2000) + (alpha * PAR / beta))
        where alpha is the initial slope of the light response curve and beta is
        the magnitude of GPP at 2000umol photons m^-2 s^-1
    Positional arguments to pass in are: 
        1) a dictionary of numpy arrays (data_dict) with the following keys:
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
               'beta_prior' - initial guess for beta
               'k_prior' - initial guess for k
    Note - no filtering is done here - any missing data values or nans will 
           propagate!!!
    """
    
    # Initialise error state variable
    error_state = 0        

    # Get drivers and response
    drivers_d = {driver: data_dict[driver] for driver in ['PAR', 'TempC', 'VPD']}
    response_array = data_dict['NEE']      

    # Prepare the beta variations for iteration
    beta_dict = {'beta_1': params_dict['beta_prior'], 
                 'beta_2': 0.5 * params_dict['beta_prior'], 
                 'beta_3': 2 * params_dict['beta_prior']}
    beta_params_dict = {}
    beta_RMSE_dict = {}

    for key in beta_dict.keys():

        # Set beta
        params_dict['beta_prior'] = beta_dict[key]

        # Do the initial fit with largest number of parameters        
        try:
            params = curve_fit(lambda x, c, d, e:
                               LRF(x,
                                   params_dict['Eo_default'], 
                                   params_dict['rb_default'],
                                   c, d, e), 
                               drivers_d, 
                               response_array, 
                               p0 = [params_dict['alpha_prior'], 
                                     params_dict['beta_prior'], 
                                     params_dict['k_prior']])[0]
        except RuntimeError:
            params = [np.nan, np.nan, np.nan]
        alpha, beta, k = params[0], params[1], params[2]
               
        # If nan returned or a negative VPD coefficient, rerun with VPD set to zero
        if ((k == np.nan) | (k < 0)):
            error_state = 1            
            k = params_dict['k_default']
            try:
                params = curve_fit(lambda x, c, d: 
                                   LRF(x,
                                       params_dict['Eo_default'], 
                                       params_dict['rb_default'],
                                       c, d, 
                                       params_dict['k_default']), 
                                   drivers_d, 
                                   response_array, 
                                   p0 = [params_dict['alpha_prior'], 
                                         params_dict['beta_prior']])[0]
            except RuntimeError:
                params = [np.nan, np.nan]
            alpha, beta = params[0], params[1]
    
        # If a nan, positive or otherwise out of range value of alpha was returned,
        # rerun with previous value of alpha if available, otherwise set to zero
        if ((alpha == np.nan) | (alpha > 0) | (alpha < -0.22)):
            error_state = 2   
            k = params_dict['k_default']
            alpha = params_dict['alpha_default']
            try:            
                params = curve_fit(lambda x, d:
                                   LRF(x,
                                       params_dict['Eo_default'],
                                       params_dict['rb_default'],
                                       params_dict['alpha_default'],
                                       d, 
                                       params_dict['k_default']),
                                   drivers_d, 
                                   response_array, 
                                   p0 = [params_dict['beta_prior']])[0]
            except RuntimeError:
                error_state = 3
                params = [np.nan]
            beta = params[0]
        
        # If valid beta was obtained, calculate RMSE
        if not np.isnan(beta):
            
            beta_params_dict[key] = [alpha, beta, k]
            NEE_est = LRF(drivers_d, 
                          params_dict['Eo_default'], params_dict['rb_default'], 
                          alpha, beta, k)
            beta_RMSE_dict[key] = np.sqrt(((response_array - NEE_est[0])**2).mean())

    # If the beta RMSE dictionary is not empty, find the lowest value, 
    # get its key and then get the corresponding parameters
    for key in beta_RMSE_dict:
        if beta_RMSE_dict[key] == min(beta_RMSE_dict.values()):
            alpha, beta, k = (beta_params_dict[key][0], 
                              beta_params_dict[key][1], 
                              beta_params_dict[key][2])
            break
        
    params = np.array([alpha, beta, k])
    
    print params    
    
    return params, error_state