# -*- coding: utf-8 -*-
"""
Created on Wed May 27 09:28:52 2015

@author: imchugh
"""
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import datetime as dt
import os
import sys
import pdb

# Main code starts here
def main(data_dict, configs_dict):

    # Assign configs to local vars
    window = configs_dict['window_size_days']

    # If user wants individual window plots, check whether output directories
    # are present, and create if not
    if configs_dict['output_plots']:
        output_path = configs_dict['plot_output_path']
        if not os.path.isdir(output_path): 
            os.makedirs(output_path)    
    
    # Create a dictionary containing initial guesses for each parameter
    prior_params_dict = make_initial_guess_dict(data_dict)
    
    # Create and initialise a dictionary containing default parameters
    # (used when optimisation fails)
    default_params_dict = {'alpha': 0, 'k': 0}
    
    # Create dictionary containing error messages with error codes as keys
    msg_dark_dict = make_error_code_dict(noct_flag = True)
    msg_light_dict = make_error_code_dict(noct_flag = False)
    
    # Get arrays of all datetimes, all dates and stepped dates
    (working_datetime_array, 
     datetime_array, 
     all_dates_array, 
     step_dates_array) = get_dates(data_dict.pop('date_time'), configs_dict)    
    
    # Create variable name lists for results output
    param_rslt_list = ['Eo', 'rb_noct', 'rb_day', 'alpha', 'Aopt', 'k', 
                       'Eo_error_code', 'rb_noct_error_code', 'light_error_code',
                       'noct_fill_code (0 = calculated; 1 = interpolated)',
                       'day_fill_code (0 = calculated; 1 = interpolated)']
    
    series_rslt_list = ['Re_noct', 'Re_day',
                        'GPP_est',
                        'GPP_calc_noct (NEE_obs - Re_noct)', 
                        'GPP_calc_day (NEE_obs - Re_day)',
                        'NEE_obs', 'NEE_est'] 
    
    # Initialise results arrays
    param_rslt_array = np.empty([len(all_dates_array), 11])
    param_rslt_array[:] = np.nan
    series_rslt_array = np.empty([len(working_datetime_array), 7])
    series_rslt_array[:] = np.nan

    # Get the annual estimates of Eo and write to the parameter results array
    (param_rslt_array[:,0], 
     param_rslt_array[:,6]) = optimise_annual_Eo(data_dict, prior_params_dict, 
                                                 configs_dict, 
                                                 working_datetime_array, 
                                                 all_dates_array)    
    
    # Do optimisation for each window
    for date in step_dates_array:

        # Subset the data
        sub_dict, series_index = subset_window(data_dict, working_datetime_array, 
                                               date, window)
        
        # Get Eo for the relevant year and write to the parameters dictionary
        default_params_dict['Eo'] = param_rslt_array[np.where(all_dates_array == date), 0].item()
        
        # Get the parameters and write to the parameter results array
        param_index = np.where(all_dates_array == date)
        dark_rb_param, dark_rb_error_state = optimise_dark(sub_dict, 
                                                           default_params_dict, 
                                                           prior_params_dict, 
                                                           configs_dict)        
        param_rslt_array[param_index, 1] = dark_rb_param
        param_rslt_array[param_index, 7] = dark_rb_error_state
        light_params, light_error_state = optimise_light(sub_dict, 
                                                         default_params_dict, 
                                                         prior_params_dict, 
                                                         configs_dict)
        param_rslt_array[param_index, 2:6] = light_params
        param_rslt_array[param_index, 8] = light_error_state

        # Print error messages if any
        if dark_rb_error_state != 0 or light_error_state != 0:
            print 'For ' + dt.datetime.strftime(date, '%Y-%m-%d') + ':'
            if dark_rb_error_state != 0:            
                print msg_dark_dict[dark_rb_error_state]
            if light_error_state != 0:
                print msg_light_dict[light_error_state]

        # Write current parameters to estimated parameters dictionary
        current_params = param_rslt_array[np.where(all_dates_array == date), :6].reshape(6)
        est_params_d = {param_rslt_list[i]: current_params[i] for i in range(6)}
        
        # Estimate the time series data and plot
        est_series_d = estimate_Re(sub_dict, est_params_d)
        if est_series_d != None and configs_dict['output_plots']:
            combine_d = dict(sub_dict, **est_series_d)
            plot_windows(combine_d, configs_dict, date)

    # Create a binary var indicating calculated or interpolated, then interpolate
    param_rslt_array[:, 9] = np.where(~np.isnan(param_rslt_array[:, 1]), 1, 0)
    param_rslt_array[:, 10] = np.where(~np.isnan(param_rslt_array[:, 2]), 1, 0)
    param_rslt_array[:, :6] = interp_params(param_rslt_array)

    # Loop through interpolated data, construct the time series and do plotting
    for ind, date in enumerate(all_dates_array):
        
        # Subset the data (single day)
        sub_dict, series_index = subset_window(data_dict, working_datetime_array, date, 1)

        # Write current parameters to estimated parameters dictionary
        current_params = param_rslt_array[ind, :6]
        est_params_dict = {param_rslt_list[i]: current_params[i] for i in range(6)}
        
        # Estimate the time series data and write to the results 
        est_series_dict = estimate_Re(sub_dict, est_params_dict)
        if est_series_dict != None:
            for i, var in enumerate(series_rslt_list):
                series_rslt_array[series_index, i] = est_series_dict[var]
    
    # Make dictionaries
    params_dict = {var: param_rslt_array[:, i] for i, var in enumerate(param_rslt_list)}
    params_dict['date_time'] = all_dates_array
    series_dict = {var: series_rslt_array[:, i] for i, var in enumerate(series_rslt_list)}
    series_dict['date_time'] = datetime_array    
    
    return params_dict, series_dict

#------------------------------------------------------------------------------
# Data optimisation algorithms
    
def TRF(data_d,Eo,rb):
    return rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))

def make_TRF(Eo):
    def TRF(data_d,rb):
        return rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))
    return TRF

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

#------------------------------------------------------------------------------        
# Estimate Re (from nocturnal and daytime parameters) and GPP
def estimate_Re(sub_d, params_d):

    if np.any(np.isnan(params_d.values())):
        
        return

    else:
        
        # Estimate GPP and Re using nocturnal rb
        GPP, Re_noct = LRF(sub_d, params_d['Eo'], params_d['rb_noct'],
                           params_d['alpha'], params_d['Aopt'], params_d['k'])
    
        # Estimate GPP and Re using daytime rb
        GPP, Re_day = LRF(sub_d, params_d['Eo'], params_d['rb_day'],
                           params_d['alpha'], params_d['Aopt'], params_d['k'])

        return_d = {'Re_noct': Re_noct, 'Re_day': Re_day,
                    'GPP_est': GPP,
                    'GPP_calc_noct (NEE_obs - Re_noct)': sub_d['NEE'] - Re_noct,
                    'GPP_calc_day (NEE_obs - Re_day)': sub_d['NEE'] - Re_day,
                    'NEE_obs': sub_d['NEE'],     
                    'NEE_est': GPP + Re_day}
        
        return return_d

#------------------------------------------------------------------------------        
# Get dates
def get_dates(date_array, configs_dict):
    
    # Assign configs to local vars
    window = configs_dict['window_size_days']
    step = configs_dict['step_size_days']
    msmt_interval = configs_dict['msmt_interval_hrs']

    # Lag by one measurement interval in minutes from the date series because
    # the last valid case in the 24 hours of data occurs at midnight (assuming that the
    # naming convention is for the timestamp signifying the end of the averaging interval).
    # So, for example, assuming half hourly interval, the first valid case is 0030 
    # (average of 0000-0030) and final valid case is 0000 (average of 2330 to 0000).
    working_date_array = date_array - dt.timedelta(minutes = 60 * msmt_interval)

    # Create a series of continuous whole day dates that will be used for output
    # (parameter series will be interpolated between window centres)
    start_date = working_date_array[0].date()
    end_date = working_date_array[-1].date()
    num_days = (end_date - start_date).days + 1 # Add 1 so is inclusive of both end members
    all_dates_array = np.array([start_date + dt.timedelta(i) for i in xrange(num_days)])
    
    # Check that first and last days are complete and revise start and end dates if required
    all_dates = np.array([i.date() for i in working_date_array])
    recs_required = 24 * (1 / msmt_interval)
    recs_count = 0
    loop_count = 0
    while recs_count != recs_required:
        recs_count = len(np.where(all_dates == all_dates_array[loop_count])[0])
        loop_count = loop_count + 1
    start_date = all_dates_array[loop_count - 1]
    recs_count = 0
    loop_count = -1
    while recs_count != recs_required:
        recs_count = len(np.where(all_dates == all_dates_array[loop_count])[0])
        loop_count = loop_count - 1
    end_date = all_dates_array[loop_count + 1]
    
    # Calculate the dates that represent the centre of the window for each step
    num_days = (end_date - start_date).days + 1 - window # Add 1 so is inclusive of both end members
    first_fit_day = start_date + dt.timedelta(window / 2)
    step_days = np.arange(0, num_days, step)
    step_dates_array = [first_fit_day + dt.timedelta(i) for i in step_days]
    
    return working_date_array, date_array, all_dates_array, step_dates_array

#------------------------------------------------------------------------------
# Interpolate parameters to create complete series
def interp_params(param_rslt_array):
    
    arr = param_rslt_array[:, :6].copy()
    
    xp = np.arange(len(param_rslt_array))
    for i in range(6):
        fp = arr[:, i]
        nan_index = np.isnan(fp)
        fp[nan_index] = np.interp(xp[nan_index], xp[~nan_index], fp[~nan_index])
        arr[:, i] = fp
    return arr

#------------------------------------------------------------------------------
# Make dictionaries

# Write error messages to dictionary with codes as keys
def make_error_code_dict(noct_flag):
    
    if noct_flag:
        
        d = {1:'Value of nocturnal rb has wrong sign - rejecting',
             10:'Data did not pass minimum percentage threshold - ' \
                'skipping optimisation'}
        
    else:
    
        d = {1:'Value of k failed range check - setting to zero and ' \
               'recalculating other parameters',
             2:'Value of alpha failed range check - using previous ' \
               'estimate (zero if unavailable) and recalculating other ' \
               'parameters',
             3:'Optimisation reached maximum number of interations' \
               'without convergence',
             4:'Value of Aopt and rb have wrong sign - ' \
               'rejecting all parameters',
             5:'Value of Aopt has wrong sign - rejecting all parameters',
             6:'Value of daytime rb has wrong sign - rejecting all parameters',
             10:'Data did not pass minimum percentage threshold - ' \
                'skipping optimisation'}
    
    return d

# Create a dictionary with initial guesses for parameters
def make_initial_guess_dict(data_d):

    d = {}
    d['Eo'] = 100
    d['k'] = 0
    d['alpha'] = -0.01
    index = np.where(data_d['PAR'] < 10)[0]
    d['rb'] = np.nanmean(data_d['NEE'][index])
    index = np.where(data_d['PAR'] > 10)[0]
    d['Aopt'] = (np.nanpercentile(data_d['NEE'][index], 5) - 
                 np.nanpercentile(data_d['NEE'][index], 95))
    return d

#------------------------------------------------------------------------------
# Data optimisation procedures

# Annual Eo
def optimise_annual_Eo(data_dict, prior_params_dict, configs_dict, 
                       datetime_array, all_dates_array):
    
    # Initialise local variables with configurations
    min_pct = configs_dict['minimum_pct_annual']    
    
    # Create a list of the number of years
    year_array = np.array([i.year for i in datetime_array])
    year_list = list(set(year_array))
    
    # Get Eo for each year and compile dictionary
    yearsEo_d = {}
    Eo_pass_keys = []
    Eo_range_fail_keys = []
    Eo_nan_fail_keys = []
    print 'Eo optimised using whole year is as follows:'
    for yr in year_list:

        # Subset data        
        index = np.where(year_array == yr)
        year_d = {}
        for item in data_dict.keys():
            year_d[item] = data_dict[item][index]
        noct_flag = True
        sub_d, pct = subset_daynight(year_d, noct_flag)
        
        # Fit L&T parameters if minimum data criterion satisfied, otherwise nan
        if pct > min_pct:
            drivers_d = {driver: sub_d[driver] for driver in ['TempC']}
            response_array = sub_d['NEE']
            try:
                params = curve_fit(TRF, drivers_d, response_array, 
                                   p0 = [prior_params_dict['Eo'], 
                                         prior_params_dict['rb']])[0]
            except RuntimeError:
                params = [np.nan, np.nan]
        else:
            params = [np.nan, np.nan]                                           

        # Assign year to pass, range_fail or nan_fail list for subsequent QC and fill
        Eo = params[0]
        yearsEo_d[yr] = Eo
        if np.isnan(Eo):
            Eo_nan_fail_keys.append(yr)
        elif ((Eo < 50) | (Eo > 400)):
            Eo_range_fail_keys.append(yr)
        else:
            Eo_pass_keys.append(yr)

        print '    - ' + str(yr) + ': ' + str(round(params[0]))
    
    # Do QC on Eo
    yearsQC_d = {yr: 0 for yr in year_list}
    if len(Eo_pass_keys) != len(yearsEo_d):
        if len(Eo_nan_fail_keys) == len(yearsEo_d):
            print 'Could not find any values of Eo for any years! Exiting...'
            sys.exit()
        elif len(Eo_pass_keys) != 0:
            try:
                Eo_mean = np.array([yearsEo_d[i] for i in Eo_pass_keys]).mean()
            except KeyError:
                pdb.set_trace()
            for i in (Eo_range_fail_keys + Eo_nan_fail_keys):
                yearsEo_d[i] = Eo_mean
                yearsQC_d[i] = 1
            print 'Eo optimisation failed for the following years: '
            print [i for i in (Eo_range_fail_keys + Eo_nan_fail_keys)]
            print 'Eo for these years estimated from the mean of all other years'
        else:
            for i in Eo_range_fail_keys:
                yearsEo_d[i] == 50 if yearsEo_d[i] < 50 else 400
                yearsQC_d[i] = 2
            Eo_mean = [yearsEo_d[i] for i in Eo_range_fail_keys].mean()
            for i in Eo_nan_fail_keys:
                yearsEo_d[i] = Eo_mean
                yearsQC_d[i] = 3
            print 'Warning! Eo estimates were out of range for all years'
            print 'Low estimates have been set to lower limit (50);'
            print 'High estimates have been set to upper limit (400);'
            print 'Parameter estimates are unlikely to be robust!'
    else:
        print 'Eo estimates passed QC for all years'
    
    # Write to arrays of same length as result arrays
    Eo_array = np.empty([len(all_dates_array)])
    QC_array = np.empty([len(all_dates_array)])
    year_array = np.array([i.year for i in all_dates_array])
    for yr in year_list:
        index = np.where(year_array == yr)
        Eo_array[index] = yearsEo_d[yr]
        QC_array[index] = yearsQC_d[yr]
    
    return Eo_array, QC_array

# Night rb    
def optimise_dark(data_dict, default_params_dict, prior_params_dict, configs_dict):

    # Initialise local variables with configurations
    min_pct = configs_dict['minimum_pct_noct_window']

    # Get dark subset
    noct_flag = True
    sub_dict, pct = subset_daynight(data_dict, noct_flag)
    
    # If minimum data criterion satisfied, fit L&T parameters
    if pct > min_pct:

        # Initialise error state variable
        error_state = 0              
        
        # Get drivers and response
        drivers_d = {driver: sub_dict[driver] for driver in ['TempC']}
        response_array = sub_dict['NEE']        
        
        try:
            params = curve_fit(make_TRF(default_params_dict['Eo']), 
                               drivers_d, response_array, 
                               p0 = [prior_params_dict['rb']])[0]
        except RuntimeError:
            params = [np.nan]
    
        # If negative rb returned, set to nan
        if params[0] < 0:
            error_state = 1
            params = [np.nan]

    # If minimum data criterion not satisfied 
    else:
        error_state = 10       
        params = [np.nan]
        
    return params, error_state

# Daytime rb, alpha, Aopt, k
def optimise_light(data_dict, default_params_dict, prior_params_dict, configs_dict):

    # Initialise local variables with configurations
    min_pct = configs_dict['minimum_pct_day_window']
        
    # Get light subset 
    noct_flag = False
    sub_d, pct = subset_daynight(data_dict, noct_flag)
    
    # If minimum data criterion satisfied, fit light response and L&T parameters
    if pct > min_pct:
        
        # Initialise error state variable
        error_state = 0        

        # Get drivers and response
        drivers_d = {driver: sub_d[driver] for driver in ['PAR', 'TempC', 'VPD']}
        response_array = sub_d['NEE']      
        
        try:
            params = curve_fit(make_LRF_1(default_params_dict['Eo']), 
                               drivers_d, response_array, 
                               p0 = [prior_params_dict['rb'], 
                                     prior_params_dict['alpha'], 
                                     prior_params_dict['Aopt'], 
                                     prior_params_dict['k']])[0] 
        except RuntimeError:
            params = [np.nan, np.nan, np.nan, np.nan]
        rb_day, alpha, Aopt, k = params[0], params[1], params[2], params[3]
               
        # If nan returned or a negative VPD coefficient, rerun with VPD set to zero
        if ((params[3] == np.nan) | (params[3] < 0)):
            error_state = 1            
            k = 0
            try:
                params = curve_fit(make_LRF_2(default_params_dict['Eo'], default_params_dict['k']), 
                                   drivers_d, response_array, 
                                   p0 = [prior_params_dict['rb'], 
                                         prior_params_dict['alpha'], 
                                         prior_params_dict['Aopt']])[0]
            except RuntimeError:
                params = [np.nan, np.nan, np.nan]
            rb_day, alpha, Aopt = params[0], params[1], params[2]
    
        # If a nan, positive or otherwise out of range value of alpha was returned,
        # rerun with previous value of alpha if available, otherwise set to zero
        if ((params[1] == np.nan) | (params[1] > 0) | (params[1] < -0.22)):
            error_state = 2   
            k = 0
            alpha = default_params_dict['alpha']
            try:            
                params = curve_fit(make_LRF_3(default_params_dict['Eo'], default_params_dict['k'], 
                                              default_params_dict['alpha']), 
                                   drivers_d, response_array, 
                                   p0 = [prior_params_dict['rb'], prior_params_dict['Aopt']])[0]
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
        
    # If minimum data criterion not satisfied    
    else:
        
        error_state = 8
        rb_day, alpha, Aopt, k = np.nan, np.nan, np.nan, np.nan
    
    params = [rb_day, alpha, Aopt, k]        
    
    # If a valid estimate of alpha is found, write to default params dictionary
    if not np.isnan(alpha):
        default_params_dict['alpha'] = alpha
    else:
        default_params_dict['alpha'] = 0    

    return np.array(params), error_state

#------------------------------------------------------------------------------
# Plots 

# Daytime and nocturnal fits for each window
def plot_windows(data_d, paths_d, options_d, date):

    # Set parameters from dicts
    path = paths_d['plot_output_path']
    window = options_d['window_size_days']
    
    for i in range(2):
        noct_flag = i == False
        sub_d = subset_daynight(data_d, noct_flag)[0]
        if noct_flag:
            daynight_ind = 'noct'
            x_lab = r'Temperature ($^{o}C$)'
            x_var = sub_d['TempC']
            y_var1 = sub_d['NEE']
            y_var2 = sub_d['Re_noct']
        else:            
            daynight_ind = 'day'
            x_lab = r'PAR ($\mu mol\/photons\/m^{-2}s^{-1}$)'
            x_var = sub_d['PAR']
            y_var1 = sub_d['NEE']
            y_var2 = sub_d['NEE_est']
              
        # Plot
        date_str = dt.datetime.strftime(date,'%Y-%m-%d')
        fig = plt.figure(figsize = (12,8))
        fig.patch.set_facecolor('white')
        plt.plot(x_var, y_var1, 'bo' , label = 'NEE_obs')
        plt.plot(x_var, y_var2, 'ro', label = 'NEE_est')
        plt.title('Fit for ' + str(window) + ' day window centred on ' + 
                  date_str + '\n', fontsize = 22)
        plt.xlabel(x_lab, fontsize = 16)
        plt.ylabel(r'Fc ($\mu mol C\/m^{-2} s^{-1}$)', fontsize = 16)
        plt.axhline(y = 0, color = 'black')
        plot_out_name = daynight_ind + '_' + date_str + '.jpg'
        plt.tight_layout()
        fig.savefig(os.path.join(path, plot_out_name))
        plt.close(fig)
        
    return

#------------------------------------------------------------------------------
# Subsetting of day and night (remove all records with bad data in ANY variable)
def subset_daynight(data_dict, noct_flag):
    
    # Turn dictionary into an array
    temp_array = np.empty([len(data_dict['NEE']), len(data_dict)])
    for i, var in enumerate(data_dict.keys()):
        temp_array[:, i] = data_dict[var]

    # Create night / day subsetting index, subset data and count records
    if noct_flag:
        daynight_index = np.where(data_dict['PAR'] < 10)[0]
    else:
        daynight_index = np.where(data_dict['PAR'] > 10)[0]
    temp_array = temp_array[daynight_index]
    num_records = len(temp_array)
    
    # Create nan subsetting index, subset data and count
    QCdata_index = np.where(np.all(~np.isnan(temp_array), axis=1))    
    temp_array = temp_array[QCdata_index]
    valid_records = len(temp_array)
    percent_avail = round(float(valid_records) / num_records * 100, 1)

    sub_dict = {var: temp_array[:, i] for i, var in enumerate(data_dict.keys())}

    return sub_dict, percent_avail

# Subsetting of date window
def subset_window(data_dict, date_array, date, window):

    # Find bracketing dates
    date_time = (dt.datetime.combine(date, dt.datetime.min.time()) 
                 + dt.timedelta(hours = 12))
    start_date = date_time - dt.timedelta(window / 2.0)
    end_date = date_time + dt.timedelta(window / 2.0)
    
    # Get index for right dates, then subset the arrays
    index = np.where((date_array >= start_date) & (date_array < end_date))
    sub_dict = {}
    for i in data_dict.keys():
        sub_dict[i] = data_dict[i][index]
    
    return sub_dict, index

#------------------------------------------------------------------------------