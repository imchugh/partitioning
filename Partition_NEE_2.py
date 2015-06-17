# -*- coding: utf-8 -*-
"""
Created on Wed May 27 09:28:52 2015

@author: imchugh
"""
import Tkinter, tkFileDialog
from configobj import ConfigObj
import numpy as np
from scipy.optimize import curve_fit
import datetime as dt
import os
import sys
import pdb

#------------------------------------------------------------------------------
# Data optimisation algorithms
    
def TRF(data_d,Eo,rb):
    return rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))

def make_TRF(Eo):
    def TRF(data_d,rb):
        return rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))
    return TRF

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
# Data optimisation procedures
    
def optimise_dark(sub_d, default_params_d, prior_params_d, options_d):

    # Get dark subset
    noct_flag = True
    drivers_d, response, pct = subset_daynight(sub_d, ['TempC'], 'NEE', noct_flag)
    
    # If minimum data criterion satisfied, fit L&T parameters
    if pct > options_d['minimum_pct_noct_window']:

        # Initialise error state variable
        error_state = 0              
        
        try:
            params = curve_fit(make_TRF(default_params_d['Eo']), drivers_d, response, 
                               p0 = [prior_params_d['rb']])[0]
        except RuntimeError:
            params = [np.nan]
    
        # If negative rb returned, set to nan
        if params[0] < 0:
            error_state = 6
            params = [np.nan]

    # If minimum data criterion not satisfied 
    else:
        error_state = 7       
        params = [np.nan]
        
    return params, error_state

# Nighttime for establishment of annual Eo
def optimise_dark_annual(sub_d, prior_params_d, options_d):
    
    # Get dark subset
    noct_flag = True
    drivers_d, response, pct = subset_daynight(sub_d, ['TempC'], 'NEE', noct_flag)
    
    # If minimum data criterion satisfied, fit L&T parameters
    if pct > options_d['minimum_pct_annual']:
        try:
            params = curve_fit(TRF, drivers_d, response, 
                               p0 = [prior_params_d['Eo'], 
                                     prior_params_d['rb']])[0]
        except RuntimeError:
            params = [np.nan, np.nan]
    
    # If minimum data criterion not satisfied 
    else:
        
        params = [np.nan, np.nan]
    
    return params

# Daytime
def optimise_light(sub_d, default_params_d, prior_params_d, options_d):
        
    # Get light subset 
    noct_flag = False
    drivers_d, response, pct = subset_daynight(sub_d, ['PAR', 'TempC', 'VPD'], 
                                               'NEE', noct_flag)
    
    # If minimum data criterion satisfied, fit light response and L&T parameters
    if pct > options_d['minimum_pct_day_window']:
        
        # Initialise error state variable
        error_state = 0        
        
        try:
            params = curve_fit(make_LRF_1(default_params_d['Eo']), 
                               drivers_d, response, 
                               p0 = [prior_params_d['rb'], 
                                     prior_params_d['alpha'], 
                                     prior_params_d['Aopt'], 
                                     prior_params_d['k']])[0] 
        except RuntimeError:
            params = [np.nan, np.nan, np.nan, np.nan]
        rb_day, alpha, Aopt, k = params[0], params[1], params[2], params[3]
               
        # If nan returned or a negative VPD coefficient, rerun with VPD set to zero
        if ((params[3] == np.nan) | (params[3] < 0)):
            error_state = 1            
            k = 0
            try:
                params = curve_fit(make_LRF_2(default_params_d['Eo'], default_params_d['k']), 
                                   drivers_d, response, 
                                   p0 = [prior_params_d['rb'], 
                                         prior_params_d['alpha'], 
                                         prior_params_d['Aopt']])[0]
            except RuntimeError:
                params = [np.nan, np.nan, np.nan]
            rb_day, alpha, Aopt = params[0], params[1], params[2]
    
        # If a nan, positive or otherwise out of range value of alpha was returned,
        # rerun with previous value of alpha if available, otherwise set to zero
        if ((params[1] == np.nan) | (params[1] > 0) | (params[1] < -0.22)):
            error_state = 2   
            k = 0
            alpha = default_params_d['alpha']
            try:            
                params = curve_fit(make_LRF_3(default_params_d['Eo'], default_params_d['k'], 
                                              default_params_d['alpha']), 
                                   drivers_d, response, 
                                   p0 = [prior_params_d['rb'], prior_params_d['Aopt']])[0]
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
        default_params_d['alpha'] = alpha
    else:
        default_params_d['alpha'] = 0    

    return np.array(params), error_state

#------------------------------------------------------------------------------
# Subsetting of day and night (remove all records with bad data in ANY variable)
def subset_daynight(data_d, drivers, NEE_name, noct_flag):
    
    # Turn dictionary into an array
    if not type(drivers) == list:
        drivers = [drivers]
    temp_array = np.empty([len(data_d[NEE_name]), len(drivers) + 1])
    for i, var in enumerate(drivers):
        temp_array[:, i] = data_d[var]
    temp_array[:, -1] = data_d[NEE_name]

    # Create night / day subsetting index, subset data and count records
    if noct_flag:
        daynight_index = np.where(data_d['PAR'] < 10)[0]
    else:
        daynight_index = np.where(data_d['PAR'] > 10)[0]
    temp_array = temp_array[daynight_index]
    num_records = len(temp_array)
    
    # Create nan subsetting index, subset data and count
    QCdata_index = np.where(np.all(~np.isnan(temp_array), axis=1))    
    temp_array = temp_array[QCdata_index]
    valid_records = len(temp_array)
    percent_avail = round(float(valid_records) / num_records * 100, 1)
    
    # Build dictionary for driver, and array for response
    drivers_d = {var: temp_array[:, i] for i, var in enumerate(drivers)}
    response = temp_array[:, -1]
    
    return drivers_d, response, percent_avail

# Subsetting of date window
def subset_window(data_d, date_array, date, options_d):

    # Assign configs to local vars
    window = options_d['window_size_days']    
    
    # Find bracketing dates
    date_time = (dt.datetime.combine(date, dt.datetime.min.time()) 
                 + dt.timedelta(hours = 12))
    start_date = date_time - dt.timedelta(window / 2.0)
    end_date = date_time + dt.timedelta(window / 2.0)
    
    # Get index for right dates, then subset the arrays
    index = np.where((date_array > start_date) & (date_array <= end_date))
    sub_d = {}
    for i in data_d.keys():
        sub_d[i] = data_d[i][index]
    
    return sub_d
    
#------------------------------------------------------------------------------
# Open configuration and build dictionaries of config file contents
def get_configs():
    
    # Prompt user for configuration file and get it
    root = Tkinter.Tk(); root.withdraw()
    cfName = tkFileDialog.askopenfilename(initialdir='')
    root.destroy()
    cf=ConfigObj(cfName)
    
    # Map input file variable names to hard-coded names used in this script
    vars_d = dict(cf['variables'])
    temp_d = {'carbon_flux': 'NEE', 
              'solar_radiation': 'PAR', 
              'temperature': 'TempC', 
              'vapour_pressure_deficit': 'VPD',
              'date_time': 'date_time'}
    vars_d = {temp_d[i]: vars_d['data'][i] for i in temp_d.keys()}
    
    # Prepare dictionary of user settings - drop strings or change to int / float
    options_d = dict(cf['options'])
    for key in options_d:
        if options_d[key].isdigit():
            options_d[key] = int(options_d[key])
        else:
            try:
                options_d[key] = float(options_d[key])
            except ValueError:
                continue
    
    # Set input file and output path
    paths_d = {}    
    paths_d['file_in'] = os.path.join(cf['files']['input_path'],cf['files']['input_file'])
    paths_d['results_output_path'] = os.path.join(cf['files']['output_path'],'Results')
    paths_d['plot_output_path'] = os.path.join(cf['files']['output_path'],'Plots')
        
    return paths_d, vars_d, options_d

def get_data(paths_d, vars_d):
    
    data_arr = np.load(paths_d['file_in'])
    data_d = {varName: data_arr[vars_d[varName]] for varName in vars_d.keys()}
    return data_d

# Create a dictionary with initial guesses for parameters
def make_initial_guess_dict(data_d):

    d = {}
    d['Eo'] = 100
    d['k'] = 0
    d['alpha'] = -0.01
    index = np.where(data_d['PAR'] < 10)[0]
    d['rb'] = data_d['NEE'][index].mean()
    index = np.where(data_d['PAR'] > 10)[0]
    d['Aopt'] = (np.nanpercentile(data_d['NEE'][index], 5) - 
                 np.nanpercentile(data_d['NEE'][index], 95))
    return d

def annual_Eo(data_d, prior_params_d, options_d, date_array):
    
    # Create a list of the number of years
    year_array = np.array([i.year for i in date_array])
    year_list = list(set(year_array))
    
    # Get Eo for each year and compile dictionary
    yearsEo_d = {}
    print 'Eo optimised using whole year is as follows:'
    for yr in year_list:
        index = np.where(year_array == yr)
        sub_d = {}
        for item in data_d.keys():
            sub_d[item] = data_d[item][index]
        params = optimise_dark_annual(sub_d, prior_params_d, options_d)
        yearsEo_d[str(yr)] = params[0]
        print '    - ' + str(yr) + ': ' + str(round(params[0]))
    
    # Do QC on Eo
    Eo_pass_keys = []
    Eo_range_fail_keys = []
    Eo_nan_fail_keys = []
    for yr in yearsEo_d.keys():
        if np.isnan(yearsEo_d[yr]):
            Eo_nan_fail_keys.append(yr)
        elif ((yearsEo_d[yr] < 50) | (yearsEo_d[yr] > 400)):
            Eo_range_fail_keys.append(yr)
        else:
            Eo_pass_keys.append(yr)
    if len(Eo_pass_keys) != len(yearsEo_d):
        if len(Eo_nan_fail_keys) == len(yearsEo_d):
            print 'Could not find any values of Eo for any years! Exiting...'
            sys.exit()
        elif len(Eo_pass_keys) != 0:
            Eo_mean = [yearsEo_d[i] for i in Eo_pass_keys].mean()
            for i in (Eo_range_fail_keys + Eo_nan_fail_keys):
                yearsEo_d[i] = Eo_mean
            print 'Eo optimisation failed for the following years: '
            print [i for i in (Eo_range_fail_keys + Eo_nan_fail_keys)]
            print 'Eo for these years estimated from the mean of all other years'
        else:
            for i in Eo_range_fail_keys:
                yearsEo_d[i] == 50 if yearsEo_d[i] < 50 else 400
            print 'Warning! Eo estimates were out of range for all years'
            print 'Low estimates have been set to lower limit (50);'
            print 'High estimates have been set to upper limit (400);'
            print 'Parameter estimates are unlikely to be robust!'
    else:
        print 'Eo estimates passed QC for all years'
    
    return yearsEo_d

def get_dates(dateStr_array, options_d):
    
    # Assign configs to local vars
    window = options_d['window_size_days']
    step = options_d['step_size_days']
    msmt_interval = options_d['msmt_interval_hrs']

    # Create a datetime array from string array
    # Note that we lag by one measurement interval in minutes from the date series because
    # the last valid case in the 24 hours of data occurs at midnight (assuming that the
    # naming convention is for the timestamp signifying the end of the averaging interval).
    # So, for example, assuming half hourly interval, the first valid case is 0030 
    # (average of 0000-0030) and final valid case is 0000 (average of 2330 to 0000).
    date_array = np.array([dt.datetime.strptime(i, '%Y-%m-%d %H:%M:%S')
                           for i in dateStr_array])    
    date_array = date_array - dt.timedelta(minutes = 60 * msmt_interval)

    # Create a series of continuous whole day dates that will be used for output
    # (parameter series will be interpolated between window centres)
    start_date = date_array[0].date()
    end_date = date_array[-1].date()
    num_days = (end_date - start_date).days + 1 # Add 1 so is inclusive of both end members
    all_dates_array = np.array([start_date + dt.timedelta(i) for i in xrange(num_days)])
    
    # Check that first and last days are complete and revise start and end dates if required
    all_dates = np.array([i.date() for i in date_array])
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
    
    return date_array, all_dates_array, step_dates_array

def make_error_code_dict():
    
    msg_d = {1:'Value of k failed range check - setting to zero and ' \
               'recalculating other parameters',
             2:'Value of alpha failed range check - using previous ' \
               'estimate and recalculating other parameters',
             3:'Optimisation reached maximum number of interations' \
               'without convergence',
             4:'Value of Aopt and rb have wrong sign - ' \
               'rejecting all parameters',
             5:'Value of Aopt has wrong sign - rejecting all parameters',
             6:'Value of daytime rb has wrong sign - rejecting all parameters',
             7:'Value of nocturnal rb has wrong sign - rejecting',
             8:'Data did not pass minimum percentage threshold - ' \
               'skipping optimisation'}
    
    return msg_d

# Main code starts here
def main():
    
    # Create dictionaries from configuration file
    paths_d, vars_d, options_d = get_configs()

    # Get the data and enter into dictionary
    data_d = get_data(paths_d, vars_d)
    
    # Create a dictionary containing initial guesses for each parameter
    prior_params_d = make_initial_guess_dict(data_d)
    
    # Create and initialise a dictionary containing default parameters
    # (used when optimisation fails)
    default_params_d = {'alpha': 0, 'k': 0}
    
    # Create dictionary containing error messages with error codes as keys
    msg_d = make_error_code_dict()
    
    # Get arrays of all datetimes, all dates and stepped dates
    datetime_array, all_dates_array, step_dates_array = get_dates(data_d.pop('date_time'), 
                                                                  options_d)    

    # Initialise result array
    rslt_array = np.empty([len(all_dates_array), 9])
    rslt_array[:,1:] = np.nan    

    # Get the annual estimates of Eo
    yearsEo_d = annual_Eo(data_d, prior_params_d, options_d, datetime_array)
#    annual_Eo(data_d, prior_params_d, options_d, datetime_array, rslt_array)

    # Do optimisation for each window and write to result array
    for date in step_dates_array:
        
        # Subset the data
        sub_d = subset_window(data_d, datetime_array, date, options_d)
        
        # Get Eo for the relevant year and write to the parameters dictionary
        Eo_current_year = yearsEo_d[str(date.year)]
        default_params_d['Eo'] = Eo_current_year
        
        # Get the parameters and write to the results array
        index = np.where(all_dates_array == date)
        rslt_array[index, 0] = Eo_current_year
        dark_rb_param, dark_rb_error_state = optimise_dark(sub_d, default_params_d, 
                                                      prior_params_d, options_d)        
        rslt_array[index, 1] = dark_rb_param
        # Dark Eo error state at [,6]
        rslt_array[index, 7] = dark_rb_error_state
        light_params, light_error_state = optimise_light(sub_d, default_params_d, 
                                                         prior_params_d, options_d)
        rslt_array[index, 2:6] = light_params
        rslt_array[index, 8] = light_error_state
        
        # Print error messages if any
        if dark_rb_error_state != 0 or light_error_state != 0:
            print 'For ' + dt.datetime.strftime(date, '%Y-%m-%d') + ':'
            if dark_rb_error_state != 0:            
                print msg_d[dark_rb_error_state]
            if light_error_state != 0:
                print msg_d[light_error_state]
     
    return rslt_array