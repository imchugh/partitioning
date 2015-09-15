# -*- coding: utf-8 -*-
"""
Created on Wed May 27 09:28:52 2015

@author: imchugh
"""
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import os
import sys
import calendar
import pdb

import dark_T_response_functions as dark
import light_and_T_response_functions as light

reload(dark)
reload(light)

# Main code starts here
def main(data_dict, configs_dict):

    # If user wants individual window plots, check whether output directories
    # are present, and create if not
    if configs_dict['configs']['output_plots']:
        output_path = configs_dict['configs']['plot_output_path']
        if not os.path.isdir(output_path): 
            os.makedirs(output_path)

    # Create a dictionary containing initial guesses for each parameter
    params_dict = make_initial_guess_dict(data_dict)
    
    # Create and initialise a dictionary containing default parameters
    # (used when optimisation fails)
    params_dict['alpha_default'] = 0
    params_dict['k_default'] = 0
    
    # Create dictionary containing error messages with error codes as keys
    msg_dark_dict = make_error_code_dict(noct_flag = True)
    msg_light_dict = make_error_code_dict(noct_flag = False)
    
    # Get arrays of all datetimes, all dates and stepped dates
    datetime_array = data_dict.pop('date_time')
    (step_date_index_dict, 
     all_date_index_dict,
     year_index_dict) = get_dates(datetime_array, configs_dict)
    pdb.set_trace()
    date_array = np.array(all_date_index_dict.keys())
    date_array.sort()
    step_date_array = np.array(step_date_index_dict.keys())
    step_date_array.sort()

    # Create variable name lists for results output
    param_rslt_list = ['Eo', 'rb', 'alpha', 'Aopt', 'k']
    error_code_rslt_list = ['Eo_noct', 'rb_noct', 'rb_noct_day', 'rb_day_day']    
       
    series_rslt_list = ['Re_noct', 'GPP_w_Re_noct', 'Re_day', 'GPP_w_Re_day']

    # Initialise results arrays
    param_dk_rslt_array = np.empty([len(date_array), 5])
    param_dk_rslt_array[:] = np.nan
    param_lt_rslt_array = np.empty([len(date_array), 5])
    param_lt_rslt_array[:] = np.nan
    param_QC_array = np.empty([len(date_array), 4])
    param_QC_array[:] = np.nan
    series_rslt_array = np.empty([len(datetime_array), 4])
    series_rslt_array[:] = np.nan

    # Get the annual estimates of Eo
    Eo_dict, EoQC_dict = optimise_annual_Eo(data_dict, 
                                            params_dict, 
                                            configs_dict,
                                            year_index_dict)

    # Write to result arrays
    year_array = np.array([i.year for i in date_array])
    for yr in year_array:
        index = np.where(year_array == yr)
        param_dk_rslt_array[index, 0] = Eo_dict[yr]
        param_lt_rslt_array[index, 0] = Eo_dict[yr]
        param_QC_array[index, 0] = EoQC_dict[yr]
                                                    
    # Do nocturnal optimisation for each window
    for date in step_date_array:

        # Get Eo for the relevant year and write to the parameters dictionary
        param_index = np.where(date_array == date)
        params_dict['Eo_default'] = param_dk_rslt_array[param_index, 0].item()

        # Subset the data and check length
        sub_dict = subset_window(data_dict, step_date_index_dict[date])
        noct_dict = subset_daynight(sub_dict, noct_flag = True)
        len_all_noct = len(noct_dict['NEE'])
        noct_dict = subset_nan(noct_dict)
        len_valid_noct = len(noct_dict['NEE'])
    
        # Do optimisation only if data passes minimum threshold
        if round(float(len_valid_noct) / len_all_noct * 100) > \
        configs_dict['min_pct_noct_window']:
     
            dark_rb_param, dark_rb_error_state = dark.optimise_rb(noct_dict, 
                                                                  params_dict)
        else:
            
            dark_rb_param, dark_rb_error_state = [np.nan], 10                                                      

        # Send data to the results arrays                                                      
        param_dk_rslt_array[param_index, 1] = dark_rb_param
        param_QC_array[param_index, 1] = dark_rb_error_state
    
    # Interpolate missing values
    param_dk_rslt_array[:, 1] = interp_params(param_dk_rslt_array[:, 1])

    # Do daytime optimisation with nocturnally fixed rb for each window    
    for date in step_date_array:

        # Get Eo and rb for the relevant year and write to the parameters dictionary
        param_index = np.where(date_array == date)
        params_dict['Eo_default'] = param_dk_rslt_array[param_index, 0].item()
        params_dict['rb_default'] = param_dk_rslt_array[param_index, 1].item()
        
        # Subset the data and check length
        sub_dict = subset_window(data_dict, step_date_index_dict[date])
        day_dict = subset_daynight(sub_dict, noct_flag = False)
        len_all_day = len(day_dict['NEE'])
        day_dict = subset_nan(day_dict)
        len_valid_day = len(day_dict['NEE'])        

        # Do optimisation only if data passes minimum threshold
        if round(float(len_valid_day) / len_all_day * 100) > \
        configs_dict['min_pct_noct_window']:
                
            light_params, light_error_state = light.optimise_no_rb \
                                              (day_dict, params_dict)

        else:
            
            light_params, light_error_state = [np.nan, np.nan, np.nan], 10
        
        # Send data to the results array                         
        param_dk_rslt_array[param_index, 2:5] = light_params
        param_QC_array[param_index, 2] = light_error_state
        
        # Write alpha default values to default dictionary if valid
        if not np.isnan(light_params[0]): params_dict['alpha_default'] = light_params[0]
    
    # Interpolate missing values
    param_dk_rslt_array[:, 2:5] = interp_params(param_dk_rslt_array[:, 2:5])

    # Do daytime optimisation with free rb for each window    
    for date in step_date_array:

        # Get Eo and rb for the relevant year and write to the parameters dictionary
        param_index = np.where(date_array == date)
        params_dict['Eo_default'] = param_lt_rslt_array[param_index, 0].item()
                
        # Subset the data and check length
        sub_dict = subset_window(data_dict, step_date_index_dict[date])
        day_dict = subset_daynight(sub_dict, noct_flag = False)
        len_all_day = len(day_dict['NEE'])
        day_dict = subset_nan(day_dict)
        len_valid_day = len(day_dict['NEE'])        

        # Do optimisation only if data passes minimum threshold
        if round(float(len_valid_day) / len_all_day * 100) > \
        configs_dict['min_pct_noct_window']:
                
            light_params, light_error_state = light.optimise_with_rb \
                                              (day_dict, params_dict)
         
        else:
            
            light_params, light_error_state = [np.nan, np.nan, np.nan, np.nan], 10
        
        # Send data to the results array                   
        param_lt_rslt_array[param_index, 1:5] = light_params
        param_QC_array[param_index, 3] = light_error_state
        
        # Write alpha default values to default dictionary if valid
        if not np.isnan(light_params[1]): params_dict['alpha_default'] = light_params[1]


#        # Print error messages if any
#        if dark_rb_error_state != 0 or light_error_state != 0:
#            print 'For ' + dt.datetime.strftime(date, '%Y-%m-%d') + ':'
#            if dark_rb_error_state != 0:            
#                print msg_dark_dict[dark_rb_error_state]
#            if light_error_state != 0:
#                print msg_light_dict[light_error_state]

#        # Write current parameters to estimated parameters dictionary
#        current_params = param_rslt_array[np.where(date_array == date), :6].reshape(6)
#        est_params_d = {param_rslt_list[i]: current_params[i] for i in range(6)}
#        
#        # Estimate the time series data and plot
#        est_series_d = estimate_Re(sub_dict, est_params_d)
#        if est_series_d != None and configs_dict['output_plots']:
#            combine_d = dict(sub_dict, **est_series_d)
#            plot_windows(combine_d, configs_dict, date)

    # Interpolate missing values
    param_lt_rslt_array[:, 1:5] = interp_params(param_lt_rslt_array[:, 1:5])

    # Loop through interpolated data and construct the time series
    for ind, date in enumerate(date_array):
        
        # Subset the data (single day)
        sub_dict = subset_window(data_dict, all_date_index_dict[date])

        # Get the indices for inserting the calculated data into the array
        start_ind = all_date_index_dict[date][0]
        end_ind = all_date_index_dict[date][1]

        # For each parameter set, write current parameters to estimated 
        # parameters dictionary, then estimate the time series
        for i, param_array in enumerate([param_dk_rslt_array, param_lt_rslt_array]):
            current_params = param_array[ind, :5]
            params_dict = {param_rslt_list[i]: current_params[i] for i in range(5)}
            series_dict = estimate_Re(sub_dict, params_dict)
            series_rslt_array[start_ind: end_ind + 1, i * 2] = series_dict['Re']
            series_rslt_array[start_ind: end_ind + 1, i * 2 + 1] = series_dict['GPP']
       
    # Make dictionaries
    params_dk_dict = {var + '_dk': param_dk_rslt_array[:, i] for i, var in enumerate(param_rslt_list)}
    params_lt_dict = {var + '_lt': param_lt_rslt_array[:, i] for i, var in enumerate(param_rslt_list)}
    params_dict = dict(params_dk_dict, **params_lt_dict)
    params_dict['date_time'] = date_array
    series_dict = {var: series_rslt_array[:, i] for i, var in enumerate(series_rslt_list)}
    series_dict['date_time'] = datetime_array    
    
    return params_dict, series_dict

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------        
# Estimate Re (from nocturnal and daytime parameters) and GPP
def estimate_Re(sub_d, params_d):

    # Estimate GPP and Re using nocturnal rb
    GPP, Re = light.LRF(sub_d, params_d['Eo'], params_d['rb'],
                        params_d['alpha'], params_d['Aopt'], params_d['k'])

    return {'Re': Re, 'GPP': GPP}

#------------------------------------------------------------------------------        
# Get dates
def get_dates(datetime_array, configs_dict):
    pdb.set_trace()
    # Assign configs to local vars
    window = configs_dict['window_size_days']

    # Create a series of continuous whole day dates that will be used for output
    # (parameter series will be interpolated between window centres)
    start_date = datetime_array[0].date()
    end_date = datetime_array[-1].date()
    num_days = (end_date - start_date).days + 1 # Add 1 so is inclusive of both end members
    all_dates_array = np.array([start_date + dt.timedelta(i) for i in xrange(num_days)])

    # Create a shifted array
    shift_mins = 60 * configs_dict['measurement_interval']
    shift_datetime_array = datetime_array - dt.timedelta(minutes = shift_mins)

    # Check that first and last days are complete and revise start and end dates if required
    temp_date = dt.datetime.combine((shift_datetime_array[0] + dt.timedelta(1)).date(), 
                                    dt.datetime.min.time())
    num_obs = len(np.where(shift_datetime_array < temp_date)[0])
    if num_obs < 24 * (1 / configs_dict['measurement_interval']):
        start_date = start_date + dt.timedelta(1)
    temp_date = dt.datetime.combine(shift_datetime_array[-1].date(), 
                                    dt.datetime.min.time())
    num_obs = len(np.where(shift_datetime_array >= temp_date)[0])
    if num_obs < 24 * (1 / configs_dict['measurement_interval']):
        end_date = end_date - dt.timedelta(1)
    
    # Calculate the dates that represent the centre of the window for each step
    num_days = (end_date - start_date).days + 1 - window # Add 1 so is inclusive of both end members
    first_fit_day = start_date + dt.timedelta(window / 2)
    step_days = np.arange(0, num_days, configs_dict['step_size_days'])
    step_dates_array = [first_fit_day + dt.timedelta(i) for i in step_days]

    # Make an index dictionary for step dates    
    step_dates_index_dict = {}
    for date in step_dates_array:    
        date_time = (dt.datetime.combine(date, dt.datetime.min.time()) 
                     + dt.timedelta(hours = 12))
        start_date = date_time - dt.timedelta(window / 2.0)
        end_date = date_time + dt.timedelta(window / 2.0)
        start_ind = np.where(datetime_array == start_date)[0].item() + 1
        end_ind = np.where(datetime_array == end_date)[0].item()
        step_dates_index_dict[date] = [start_ind, end_ind]
    
    # Make an index dictionary for all dates
    all_dates_index_dict = {}
    for date in all_dates_array:
        date_time = dt.datetime.combine(date, dt.datetime.min.time())
        if date == all_dates_array[0]:
            start_ind = 0
        else:
            start_date = date_time + dt.timedelta(hours = configs_dict['measurement_interval'])
            start_ind = np.where(datetime_array == start_date)[0].item()
        if date == all_dates_array[-1]:
            end_ind = len(datetime_array)
        else:
            end_date = date_time + dt.timedelta(1)
            end_ind = np.where(datetime_array == end_date)[0].item()
        all_dates_index_dict[date] = [start_ind, end_ind]
    
    # Make an index dictionary for years
    years_index_dict = {}
    year_array = np.array([i.year for i in shift_datetime_array])
    year_list = list(set(year_array))
    for yr in year_list:
        index = np.where(year_array == yr)[0]
        years_index_dict[yr] = [index[0], index[-1]]
    
    return step_dates_index_dict, all_dates_index_dict, years_index_dict

#------------------------------------------------------------------------------
# Interpolate parameters to create complete series
def interp_params(param_rslt_array):

    def do_interp(array_1D):
        xp = np.arange(len(arr))
        fp = array_1D[:]
        nan_index = np.isnan(fp)
        fp[nan_index] = np.interp(xp[nan_index], xp[~nan_index], fp[~nan_index])
        return fp   
    
    arr = param_rslt_array.copy()    
    num_vars = np.shape(arr)
    if len(num_vars) == 1:
        arr = do_interp(arr)
    else:
        num_vars = num_vars[1]
        for i in range(num_vars):
            arr[:, i] = do_interp(arr[:, i])

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
             7:'Value of nocturnal rb out of range - rejecting',
             8:'Value of Eo out of range - set to mean of other years or ' \
               'nearest range limit (50-400)',
             10:'Data did not pass minimum percentage threshold - ' \
                'skipping optimisation'}
    
    return d

# Create a dictionary with initial guesses for parameters
def make_initial_guess_dict(data_d):

    d = {}
    d['Eo_prior'] = 100
    d['k_prior'] = 0
    d['alpha_prior'] = -0.01
    index = np.where(data_d['PAR'] < 10)[0]
    d['rb_prior'] = np.nanmean(data_d['NEE'][index])
    index = np.where(data_d['PAR'] > 10)[0]
    d['Aopt_prior'] = (np.nanpercentile(data_d['NEE'][index], 5) - 
                 np.nanpercentile(data_d['NEE'][index], 95))
    return d

#------------------------------------------------------------------------------
# Data optimisation procedures

# Annual Eo
def optimise_annual_Eo(data_dict, params_dict, configs_dict, year_index_dict):
    
    # Initialise local variables with configurations
    min_pct = configs_dict['min_pct_annual']
    msmt_int = configs_dict['measurement_interval']
    
    # Get Eo for each year and compile dictionary
    yearsEo_dict = {}
    Eo_pass_keys = []
    Eo_range_fail_keys = []
    Eo_nan_fail_keys = []
    year_list = year_index_dict.keys()
    print 'Eo optimised using whole year is as follows:'
    for yr in year_list:

        # Calculate number of recs for year
        days = 366 if calendar.isleap(yr) else 365
        recs = days * (24 / msmt_int) / 2

        # Subset data        
        sub_dict = subset_window(data_dict, year_index_dict[yr])
        sub_dict = subset_nan(sub_dict)
        noct_flag = True
        sub_dict = subset_daynight(sub_dict, noct_flag)
        
        # Calculate percent of potential annual data that the subset contains
        pct = round(float(len(sub_dict['NEE'])) / recs * 100)

        # Fit L&T parameters if minimum data criterion satisfied, otherwise nan
        if pct > min_pct:
            params, error_code = dark.optimise_all(sub_dict, params_dict)
        else:
            params, error_code = [np.nan, np.nan], 10                                         

        # Assign year to pass, range_fail or nan_fail list for subsequent QC and fill
        Eo = params[0]
        yearsEo_dict[yr] = Eo
        if np.isnan(Eo):
            Eo_nan_fail_keys.append(yr)
        elif ((Eo < 50) | (Eo > 400)):
            Eo_range_fail_keys.append(yr)
        else:
            Eo_pass_keys.append(yr)

        print '    - ' + str(yr) + ': ' + str(round(params[0], 1))
    
    # Do QC on Eo
    yearsQC_dict = {yr: 0 for yr in year_list}
    if len(Eo_pass_keys) != len(yearsEo_dict):
        if len(Eo_nan_fail_keys) == len(yearsEo_dict):
            print 'Could not find any values of Eo for any years! Exiting...'
            sys.exit()
        elif len(Eo_pass_keys) != 0:
            Eo_mean = np.array([yearsEo_dict[i] for i in Eo_pass_keys]).mean()
            all_fail_keys = Eo_range_fail_keys + Eo_nan_fail_keys
            for i in (all_fail_keys):
                yearsEo_dict[i] = Eo_mean
                yearsQC_dict[i] = 1
            all_fail_keys = [str(key) for key in all_fail_keys]
            if len(all_fail_keys) > 1:
                all_fail_str = ', '.join(all_fail_keys)
            else:
                all_fail_str = all_fail_keys[0]
            print 'Eo optimisation failed for the following years: ' + all_fail_str
            print 'Eo for these years estimated from the mean of all other years'
        else:
            for i in Eo_range_fail_keys:
                yearsEo_dict[i] == 50 if yearsEo_dict[i] < 50 else 400
                yearsQC_dict[i] = 2
            Eo_mean = [yearsEo_dict[i] for i in Eo_range_fail_keys].mean()
            for i in Eo_nan_fail_keys:
                yearsEo_dict[i] = Eo_mean
                yearsQC_dict[i] = 3
            print 'Warning! Eo estimates were out of range for all years'
            print 'Low estimates have been set to lower limit (50);'
            print 'High estimates have been set to upper limit (400);'
            print 'Parameter estimates are unlikely to be robust!'
    else:
        print 'Eo estimates passed QC for all years'
        
    return yearsEo_dict, yearsQC_dict

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

    # Create night / day subsetting index and subset data
    if noct_flag:
        daynight_index = np.where(data_dict['PAR'] < 10)[0]
    else:
        daynight_index = np.where(data_dict['PAR'] > 10)[0]
    temp_array = temp_array[daynight_index]
    
    sub_dict = {var: temp_array[:, i] for i, var in enumerate(data_dict.keys())}

    return sub_dict

def subset_nan(data_dict):
    
    # Turn dictionary into an array
    temp_array = np.empty([len(data_dict['NEE']), len(data_dict)])
    for i, var in enumerate(data_dict.keys()):
        temp_array[:, i] = data_dict[var]

    # Create nan subsetting index and subset data and count
    QCdata_index = np.where(np.all(~np.isnan(temp_array), axis=1))    
    temp_array = temp_array[QCdata_index]

    sub_dict = {var: temp_array[:, i] for i, var in enumerate(data_dict.keys())}

    return sub_dict

# Subsetting of date window
def subset_window(data_dict, index_list):

    # Subset the arrays on basis of index list
    sub_dict = {}
    for i in data_dict.keys():
        sub_dict[i] = data_dict[i][index_list[0]: index_list[1] + 1]
    
    return sub_dict

#------------------------------------------------------------------------------