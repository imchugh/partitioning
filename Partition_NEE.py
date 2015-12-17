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
sys.path.append('../Analysis_tools')
import DataIO as io

reload(dark)
reload(light)

# Main code starts here
def main(data_dict, configs_dict):

    # If user wants individual window plots, check whether output directories
    # are present, and create if not
    if configs_dict['output_plots']:
        output_path = os.path.join(configs_dict['output_path'], 'Window plots')
        configs_dict['window_plot_output_path'] = output_path
        if not os.path.isdir(output_path): 
            os.makedirs(output_path)
    
    # Get arrays of all datetimes, all dates and stepped dates
    datetime_array = data_dict.pop('date_time')
    (step_date_index_dict, 
     all_date_index_dict,
     year_index_dict) = get_dates(datetime_array, configs_dict)
    date_array = np.array(all_date_index_dict.keys())
    date_array.sort()
    step_date_array = np.array(step_date_index_dict.keys())
    step_date_array.sort()

    # Create variable name lists for results output
    series_rslt_list = ['Nocturnally derived Re', 'GPP from nocturnal derived Re', 
                        'Daytime derived Re', 'GPP from daytime derived Re']

    new_param_list = ['Eo', 'rb_noct', 'rb_day', 'alpha_fixed_rb', 
                      'alpha_free_rb', 'beta_fixed_rb', 'beta_free_rb', 
                      'k_fixed_rb', 'k_free_rb', 'Eo error code', 
                      'Nocturnal rb error code', 
                      'Light response parameters + fixed rb error code', 
                      'Light response parameters + free rb error code']

    # Create dictionaries for results
    # First the parameter estimates and error codes...
    empty_array = np.empty([len(date_array)])
    empty_array[:] = np.nan
    opt_params_dict = {var: empty_array.copy() for var in new_param_list}
    opt_params_dict['date'] = date_array  
    # Then the time series estimation
    empty_array = np.empty([len(datetime_array)])
    empty_array[:] = np.nan
    series_est_dict = {var: empty_array.copy() for var in series_rslt_list}
    series_est_dict['date_time'] = datetime_array

    # Create a dictionary containing initial guesses for each parameter
    params_dict = make_initial_guess_dict(data_dict)

    # Get the annual estimates of Eo
    print 'Optimising fit for Eo for each year...'
    Eo_dict, EoQC_dict = optimise_annual_Eo(data_dict, 
                                            params_dict, 
                                            configs_dict,
                                            year_index_dict)

    print 'Done!'

    # Write to result arrays
    year_array = np.array([i.year for i in date_array])
    for yr in year_array:
        index = np.where(year_array == yr)
        opt_params_dict['Eo'][index] = Eo_dict[yr]
        opt_params_dict['Eo error code'][index] = EoQC_dict[yr]

    # Rewrite the parameters dictionary so that there will be one set of 
    # defaults for the free and one set of defaults for the fixed parameters
    params_dict = {'fixed_rb': make_initial_guess_dict(data_dict),
                   'free_rb': make_initial_guess_dict(data_dict)}

    # Do nocturnal optimisation for each window
    print 'Optimising fit for rb using nocturnal data...'
    for date in step_date_array:

        # Get Eo for the relevant year and write to the parameters dictionary
        param_index = np.where(date_array == date)
        params_dict['fixed_rb']['Eo_default'] = opt_params_dict['Eo']\
                                                [param_index]

        # Subset the data and check length
        sub_dict = subset_window(data_dict, step_date_index_dict[date])
        
        # Subset again to remove daytime and then nan
        noct_dict = subset_daynight(sub_dict, noct_flag = True)
        len_all_noct = len(noct_dict['NEE'])
        noct_dict = subset_nan(noct_dict)
        len_valid_noct = len(noct_dict['NEE'])
    
        noct_dict['Fc_series'] = noct_dict['NEE']        
            
        # Do optimisation only if data passes minimum threshold
        if round(float(len_valid_noct) / len_all_noct * 100) > \
        configs_dict['minimum_pct_noct_window']:
            params, error_state = dark.optimise_rb(noct_dict, 
                                                   params_dict['fixed_rb'])
        else:
            params, error_state = [np.nan], 10                                                      

        # Send data to the results dict
        opt_params_dict['rb_noct'][param_index] = params
        opt_params_dict['Nocturnal rb error code'][param_index] = error_state
        
        # Estimate time series and plot if requested
        if error_state == 0 and configs_dict['output_plots']:        
            this_params_dict = {'Eo': opt_params_dict['Eo'][param_index],
                                'rb': opt_params_dict['rb_noct'][param_index]}
            est_series_dict = estimate_Re_GPP(sub_dict, this_params_dict)
            combine_dict = dict(sub_dict, **est_series_dict)
            plot_windows(combine_dict, configs_dict, date, noct_flag = True)

    # Interpolate
    opt_params_dict['rb_noct'] = interp_params(opt_params_dict['rb_noct'])

    print 'Done!'

    print 'Optimising fit for light response parameters using fixed and free rb...'

    # Now do daytime
    for date in step_date_array:

        # Get Eo for the relevant year and write to the parameters dictionary
        param_index = np.where(date_array == date)

        # Subset the data and check length
        sub_dict = subset_window(data_dict, step_date_index_dict[date])

        # Subset the data and check length
        sub_dict = subset_window(data_dict, step_date_index_dict[date])
        day_dict = subset_daynight(sub_dict, noct_flag = False)
        len_all_day = len(day_dict['NEE'])
        day_dict = subset_nan(day_dict)
        len_valid_day = len(day_dict['NEE'])        

        #----------------------------------------------------------------------
        # With fixed rb

        # Get Eo and rb for the relevant date
        params_dict['fixed_rb']['Eo_default'] = opt_params_dict['Eo']\
                                                [param_index]
        params_dict['fixed_rb']['rb_default'] = opt_params_dict['rb_noct']\
                                                [param_index]

        # Do optimisation only if data passes minimum threshold
        if round(float(len_valid_day) / len_all_day * 100) > \
        configs_dict['minimum_pct_noct_window']:
            params, error_state = light.optimise_fixed_rb(day_dict, 
                                                          params_dict['fixed_rb'])
        else:
            params, error_state = [np.nan, np.nan, np.nan], 10

        # Send data to the results dict
        for i, var in enumerate(['alpha_fixed_rb', 'beta_fixed_rb', 
                                 'k_fixed_rb']):
            opt_params_dict[var][param_index] = params[i]
        (opt_params_dict['Light response parameters + fixed rb error code']
         [param_index]) = error_state

        # Write alpha default values to default dictionary if valid
        if error_state < 2:
            params_dict['fixed_rb']['alpha_default'] = params[0]
        else:
            params_dict['fixed_rb']['alpha_default'] = 0    

        # Estimate time series and plot if requested
        if error_state == 0 and configs_dict['output_plots']:        
            this_params_dict = {'Eo': opt_params_dict['Eo'][param_index],
                                'rb': opt_params_dict['rb_noct'][param_index],
                                'alpha': opt_params_dict['alpha_fixed_rb'][param_index],
                                'beta': opt_params_dict['beta_fixed_rb'][param_index],
                                'k': opt_params_dict['k_fixed_rb'][param_index]}
            est_series_dict = estimate_Re_GPP(sub_dict, this_params_dict, GPP = True)
            est_series_dict['NEE'] = est_series_dict['Re'] + est_series_dict['GPP']
            combine_dict = dict(sub_dict, **est_series_dict)
            plot_windows(combine_dict, configs_dict, date, noct_flag = False)

        #----------------------------------------------------------------------
        # With free rb

        # Get Eo for the relevant date
        params_dict['free_rb']['Eo_default'] = opt_params_dict['Eo']\
                                                [param_index]
           
        # Do optimisation only if data passes minimum threshold
        if round(float(len_valid_day) / len_all_day * 100) > \
        configs_dict['minimum_pct_noct_window']:
            try:
                params, error_state = light.optimise_free_rb(day_dict, 
                                                             params_dict['free_rb'])
            except:
                pdb.set_trace()
        else:
            params, error_state = [np.nan, np.nan, np.nan, np.nan], 10
        
        # Send data to the results dicts                 
        for i, var in enumerate(['rb_day', 'alpha_free_rb', 'beta_free_rb', 
                                 'k_free_rb']):
            opt_params_dict[var][param_index] = params[i]
        (opt_params_dict['Light response parameters + free rb error code']
         [param_index]) = error_state
        
        # Write alpha default values to default dictionary if valid
        if error_state < 2:
            params_dict['free_rb']['alpha_default'] = params[1]
        else:
            params_dict['free_rb']['alpha_default'] = 0

    print 'Done!'
    
    print 'Running QC...'
    
    print 'Done!'
    
    print 'Interpolating...'

    # Interpolate missing values
    for var in new_param_list[2: 9]:
        opt_params_dict[var] = interp_params(opt_params_dict[var])

    print 'Done!'

    print 'Plotting interpolated time series of parameter estimates...'

    plot_parameter_series(opt_params_dict, configs_dict)
    
    print 'Done!'

    print 'Calculating time series from interpolated model parameters and ' \
          'drivers...'

    # Loop through interpolated data and construct the time series
    for ind, date in enumerate(date_array):
        
        # Subset the data (single day)
        sub_dict = subset_window(data_dict, all_date_index_dict[date])

        # Get the indices for inserting the calculated data into the array
        start_ind = all_date_index_dict[date][0]
        end_ind = all_date_index_dict[date][1]

        # For each parameter set, write current parameters to estimated 
        # parameters dictionary, then estimate the time series
        swap_dict = {'Eo': 'Eo',
                     'rb': 'rb_noct', 
                     'alpha': 'alpha_fixed_rb',
                     'beta': 'beta_fixed_rb',
                     'k': 'k_fixed_rb'}
        this_dict = {var: opt_params_dict[swap_dict[var]][ind]
                     for var in swap_dict.keys()}
        temp_dict = estimate_Re_GPP(sub_dict, this_dict, GPP = True)
        series_est_dict['Nocturnally derived Re'][start_ind: end_ind + 1] = \
        temp_dict['Re']
        series_est_dict['GPP from nocturnal derived Re'][start_ind: end_ind + 1] = \
        temp_dict['GPP']
        swap_dict = {'Eo': 'Eo',
                     'rb': 'rb_day', 
                     'alpha': 'alpha_free_rb',
                     'beta': 'beta_free_rb',
                     'k': 'k_free_rb'}
        this_dict = {var: opt_params_dict[swap_dict[var]][ind]
                     for var in swap_dict.keys()}
        temp_dict = estimate_Re_GPP(sub_dict, this_dict, GPP = True)
        series_est_dict['Daytime derived Re'][start_ind: end_ind + 1] = \
        temp_dict['Re']
        series_est_dict['GPP from daytime derived Re'][start_ind: end_ind + 1] = \
        temp_dict['GPP']        

    print 'Done'
    
    'Writing results to file...'
    
    # Rewrite the variable names for each of the parameter arrays to 
    # differentiate night and day, then join and drop Eo (duplicated),
    # then join with error codes and output to csv with specified key order   
    keyorder = ['Eo', 'rb_noct', 'rb_day', 'alpha_fixed_rb', 'alpha_free_rb',
                'beta_fixed_rb', 'beta_free_rb', 'k_fixed_rb', 'k_free_rb',
                'Eo error code', 'Nocturnal rb error code', 
                'Light response parameters + fixed rb error code', 
                'Light response parameters + free rb error code']
    params_fileout_name = os.path.join(configs_dict['output_path'], 
                                       'fit_parameters.csv')    
    io.array_dict_to_csv(opt_params_dict, params_fileout_name, 
                         keyorder = keyorder)
    
    # Now output the estimated time series
    ts_fileout_name = os.path.join(configs_dict['output_path'], 
                                   'estimated_Re_and_GPP.csv')   
    io.array_dict_to_csv(series_est_dict, ts_fileout_name, 
                         keyorder = series_rslt_list)
    
    'Done!'
                     
    return opt_params_dict, series_est_dict

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------        
# Estimate Re (from nocturnal and daytime parameters) and GPP
def estimate_Re_GPP(sub_d, params_d, GPP = False):

    return_dict = {}
    if GPP:
        GPP, Re = light.LRF_part(sub_d, params_d['Eo'], params_d['rb'],
                                 params_d['alpha'], params_d['beta'], 
                                 params_d['k'])
        return_dict['Re'] = Re
        return_dict['GPP'] = GPP
    else:
        Re = dark.TRF(sub_d, params_d['Eo'], params_d['rb'])
        return_dict['Re'] = Re
    return return_dict


#------------------------------------------------------------------------------        
# Get dates
def get_dates(datetime_array, configs_dict):

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
    
    d = {0:'Optimisation successful',
         1:'Value of k failed range check - setting to zero and ' \
           'recalculating other parameters',
         2:'Value of alpha failed range check - using previous ' \
           'estimate (zero if unavailable) and recalculating other ' \
           'parameters',
         3:'Optimisation reached maximum number of iterations' \
           'without convergence',
         4:'Value of beta and rb have wrong sign - ' \
           'rejecting all parameters',
         5:'Value of beta has wrong sign - rejecting all parameters',
         6:'Value of daytime rb has wrong sign - rejecting all parameters',
         7:'Value of nocturnal rb out of range - rejecting',
         8:'Value of Eo out of range - set to mean of other years or ' \
           'nearest range limit (50-400)',
         9:'Value of nocturnal rb has wrong sign - rejecting',
         10:'Data did not pass minimum percentage threshold - ' \
            'skipping optimisation'}
    
    return d

# Create a dictionary with initial guesses for parameters
def make_initial_guess_dict(data_d):

    # Calculate the parameter values that are intialised from data
    index = np.where(data_d['PAR'] < 10)[0]
    daytime_NEE_mean = np.nanmean(data_d['NEE'][index])
    daytime_NEE_range = (np.nanpercentile(data_d['NEE'][index], 5) - 
                         np.nanpercentile(data_d['NEE'][index], 95))

    params_dict = {'Eo_prior': 100,
                   'k_prior': 0,
                   'alpha_prior': -0.01,
                   'rb_prior': daytime_NEE_mean,
                   'beta_prior': daytime_NEE_range,
                   'alpha_default': 0,
                   'beta_default': 0,
                   'k_default': 0 }
    
    return params_dict

#------------------------------------------------------------------------------
# Data optimisation procedures

# Annual Eo
def optimise_annual_Eo(data_dict, params_dict, configs_dict, year_index_dict):
    
    # Initialise local variables with configurations
    min_pct = configs_dict['minimum_pct_annual']
    msmt_int = configs_dict['measurement_interval']
    
    # Get Eo for each year and compile dictionary
    yearsEo_dict = {}
    yearsQC_dict = {}
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
        sub_dict['Fc_series'] = sub_dict['NEE']
        if pct > min_pct:
            params, error_code = dark.optimise_all(sub_dict, params_dict)
        else:
            params, error_code = [np.nan, np.nan], 10                                         

        # Assign year to pass, range_fail or nan_fail list for subsequent QC and fill
        Eo = params[0]
        yearsEo_dict[yr] = Eo
        yearsQC_dict[yr] = error_code
        if np.isnan(Eo):
            Eo_nan_fail_keys.append(yr)
        elif ((Eo < 50) | (Eo > 400)):
            Eo_range_fail_keys.append(yr)
        else:
            Eo_pass_keys.append(yr)

        print '    - ' + str(yr) + ': ' + str(round(params[0], 1))
    
    # Do QC on Eo
    if len(Eo_pass_keys) != len(yearsEo_dict):
        if len(Eo_nan_fail_keys) == len(yearsEo_dict):
            print 'Could not find any values of Eo for any years! Exiting...'
            sys.exit()
        elif len(Eo_pass_keys) != 0:
            Eo_mean = np.array([yearsEo_dict[i] for i in Eo_pass_keys]).mean()
            all_fail_keys = Eo_range_fail_keys + Eo_nan_fail_keys
            for i in (all_fail_keys):
                yearsEo_dict[i] = Eo_mean
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
            Eo_mean = [yearsEo_dict[i] for i in Eo_range_fail_keys].mean()
            for i in Eo_nan_fail_keys:
                yearsEo_dict[i] = Eo_mean
            print 'Warning! Eo estimates were out of range for all years'
            print 'Low estimates have been set to lower limit (50);'
            print 'High estimates have been set to upper limit (400);'
            print 'Parameter estimates are unlikely to be robust!'
    else:
        print 'Eo estimates passed QC for all years'
        
    return yearsEo_dict, yearsQC_dict

#------------------------------------------------------------------------------
# Plots 

# Parameter estimates
def plot_parameter_series(params_dict, configs_dict):
    
    # Get years
    years = list(set([i.year for i in params_dict['date']]))
    year_dates = [dt.date(year, 1, 1) for year in years]
    year_dates = [date for date in year_dates if date in params_dict['date']]
    year_loc_dates = [dt.date(year, 7, 1) for year in years]
    year_loc_dates = [date for date in year_loc_dates if date 
                      in params_dict['date']]
    year_labels = [str(date.year) for date in year_loc_dates]

    # Get months
    markers_dict = {1: 'J', 2: 'F', 3: 'M', 4: 'A', 5: 'M', 6: 'J', 
                    7: 'J', 8: 'A', 9: 'S', 10: 'O', 11: 'N', 12: 'D'}
    month_dates = params_dict['date'][np.array([i.day == 1 for i in 
                                                params_dict['date']])]
    month_labels = [markers_dict[date.month] for date in month_dates]
    
    # Make dicts for each axes instance
    labels_dict = {0: ['rb_noct', 'rb_day', 'Nocturnal rb error code',
                       'Light response parameters + free rb error code',
                       r'$rb\/(\mu mol C\/m^{-2} s^{-1})$'],
                   1: ['alpha_fixed_rb','alpha_free_rb',
                       'Light response parameters + fixed rb error code',
                       'Light response parameters + free rb error code',
                       r'$\alpha\/(\mu mol C\/ photon^{-1} m^{-2} s^{-1})$'],
                   2: ['beta_fixed_rb','beta_free_rb',
                       'Light response parameters + fixed rb error code',
                       'Light response parameters + free rb error code',
                       r'$\beta\/(\mu mol C\/ m^{-2} s^{-1})$'],
                   3: ['k_fixed_rb','k_free_rb',
                       'Light response parameters + fixed rb error code',
                       'Light response parameters + free rb error code',
                       '$k\/(unitless)$']}

    # Make axes
    fig = plt.figure(figsize = (12, 32))
    fig.patch.set_facecolor('white')
    ax1 = fig.add_subplot(411)
#    plt.setp(ax1.get_xticklabels(), visible = False)
    ax2 = fig.add_subplot(412, sharex = ax1)
#    plt.setp(ax2.get_xticklabels(), visible = False)
    ax3 = fig.add_subplot(413, sharex = ax1)
#    plt.setp(ax3.get_xticklabels(), visible = False)
    ax4 = fig.add_subplot(414, sharex = ax1)
    ax4.set_xlabel('$Month$', fontsize = 18)
    ax4.xaxis.labelpad = 10
    
    # Do plotting (lines for interpolated series, markers for optimised estimates)
    for i, ax in enumerate([ax1, ax2, ax3, ax4]):
        if i == 0:
            ax1_twin = ax1.twiny()
            ax1_twin.set_xlim([params_dict['date'][0], params_dict['date'][-1]])
            ax1_twin.set_xticks(year_loc_dates)
            ax1_twin.set_xticklabels(year_labels, fontsize = 16)
            
        series1 = params_dict[labels_dict[i][0]]
        series1_error_code = params_dict[labels_dict[i][2]]
        series1_filt = series1.copy()        
        series1_filt[np.isnan(series1_error_code)] = np.nan
        series2 = params_dict[labels_dict[i][1]]
        series2_error_code = params_dict[labels_dict[i][3]]
        series2_filt = series2.copy()
        series2_filt[np.isnan(series2_error_code)] = np.nan        
        ax.plot(params_dict['date'], series1, color = 'black')
        ax.plot(params_dict['date'], 
                series1_filt, color = 'black', marker = 'o', 
                markerfacecolor = 'None')
        ax.plot(params_dict['date'], series2, color = 'black', 
                 linestyle = ':')
        ax.plot(params_dict['date'], 
                series2_filt, color = 'black', marker = 's', 
                markerfacecolor = 'None')
        ax.tick_params(axis = 'y', labelsize = 12)
        ax.set_ylabel(labels_dict[i][4], fontsize = 16)
        ax.set_xticks(month_dates)
        ax.set_xticklabels(month_labels, fontsize = 12)
        ax.xaxis.set_ticks_position('bottom')
        [ax.axvline(x = i, color = 'black') for i in year_dates]

    plot_out_name = 'optimised parameter time series.jpg'
    fig.savefig(os.path.join(configs_dict['output_path'], plot_out_name))
    plt.show()
    
    return
    

# Daytime and nocturnal fits for each window
def plot_windows(data_dict, configs_dict, date, noct_flag):

    # Set parameters from dicts
    path = configs_dict['window_plot_output_path']
    window = configs_dict['window_size_days']
    
    for i in range(2):
        sub_d = subset_daynight(data_dict, noct_flag)
        if noct_flag:
            daynight_ind = 'noct'
            x_lab = r'Temperature ($^{o}C$)'
            x_var = sub_d['TempC']
            y_var1 = sub_d['NEE']
            y_var2 = sub_d['Re']
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
        plt.ylabel(r'NEE ($\mu mol C\/m^{-2} s^{-1}$)', fontsize = 16)
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