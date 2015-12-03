# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 11:23:04 2015

@author: imchugh
"""
import sys
import os
import numpy as np
import pdb
import calendar

sys.path.append('../Analysis_tools')
sys.path.append('../Partitioning')
import DataIO as io
import gap_filling as gf
import dark_T_response_functions as dark
import datetime_functions as dtf
import data_filtering as filt

# Local functions
#------------------------------------------------------------------------------
# This removes daytime data then removes all records containing nan
def filtering(this_dict):
    noct_dict = filt.subset_arraydict_on_threshold(this_dict, 'Fsd', 5, '<', 
                                                   drop = True)    
    ustar_dict = filt.subset_arraydict_on_threshold(noct_dict, 'ustar', 0.42, 
                                                     '>', drop = True)
    sub_dict = filt.subset_arraydict_on_nan(ustar_dict)
    return sub_dict

# This gets the data
def get_data(configs_dict):

    # Initialise name change dictionary with new names via common keys
    vars_dict = configs_dict['variables']
    newNames_dict = {'carbon_flux':'NEE',
                     'temperature': 'TempC',
                     'solar_radiation':'Fsd',
                     'vapour_pressure_deficit': 'VPD',
                     'friction_velocity': 'ustar'}

    # Get file extension and target
    paths_dict = configs_dict['files']
    ext = os.path.splitext(paths_dict['input_file'])[1]
    data_input_target = os.path.join(paths_dict['input_path'],
                                     paths_dict['input_file'])

    # get data (screen only the Fc data to obs only)
    if ext == '.nc':
        Fc_dict = io.OzFluxQCnc_to_data_structure(data_input_target,
                                                  var_list = [vars_dict['carbon_flux']],
                                                  QC_accept_codes = [0])
        date_time = Fc_dict.pop('date_time')
        ancillary_vars = [vars_dict[var] for var in vars_dict.keys() 
                          if not var == 'carbon_flux']
        ancillary_dict, global_attr = io.OzFluxQCnc_to_data_structure(
                                          data_input_target,
                                          var_list = ancillary_vars,
                                          return_global_attr = True)
        data_dict = dict(Fc_dict, **ancillary_dict)
    elif ext == '.df':
        data_dict, global_attr = io.DINGO_df_to_data_structure(
                                     data_input_target,
                                     var_list = vars_dict.values(),
                                     return_global_attr = True)
        date_time = data_dict.pop('date_time')
                                       
    # Reconstruct data dict with standard names used by algorithms
    new_dict = {}
    for key in vars_dict.keys():
        new_dict[newNames_dict[key]] = data_dict.pop(vars_dict[key])  
    new_dict['date_time'] = date_time

    return new_dict, global_attr
    
def estimate_Re_time_series(data_dict, params_dict):
    
    return
    
    
#------------------------------------------------------------------------------

# Get configurations
configs_dict = io.config_to_dict(io.file_select_dialog())

# Get data
data_dict, attr = get_data(configs_dict)

# Drop levels in config file
configs_dict['configs']['output_path'] = configs_dict['files']['output_path']
configs_dict = configs_dict['configs']
window = configs_dict['window_size_days']
step = configs_dict['step_size_days']

# Get data step indices
datetime_array = data_dict.pop('date_time')
years_index_dict = dtf.get_year_indices(datetime_array)
dates_index_dict = dtf.get_day_indices(datetime_array)
step_dates_index_dict = dtf.get_moving_window_indices(datetime_array, 
                                                      window, step)

# Create a dates array to use as index for results
date_array = np.array(dates_index_dict.keys())
date_array.sort()
rb_array = np.empty([len(date_array)])
rb_array[:] = np.nan
Eo_array = np.empty([len(date_array)])

# Initalise dicts
all_noct_dict = filtering(temp_dict)
params_dict = {'Eo_prior': 100,
               'rb_prior': all_noct_dict['NEE'].mean()}

# Do annual fits for Eo
Eo_annual_data_dict = {}
Eo_annual_error_dict = {}
Eo_pass_keys = []
Eo_fail_keys = []
for year in years_index_dict.keys():

    # Calculate number of nocturnal recs for year
    days = 366 if calendar.isleap(year) else 365
    recs = days * (24 / configs_dict['measurement_interval']) / 2

    # Input and output indices
    

    # Calculate Eo
    indices = years_index_dict[year]
    this_dict = {var: data_dict[var][indices[0]: indices[1]] 
                 for var in data_dict.keys()}
    sub_dict = filtering(this_dict)
    data_pct = int(len(sub_dict['NEE']) / float(recs) * 100)
    if not data_pct < configs_dict['minimum_pct_annual']:
        params, error_code = dark.optimise_all(sub_dict, params_dict)
    else:
        params, error_code = [np.nan, np.nan], 10
    Eo_annual_data_dict[year] = params[0]
    Eo_annual_error_dict[year] = error_code
    if error_code == 0: 
        Eo_pass_keys.append(year)
    else:
        Eo_fail_keys.append(year)
        
    # Send to Eo results array
    out_index

# Fill any gaps
if np.all(np.isnan(Eo_annual_data_dict.values())):
    print 'Could not find any values of Eo for any years! Exiting...'
    sys.exit()
if np.any(np.isnan(Eo_annual_data_dict.values())):
    Eo_mean = np.array([Eo_annual_data_dict[year] for year in Eo_pass_keys]).mean()    
    for year in Eo_fail_keys:
        Eo_annual_data_dict[year] = Eo_mean

# Calculate rb for windows    
for date in step_dates_index_dict.keys():
    
    # Specify Eo value for the relevant year
    params_dict['Eo_default'] = Eo_annual_data_dict[date.year]
    
    # Input and output indices
    in_indices = step_dates_index_dict[date]
    out_index = np.where(date_array == date)

    this_dict = {var: data_dict[var][in_indices[0]: in_indices[1]] 
                 for var in data_dict.keys()}
    sub_dict = filtering(this_dict)
    data_pct = int(len(sub_dict['NEE']) / float((24 / configs_dict['measurement_interval']) / 2) * 100)
    if not data_pct < configs_dict['minimum_pct_noct_window']:
        params, error_code = dark.optimise_rb(sub_dict, params_dict)
    else:
        params, error_code = [np.nan, np.nan], 10
    rb_array[out_index] = params

# Interpolate    
rb_array = gf.generic_2d_linear(rb_array)

# Estimate time series Re