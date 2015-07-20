# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 17:28:19 2015

@author: imchugh
"""

import numpy as np
import os
import pandas as pd
import pdb

# Open configuration and build dictionaries of config file contents
def get_configs():
    
    configs_dict = {'minimum_temperature_spread': 5,
                    'step_size_days': 5,
                    'window_size_days': 15,
                    'min_pct_annual': 20,
                    'min_pct_noct_window': 20,
                    'min_pct_day_window': 50,
                    'output_plots': False,
                    'msmt_interval_hrs': 0.5,
                    'QC_accept_code': 0,
                    'plot_output_path': '/home/imchugh/Documents'}
                    
    return configs_dict
    
def get_data(configs_dict):

    # Specify file location and name
    data_input_path = r'E:'
    data_input_file = 'Advanced_processed_data_Whroo_v12a.df'
    data_input_target = os.path.join(data_input_path, data_input_file)
    print data_input_target
    # Initialise dictionaries
    vars_dict = {'carbon_flux':'Fc',
                 'solar_radiation':'Fsd_Con',
                 'temperature':'Ta',
                 'vapour_pressure_deficit':'VPD_Con',
                 'ustar': 'ustar'}
    newNames_dict = {'carbon_flux':'NEE',
                     'solar_radiation':'PAR',
                     'temperature':'TempC',
                     'vapour_pressure_deficit':'VPD',
                     'ustar': 'ustar'}
                     
    # Read .nc file
    d = {}
    df = pd.read_pickle(data_input_target)
    for key in vars_dict.keys():
        d[newNames_dict[key]] = np.array(df[vars_dict[key]])
    
    # Remove low u*    
    index = np.where(d[vars_dict['ustar']] < 0.4)
    d[newNames_dict['carbon_flux']][index] = np.nan

    # Add date_time variable
    d['date_time'] = df.index
    
    return d

# Get configurations and data
configs_dict = get_configs()
data_dict = get_data(configs_dict)

# Return parameter and series dictionaries
params_dict, series_dict = pt.main(data_dict, configs_dict)