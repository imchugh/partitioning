# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 09:43:24 2015

@author: imchugh
"""

import numpy as np
import os
import pdb
import datetime as dt
import netCDF4
import xlrd

import Partition_NEE_3 as pt

# Open configuration and build dictionaries of config file contents
def get_configs():
    
    configs_dict = {'nan_value': -9999,
                    'minimum_temperature_spread': 5,
                    'step_size_days': 5,
                    'window_size_days': 15,
                    'minimum_pct_annual': 30,
                    'minimum_pct_noct_window': 20,
                    'minimum_pct_day_window': 50,
                    'output_plots': False,
                    'msmt_interval_hrs': 0.5,
                    'QC_accept_code': 0,
                    'plot_output_path': '/home/imchugh/Documents'}
                    
    return configs_dict
    
def get_data(configs_dict):

    # Specify file location and name
    data_input_path = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all'
    data_input_file = 'Whroo_2011_to_2014_L6.nc'
    data_input_target = os.path.join(data_input_path, data_input_file)

    # Initialise dictionaries
    vars_dict = {'carbon_flux':'Fc',
                 'solar_radiation':'Fsd',
                 'temperature':'Ta',
                 'vapour_pressure_deficit':'VPD'}
    QC_dict = {'carbon_flux':'Fc_QCFlag',
               'solar_radiation':'Fsd_QCFlag',
               'temperature':'Ta_QCFlag',
               'vapour_pressure_deficit':'VPD_QCFlag'}
    newNames_dict = {'carbon_flux':'NEE',
                     'solar_radiation':'PAR',
                     'temperature':'TempC',
                     'vapour_pressure_deficit':'VPD'}

    # Read .nc file
    d={}
    vars_list = vars_dict.values() + QC_dict.values()
    nc_obj = netCDF4.Dataset(data_input_target)
    date_time = np.array([dt.datetime(*xlrd.xldate_as_tuple(elem,0)) 
                          for elem in nc_obj.variables['xlDateTime']])
    for i in vars_list:
        ndims=len(nc_obj.variables[i].shape)
        if ndims==3:
            d[i]=nc_obj.variables[i][:,0,0]
        elif ndims==1:    
            d[i]=nc_obj.variables[i][:]
        d[i] = np.where(d[i] == configs_dict['nan_value'], np.nan, d[i])
    nc_obj.close()

    # Estimate PAR from Fsd
    d[vars_dict['solar_radiation']] = d[vars_dict['solar_radiation']] * 0.46 * 4.6

    # Replace configured error values with NaNs and remove data with unacceptable QC codes, 
    # then drop QC flag variables and rename variable to new names
    for key in vars_dict.keys():
        d[vars_dict[key]] = np.where(d[QC_dict[key]] != configs_dict['QC_accept_code'],
                                     np.nan, d[vars_dict[key]])
        d.pop(QC_dict[key])
        try:
            d[newNames_dict[key]] = d.pop(vars_dict[key])
        except KeyError:
            pdb.set_trace()
    
    # Add date_time variable
    d['date_time'] = date_time
          
    return d
    
# Get configurations and data
configs_dict = get_configs()
data_dict = get_data(configs_dict)

# Return parameter and series dictionaries
param_dict, series_dict = pt.main(data_dict, configs_dict)



