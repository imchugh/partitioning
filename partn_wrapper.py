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

import Partition_NEE_6 as pt
import sys
sys.path.append('../Analysis_tools')
import DataIO as io


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

def main():    

    reload(pt)

    # Get configurations
    configs_dict = io.config_to_dict(io.file_select_dialog())
    
    # Get data
    data_dict, attr = get_data(configs_dict)
    configs_dict['configs']['output_path'] = configs_dict['files']['output_path']
    configs_dict = configs_dict['configs']
    
    # Estimate incoming PAR from Fsd
    data_dict['PAR'] = data_dict['Fsd'] * 0.46 * 4.6

    # Screen low ustar values (need to fix this to iterate on year!!!)
    data_dict['NEE'][(data_dict['Fsd'] < 5) & 
                    (data_dict['ustar'] < 0.4)] = np.nan

    # Return parameter and series dictionaries
    param_dict, series_dict = pt.main(data_dict, configs_dict)

    return param_dict, series_dict

