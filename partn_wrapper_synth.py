# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 09:43:24 2015

@author: imchugh
"""

import Tkinter, tkFileDialog
from configobj import ConfigObj
import numpy as np
import os
import pdb
import Partition_NEE_3 as pt
import datetime as dt

# Open configuration and build dictionaries of config file contents
def get_configs():
    
    # Prompt user for configuration file and get it
    root = Tkinter.Tk(); root.withdraw()
    cfName = tkFileDialog.askopenfilename(initialdir='')
    root.destroy()
    cf=ConfigObj(cfName)
    
    # Map input file variable names to hard-coded names used in this script
    configs_dict = {}
    extNames_dict = dict(cf['variables'])
    intNames_dict = {'carbon_flux': 'NEE', 
                     'solar_radiation': 'PAR', 
                     'temperature': 'TempC', 
                     'vapour_pressure_deficit': 'VPD',
                     'date_time': 'date_time'}
                  
    configs_dict['variables'] = {intNames_dict[i]: 
                                 extNames_dict[i] 
                                 for i in intNames_dict.keys()}
    
    # Prepare dictionary of user settings - drop strings or change to int / float
    configs_dict['options'] = dict(cf['options'])
    for key in configs_dict['options']:
        if configs_dict['options'][key].isdigit():
            configs_dict['options'][key] = int(configs_dict['options'][key])
        elif configs_dict['options'][key] == 'True':
            configs_dict['options'][key] = True
        elif configs_dict['options'][key] == 'False':
            configs_dict['options'][key] = False
        else:
            try:
                configs_dict['options'][key] = float(configs_dict['options'][key])
            except ValueError:
                continue
    
    # Set output path
    configs_dict['paths'] = {'data_input_target': 
                             os.path.join(cf['files']['data_input_path'], 
                                          cf['files']['data_input_file']), 
                             'plot_output_path': 
                             os.path.join(cf['files']['plot_output_path'], 
                                          'Plots')}

    return configs_dict

def get_data(configs_dict):

    data_arr = np.load(configs_dict['paths']['data_input_target'])
    data_d = {varName: data_arr[configs_dict['variables'][varName]] 
              for varName in configs_dict['variables'].keys()}

    data_d['date_time'] = np.array(map(lambda x: 
                                       dt.datetime.strptime(x,'%Y-%m-%d %H:%M:%S'), 
                                       data_d['date_time']))
                  
    return data_d
    
# Get configurations and data
configs_dict = get_configs()
data_dict = get_data(configs_dict)
configs_dict['paths'].pop('data_input_target')
configs_dict.pop('variables')

# Return parameter and series dictionaries
params_dict, series_dict = pt.main(data_dict, configs_dict)

