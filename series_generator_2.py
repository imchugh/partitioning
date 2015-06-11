# -*- coding: utf-8 -*-
"""
Created on Wed May 27 09:28:52 2015

@author: imchugh
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import datetime as dt
import pdb
import os
import sys

#------------------------------------------------------------------------------
# Data generation algorithms

def make_data(data_d,Eo,rb,alpha,Aopt,k):

    Aopt_VPD = Aopt * np.exp(-k * (data_d['VPD'] - 1))
    index = np.where(data_d['VPD'] <= 1)[0]
    Aopt_VPD[index] = Aopt
    GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + (alpha * data_d['PAR'] / Aopt_VPD))
    index = np.where(data_d['PAR'] < 10)[0]
    GPP[index] = 0
    Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['Ta'] + 46.02)))
    return GPP + Reco
    
def random_error(data_d, stats_d):
    
    sigma_delta = np.where(data_d['NEE'] > 0, 
                           data_d['NEE'] * stats_d['noct_slope'] + stats_d['noct_intcpt'],
                           data_d['NEE'] * stats_d['day_slope'] + stats_d['day_intcpt'])     

    return (np.random.laplace(0, sigma_delta) / 2).astype(np.float32)
    
#------------------------------------------------------------------------------
# Data optimisation algorithms
    
def TRF(data_d,Eo,rb):
    return rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['Ta'] + 46.02)))

def make_TRF(Eo):
    def TRF(data_d,rb):
        return rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['Ta'] + 46.02)))
    return TRF

def make_LRF(Eo):
        def LRF(data_d,rb,alpha,Aopt,k):
            Aopt_VPD = Aopt * np.exp(-k * (data_d['VPD'] - 1))
            index = np.where(data_d['VPD'] <= 1)[0]
            Aopt_VPD[index] = Aopt
            GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + (alpha * data_d['PAR'] / Aopt_VPD))
            index = np.where(data_d['PAR'] < 10)[0]
            GPP[index] = 0
            Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['Ta'] + 46.02)))
            return GPP + Reco
        return LRF
    
#------------------------------------------------------------------------------
# Data optimisation procedure 
    
def optimise_dark(data_d, Eo, var_to_fit, Eo_flag):

    index = np.where(data_d['PAR'] < 10)[0]
    
    temp_d = {}
    for i in data_d.keys():
        temp_d[i] = data_d[i][index]

    drivers_d = {'Ta': temp_d['Ta']}
    response_d = temp_d[var_to_fit]

    try:
        if Eo_flag:
            params = curve_fit(TRF, drivers_d, response_d, p0 = [100, 1])[0]
        else:
            params = curve_fit(make_TRF(Eo), drivers_d, response_d, p0 = [1])[0]
    except RuntimeError:
        if Eo_flag:
            params = [np.nan, np.nan]
        else:
            params = [np.nan]
    
    return params
    
def optimise_light(data_d, Eo, var_to_fit):
    
    index = np.where(data_d['PAR'] > 10)[0]    

    temp_d = {}
    for i in data_d.keys():
        temp_d[i] = data_d[i][index]

    drivers_d = {i: temp_d[i] for i in ['PAR','Ta','VPD']}
    response_d = temp_d[var_to_fit]

    try:
        params = curve_fit(make_LRF(Eo), drivers_d, response_d, p0 = [1, -1, -10, 1])[0]
    except RuntimeError:
        params = [np.nan, np.nan, np.nan, np.nan]
    
    return params
    
#------------------------------------------------------------------------------

# Choose variables and years
var_to_fit = 'NEE'
year_to_pass = 2013
window = 15
step = 5

# Initialise parameters for data generation
Eo_init = 200
rb_init = 1.4
alpha_init = -0.1
Aopt_init = -16.0
k_init = 0.2

# Specify working directories and file names
working_dir = '/home/imchugh/Analysis/Whroo/Data/Flux_and_met/'
input_file = 'drivers.npz'
random_error_file = 'random_error_stats.npz'

# Get the data and make a dict-based data structure
target = os.path.join(working_dir, input_file)
data_arr = np.load(target)
date_time = data_arr['date_time']
year = np.array([dt.datetime.strptime(i, '%Y-%m-%d %H:%M:%S').year for i in date_time])
index = np.where(year==year_to_pass)
data_d = {}
for item in data_arr.files:
    data_d[item] = data_arr[item][index]

# Get the random error data and make a dict-based structure
target = os.path.join(working_dir, random_error_file)
data_arr = np.load(target)
stats_d = {}
for var in data_arr.files:
    stats_d[var] = data_arr[var].item(0)

# Generate NEE    
data_d['NEE'] = make_data(data_d, Eo_init, rb_init, alpha_init, Aopt_init, k_init)
data_d['NEE_err'] = data_d['NEE'] + random_error(data_d, stats_d)

# Optimise Eo over whole dataset
optimise_Eo = True
params = optimise_dark(data_d, Eo_init, var_to_fit, optimise_Eo)
optimise_Eo = False
Eo_year = params[0]
print 'Eo derived from all generated nocturnal NEE data is: ' + str(round(Eo_year, 1))
print 'Original Eo parameter value is: ' + str(Eo_init)

## Optimise other parameters over whole dataset
#params = optimise_light(data_d, Eo_opt, var_to_fit)
#rb_opt, alpha_opt, Aopt_opt, k_opt = params[0], params[1], params[2], params[3] 
#print 'rb, alpha, Aopt and k derived from generated daytime NEE data are: '
#print (str(round(rb_opt, 1)) + ', ' + str(round(alpha_opt, 1)) + ', ' +
#       str(round(Aopt_opt, 1)) + ', ' + str(round(k_opt, 1)) + ', ')
#print 'Original rb, alpha, Aopt and k parameter values are: '
#print str(rb_init) + ', ' + str(alpha_init) + ', ' + str(Aopt_init) + ', ' + str(k_init)

# Calculate dates for windows and for output
num_days_year = 366 if year_to_pass % 4 == 0 else 365
date_buffer = (window + 1) / 2.0 if window % 2 == 0 else (window - 1) / 2.0
time_del_days = date_buffer % int(date_buffer)
date_buffer = int(date_buffer)
first_day = 1 + date_buffer
last_day = num_days_year - date_buffer
step_days = np.arange(first_day, last_day + 1, step) # Add one here because arange includes only half open interval
step_dates = [dt.datetime(year_to_pass, 1, 1) + 
              dt.timedelta(i - 1 + time_del_days) for i in step_days]
rslt_dates = np.array([dt.datetime(year_to_pass, 1, 1) + 
                       dt.timedelta(i) for i in xrange(num_days_year)])

# Initialise result arrays and step through data windows, find parameters for each
rslt_list = []
rslt_arr = np.empty([num_days_year, 5])
rslt_arr[:] = np.nan
data_dates = np.array([dt.datetime.strptime(i, '%Y-%m-%d %H:%M:%S') for i in data_d['date_time']])
for date_ in step_dates:
    
    # Find bracketing dates
    start_date = date_ - dt.timedelta(date_buffer)
    end_date = date_ + dt.timedelta(date_buffer)
    
    # Get index for right dates, then subset the arrays
    index = np.where((data_dates >= start_date) & (data_dates <= end_date))
    sub_d = {}
    for i in data_d.keys():
        sub_d[i] = data_d[i][index][1:]

    # Do fitting to subsets - dark...
    params = optimise_dark(sub_d, Eo_year, var_to_fit, optimise_Eo)
    rb_noct = params[0]

    # ... then light
    params = optimise_light(sub_d, Eo_year, var_to_fit)

    # Insert into results array
    params = np.append(rb_noct, params)
    index = np.where(rslt_dates == date_)
    rslt_arr[index, :] = params
    
    rslt_list.append(np.append(Eo_opt, params))
    
#rslt_arr = np.vstack(rslt_list)
