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
            params = [np.nan, np.nan]
    
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
flux_interval_hrs = 0.5

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
data_d = {item: data_arr[item] for item in data_arr.files}

# Get the random error data and make a dict-based structure
target = os.path.join(working_dir, random_error_file)
data_arr = np.load(target)
stats_d = {}
for var in data_arr.files:
    stats_d[var] = data_arr[var].item(0)

# Generate NEE    
data_d['NEE'] = make_data(data_d, Eo_init, rb_init, alpha_init, Aopt_init, k_init)
data_d['NEE_err'] = data_d['NEE'] + random_error(data_d, stats_d)

# Generate a date series
# Note that we subtract one measurement interval in minutes from the date series because
# the last valid case in the 24 hours of data occurs at midnight (assuming that the
# naming convention is for the timestamp signifying the end of the averaging interval).
# So, for example, assuming half hourly interval, the first valid case is 0030 
# (average of 0000-0030) and final valid case is 0000 (average of 2330 to 0000).
# So any whole day date-based indexing will not capture the true day unless the 
# timestamp is lagged by one measurement interval.
date_arr = np.array([dt.datetime.strptime(i, '%Y-%m-%d %H:%M:%S') for i in data_d['date_time']])
date_arr = date_arr - dt.timedelta(minutes = 60 * flux_interval_hrs)

# Get Eo for each year
yearsEo_d = {}
optimise_Eo = True
year_arr = np.array([i.year for i in date_arr])
year_list = list(set(year_arr))
print 'Eo optimised using whole year is as follows:'
for yr in year_list:
    index = np.where(year_arr == yr)
    temp_d = {}
    for item in data_d.keys():
        temp_d[item] = data_d[item][index]
    params = optimise_dark(temp_d, Eo_init, var_to_fit, optimise_Eo)
    yearsEo_d[str(yr)] = params[0]
    print '    - ' + str(yr) + ': ' + str(round(params[0]))
optimise_Eo = False

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
        print 'Warning! Eo estimates were out of range for the following years: '
        print [i for i in (Eo_range_fail_keys)]
        print 'Low estimates have been set to lower limit (50);'
        print 'High estimates have been set to upper limit (400);'
        print 'Parameter estimates are unlikely to be robust!'
else:
    print 'Eo estimates passed QC for all years'

# Create a series of continuous whole day dates that will be used for output
# (parameter series will be interpolated between window centres)
start_date = date_arr[0].date()
end_date = date_arr[-1].date()
num_days = (end_date - start_date).days + 1 # Add 1 so is inclusive of both end members
all_whole_day_dates = np.array([start_date + dt.timedelta(i) for i in xrange(num_days)])

# Check that first and last days are complete and revise start and end dates if required
all_dates = np.array([i.date() for i in date_arr])
recs_required = 24 * (1 / flux_interval_hrs)
recs_count = 0
loop_count = 0
while recs_count != recs_required:
    recs_count = len(np.where(all_dates == all_whole_day_dates[loop_count])[0])
    loop_count = loop_count + 1
start_date = all_whole_day_dates[loop_count - 1]
recs_count = 0
loop_count = -1
while recs_count != recs_required:
    recs_count = len(np.where(all_dates == all_whole_day_dates[loop_count])[0])
    loop_count = loop_count - 1
end_date = all_whole_day_dates[loop_count + 1]

# Calculate the dates that represent the centre of the window for each step
num_days = (end_date - start_date).days + 1 - window # Add 1 so is inclusive of both end members
first_fit_day = start_date + dt.timedelta(window / 2)
step_days = np.arange(0, num_days, step)
step_whole_day_dates = [first_fit_day + dt.timedelta(i) for i in step_days]

# Initialise result arrays and step through data windows, find parameters for each
rslt_list = []
rslt_arr = np.empty([len(all_whole_day_dates), 6])
rslt_arr[:] = np.nan
for date in step_whole_day_dates:
    
    # Find bracketing dates
    date_time = dt.datetime.combine(date, dt.datetime.min.time()) + dt.timedelta(hours = 12)
    start_date = date_time - dt.timedelta(window / 2.0)
    end_date = date_time + dt.timedelta(window / 2.0)
    
    # Get index for right dates, then subset the arrays
    index = np.where((date_arr >= start_date) & (date_arr < end_date))
    sub_d = {}
    for i in data_d.keys():
        sub_d[i] = data_d[i][index][1:]

    # Get Eo for the relevant year
    Eo_year = yearsEo_d[str(date_time.year)] 

    # Do fitting to subsets - dark...
    params = optimise_dark(sub_d, Eo_year, var_to_fit, optimise_Eo)
    rb_noct = params[0]

    # ... then light
    params = optimise_light(sub_d, Eo_year, var_to_fit)

    # Do QC
    

    # Insert into results array
    params = np.append(np.array([Eo_year, rb_noct]), params)
    index = np.where(all_whole_day_dates == date)
    rslt_arr[index, :] = params
    
#rslt_arr = np.vstack(rslt_list)
