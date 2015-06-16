# -*- coding: utf-8 -*-
"""
Created on Wed May 27 09:28:52 2015

@author: imchugh
"""

import numpy as np
from scipy.optimize import curve_fit
import datetime as dt
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

def make_LRF_2(Eo, k):
        def LRF(data_d,rb,alpha,Aopt):
            Aopt_VPD = Aopt * np.exp(-k * (data_d['VPD'] - 1))
            index = np.where(data_d['VPD'] <= 1)[0]
            Aopt_VPD[index] = Aopt
            GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + (alpha * data_d['PAR'] / Aopt_VPD))
            index = np.where(data_d['PAR'] < 10)[0]
            GPP[index] = 0
            Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['Ta'] + 46.02)))
            return GPP + Reco
        return LRF

def make_LRF_3(Eo, k, alpha):
        def LRF(data_d,rb,Aopt):
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
    
def optimise_dark(drivers_d, response, params_d, priors_d, algo_num):
    
    if not params_d == None:
        Eo = params_d['Eo']
    
    prior_Eo = priors_d['Eo']
    prior_rb = priors_d['rb']
    
    try:
        if algo_num == 0:
            params = curve_fit(TRF, drivers_d, response, p0 = [prior_Eo, prior_rb])[0]
        elif algo_num == 1:
            params = curve_fit(make_TRF(Eo), drivers_d, response, p0 = [prior_rb])[0]
    except RuntimeError:
        params = [np.nan, np.nan]
        
    return params

def optimise_light(drivers_d, response, params_d, priors_d, algo_num):
    
    prior_rb = priors_d['rb']
    prior_alpha = priors_d['alpha']
    prior_Aopt = priors_d['Aopt']
    prior_k = priors_d['k']
      
    try:
        if algo_num == 0:
            params = curve_fit(make_LRF(params_d['Eo']), 
                               drivers_d, response, 
                               p0 = [prior_rb, prior_alpha, prior_Aopt, prior_k])[0] 
        elif algo_num == 1:
            params = curve_fit(make_LRF_2(params_d['Eo'], params_d['k']), 
                               drivers_d, response, 
                               p0 = [prior_rb, prior_alpha, prior_Aopt])[0]
        elif algo_num == 2:
            params = curve_fit(make_LRF_3(params_d['Eo'], params_d['k'], params_d['alpha']), 
                               drivers_d, response, 
                               p0 = [prior_rb, prior_Aopt])[0]
    except RuntimeError:
        params = [np.nan, np.nan, np.nan, np.nan]

    return params

#------------------------------------------------------------------------------
# Subsetting
def subset_data(data_d, drivers, var_to_fit, noct_flag):
    
    if not type(drivers) == list:
        drivers = [drivers]
    temp_array = np.empty([len(data_d[var_to_fit]), len(drivers) + 1])
    for i, var in enumerate(drivers):
        temp_array[:, i] = data_d[var]
    temp_array[:, -1] = data_d[var_to_fit]

    if noct_flag:
        daynight_index = np.where(data_d['PAR'] < 10)[0]
    else:
        daynight_index = np.where(data_d['PAR'] > 10)[0]
        
    temp_array = temp_array[daynight_index]
    
    num_records = len(temp_array)
    
    QCdata_index = np.where(np.all(~np.isnan(temp_array), axis=1))    
    
    temp_array = temp_array[QCdata_index]
    
    valid_records = len(temp_array)
    
    percent_avail = round(float(valid_records) / num_records * 100, 1)
    
    drivers_d = {var: temp_array[:, i] for i, var in enumerate(drivers)}
    response = temp_array[:, -1]
    
    return drivers_d, response, percent_avail
    
#------------------------------------------------------------------------------
# Open configuration file and 
def get_configs():
    
    root = Tkinter.Tk(); root.withdraw()
    file_in = tkFileDialog.askopenfilename(initialdir='')
    root.destroy()
    return file_in

# Main code starts here


# Choose variables and years
var_to_fit = 'NEE_err'
window = 15
step = 5
flux_interval_hrs = 0.5

# Initialise parameters for data generation
Eo_init = 200
rb_init = 1.4
alpha_init = -0.1
Aopt_init = -16.0
k_init = 0.2

VPD_name = 'VPD'
T_name = 'Ta'
PAR_name = 'PAR'
NEE_name = 'NEE'

min_pct_annual = 30
min_pct_noct_window = 30
min_pct_day_window = 50

# Specify working directories and file names
working_dir = '/home/imchugh/Analysis/Whroo/Data/Flux_and_met/'
input_file = 'drivers.npz'

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

# Make the data array

# Create a dictionary with initial guesses for parameters
priors_d = {}
priors_d['Eo'] = 100
priors_d['k'] = 0
priors_d['alpha'] = -0.01
index = np.where(data_d['PAR'] < 10)[0]
priors_d['rb'] = data_d[var_to_fit][index].mean()
index = np.where(data_d['PAR'] > 10)[0]
priors_d['Aopt'] = (np.nanpercentile(data_d[var_to_fit][index], 5) - 
                    np.nanpercentile(data_d[var_to_fit][index], 95))

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
    drivers_d, response, pct = subset_data(temp_d, [T_name], var_to_fit, True)
    if pct > min_pct_annual:
        params = optimise_dark(drivers_d, response, None, priors_d, 0)
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
        print 'Warning! Eo estimates were out of range for all years'
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
rslt_arr = np.empty([len(all_whole_day_dates), 6])
rslt_arr[:] = np.nan
params_d = {'k': 0, 'alpha': 0}
initialise_alpha = True

for date in step_whole_day_dates:
    
    print date    
    
    if initialise_alpha:
        prev_alpha = 0
        initialise_alpha = False        
    
    # Find bracketing dates
    date_time = dt.datetime.combine(date, dt.datetime.min.time()) + dt.timedelta(hours = 12)
    start_date = date_time - dt.timedelta(window / 2.0)
    end_date = date_time + dt.timedelta(window / 2.0)
    
    # Get index for right dates, then subset the arrays
    index = np.where((date_arr >= start_date) & (date_arr < end_date))
    sub_d = {}
    for i in data_d.keys():
        sub_d[i] = data_d[i][index][1:]

    # Get Eo for the relevant year and write to the parameters dictionary
    params_d['Eo'] = yearsEo_d[str(date_time.year)] 
    Eo_year = yearsEo_d[str(date_time.year)] 

    # Do fitting to subsets - dark...
    drivers_d, response, pct = subset_data(sub_d, [T_name], var_to_fit, True)
    
    if pct > min_pct_noct_window:
        params = optimise_dark(drivers_d, response, params_d, priors_d, 1)
    else:
        params = np.nan
    rb_noct = params[0]
    
    # ... then light
    drivers_d, response, pct = subset_data(sub_d, [PAR_name, T_name, VPD_name], var_to_fit, False)
    
    if pct > min_pct_day_window:
        
        params = optimise_light(drivers_d, response, params_d, priors_d, 0)
        rb_day, alpha, Aopt, k = params[0], params[1], params[2], params[3]
               
        # If nan returned or a negative VPD coefficient, rerun with VPD set to zero
        if ((params[3] == np.nan) | (params[3] < 0)):
            print 'Value of k unacceptable - setting to zero'
            k = 0
            params = optimise_light(drivers_d, response, params_d, priors_d, 1)
            rb_day, alpha, Aopt = params[0], params[1], params[2]
    
        # If a positive or otherwise out of range value of alpha was returned,
        # rerun with previous value of alpha if available, otherwise set to zero
        if ((params[1] == np.nan) | (params[1] > 0) | (params[1] < -0.22)):
            print 'Value of alpha unacceptable - using previous estimate and recalculating other parameters'
            alpha = prev_alpha
            params_d['alpha'] = alpha
            params = optimise_light(drivers_d, response, params_d, priors_d, 2)
            rb_day, Aopt = params[0], params[1]
        
        # If Aopt or rb is of the wrong sign, reject all parameters
        if Aopt > 0 or rb_day < 0:
            print 'Value of Aopt or rb has wrong sign - rejecting all parameters'
            params = [np.nan, np.nan, np.nan, np.nan]
            rb_day, alpha, Aopt, k = params[0], params[1], params[2], params[3]
        
        params = [rb_day, alpha, Aopt, k]        
        
        # Increment the trailing alpha parameter
        prev_alpha = alpha
    
    else:
        
        params = [np.nan, np.nan, np.nan, np.nan]
        rb_day, alpha, Aopt, k = params[0], params[1], params[2], params[3]
        
    # Insert into results array
    params = np.append(np.array([Eo_year, rb_noct]), params)
    index = np.where(all_whole_day_dates == date)
    rslt_arr[index, :] = params