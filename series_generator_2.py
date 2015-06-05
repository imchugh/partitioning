# -*- coding: utf-8 -*-
"""
Created on Wed May 27 09:28:52 2015

@author: imchugh
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import pdb
import os

#------------------------------------------------------------------------------
# Data generation algorithm

def make_data(data_d,Eo,rb,alpha,Aopt,k):

    Aopt_VPD = Aopt * np.exp(-k * (data_d['VPD'] - 1))
    index = np.where(data_d['VPD'] <= 1)[0]
    Aopt_VPD[index] = Aopt
    GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + (alpha * data_d['PAR'] / Aopt_VPD))
    index = np.where(data_d['PAR'] < 10)[0]
    GPP[index] = 0
    Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['Ta'] + 46.02)))
    return GPP + Reco
    
#def random_error(data_d, stats_d):
#    
#    sigma_delta = np.where(data_d['Fc_syn'] > 0, 
#                           data_d['Fc_syn'] * stats_d['noct_slope'] + stats_d['noct_intcpt'],
#                           data_d['Fc_syn'] * stats_d['day_slope'] + stats_d['day_intcpt'])     
#
#    error = (np.random.laplace(0, sigma_delta) / 2).astype(np.float32)
#
#    return data_d['Fc_syn'] + error
    
#------------------------------------------------------------------------------
# Data optimisation algorithm
    
def TRF(data_d,Eo,rb):

    return rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['Ta'] + 46.02)))

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
    
def optimise_dark(data_d, var_to_fit):

    index = np.where(data_d['PAR'] < 10)[0]
    
    temp_d = {}
    for i in data_d.keys():
        temp_d[i] = data_d[i][index]

    drivers_d = {'Ta': temp_d['Ta']}
    response_d = temp_d[var_to_fit]

    params = curve_fit(TRF, drivers_d, response_d, p0 = [100, 1])

    return params[0]
    
def optimise_light(data_d, Eo, var_to_fit):
    
    index = np.where(data_d['PAR'] > 10)[0]    

    temp_d = {}
    for i in data_d.keys():
        temp_d[i] = data_d[i][index]

    drivers_d = {i: temp_d[i] for i in ['PAR','Ta','VPD']}
    response_d = temp_d[var_to_fit]

    params = curve_fit(make_LRF(Eo), drivers_d, response_d, p0 = [1, -1, -10, 1])
    
    return params[0]
    
#------------------------------------------------------------------------------
    
def pandasdf_to_dict(df):
    
    d = {}
    for i in df.columns:
        d[i] = np.array(df[i])
        
    return d

# Choose variables and years
var_to_fit = 'NEE'
year_to_pass = 2013
window = 10
step = 2

# Initialise parameters for data generation
Eo_init = 200
rb_init = 1.4
alpha_init = -0.1
Aopt_init = -16.0
k_init = 0.2

# Specify working directories and file names
working_dir = '/home/imchugh/Analysis/Whroo/Data/Flux_and_met/'
input_file = 'drivers.npz'

# Get the data and make a dict-based data structure
target = os.path.join(working_dir, input_file)
data_arr = np.load(target)
date_time = data_arr['date_time']
year = np.array([dt.datetime.strptime(i, '%Y-%m-%d %H:%M:%S').year for i in date_time])
index = np.where(year==year_to_pass)
data_d = {}
for item in data_arr.files:
    data_d[item] = data_arr[item][index]

# Generate NEE    
data_d['NEE'] = make_data(data_d, 200, 1.4, -0.1, -16, 0.2)

# Optimise Eo over whole dataset
params = optimise_dark(data_d, var_to_fit)
Eo_opt = params[0]
print 'Eo derived from all generated nocturnal NEE data is: ' + str(round(Eo_opt, 1))
print 'Original Eo parameter value is: ' + str(Eo_init)

# Optimise other parameters over whole dataset
params = optimise_light(data_d, Eo_opt, var_to_fit)
rb_opt, alpha_opt, Aopt_opt, k_opt = params[0], params[1], params[2], params[3] 
print 'rb, alpha, Aopt and k derived from generated daytime NEE data are: '
print (str(round(rb_opt, 1)) + ', ' + str(round(alpha_opt, 1)) + ', ' +
       str(round(Aopt_opt, 1)) + ', ' + str(round(k_opt, 1)) + ', ')
print 'Original rb, alpha, Aopt and k parameter values are: '
print str(rb_init) + ', ' + str(alpha_init) + ', ' + str(Aopt_init) + ', ' + str(k_init)

# Calculate dates for windows
num_days_year = 366 if year_to_pass % 4 == 0 else 365
date_buffer = window / 2 if window % 2 == 0 else (window - 1) / 2
time_del_days = 0.5 if window % 2 == 0 else 0
step_days = np.arange(date_buffer + 1, num_days_year - date_buffer + 1, step)
step_dates = [dt.datetime(year_to_pass, 1, 1) + dt.timedelta(i - 1 + time_del_days) for i in step_days]
 
def main():
        
    df=pd.read_pickle('/home/imchugh/Analysis/Whroo/Data/Flux_and_met/synthetic_partitioning.df')
    
    
    rslt_lst = []    
    
    for yr in ['2012','2013','2014']:
        for mn in range(12):
            indexer = yr + '-' + str(mn + 1)
            data_d = pandasdf_to_dict(df.loc[indexer])
            try:
                rslt_lst.append(optimise_light(data_d, 'Fc_syn_error'))
            except RuntimeError:
                rslt_lst.append([np.nan, np.nan, np.nan, np.nan])
    
    arr = np.vstack(rslt_lst)
        
    return arr
            


