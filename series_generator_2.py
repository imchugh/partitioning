# -*- coding: utf-8 -*-
"""
Created on Wed May 27 09:28:52 2015

@author: imchugh
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import pdb

#------------------------------------------------------------------------------

def make_data(data_d,Eo,rb,alpha,Aopt,k):

    Aopt_VPD = Aopt * np.exp(-k * (data_d['VPD'] - 1))
    index = np.where(data_d['VPD'] <= 1)[0]
    Aopt_VPD[index] = Aopt
    GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + (alpha * data_d['PAR'] / Aopt_VPD))
    index = np.where(data_d['PAR'] < 20)[0]
    GPP[index] = 0
    Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['Ta'] + 46.02)))
    return GPP + Reco
    
def random_error(data_d, stats_d):
    
    sigma_delta = np.where(data_d['Fc_syn'] > 0, 
                           data_d['Fc_syn'] * stats_d['noct_slope'] + stats_d['noct_intcpt'],
                           data_d['Fc_syn'] * stats_d['day_slope'] + stats_d['day_intcpt'])     

    error = (np.random.laplace(0, sigma_delta) / 2).astype(np.float32)

    return data_d['Fc_syn'] + error
    
#------------------------------------------------------------------------------

def TRF(data_d,Eo,rb):

    return rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['Ta'] + 46.02)))

def LRF(data_d,rb,alpha,Aopt,k):
    
    Aopt_VPD = Aopt * np.exp(-k * (data_d['VPD'] - 1))
    index = np.where(data_d['VPD'] <= 1)[0]
    Aopt_VPD[index] = Aopt
    GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + (alpha * data_d['PAR'] / Aopt_VPD))
    index = np.where(data_d['PAR'] < 20)[0]
    GPP[index] = 0
    Reco = rb * np.exp(200 * (1 / (10 + 46.02) - 1 / (data_d['Ta'] + 46.02)))
    return GPP + Reco
    
#------------------------------------------------------------------------------

def optimise_Eo(data_d, var_to_fit):

    index = np.where(data_d['PAR'] < 20)[0]
    
    temp_d = {}
    for i in data_d.keys():
        temp_d[i] = data_d[i][index]

    drivers_d = {'Ta': temp_d['Ta']}
    response_d = temp_d[var_to_fit]

    params = curve_fit(TRF, drivers_d, response_d, p0 = [100, 1])

    return params[0]
    
def optimise_light(data_d, var_to_fit):
    
    index = np.where(data_d['PAR'] > 20)[0]    

    temp_d = {}
    for i in data_d.keys():
        temp_d[i] = data_d[i][index]

    drivers_d = {i: temp_d[i] for i in ['PAR','Ta','VPD']}
    response_d = temp_d[var_to_fit]

    params = curve_fit(LRF, drivers_d, response_d, p0 = [1, -1, -10, 1])
    
    return params[0]
    
#------------------------------------------------------------------------------
    
def pandasdf_to_dict(df):
    
    d = {}
    for i in df.columns:
        d[i] = np.array(df[i])
        
    return d

def main():
    
    df = pd.read_pickle('/home/imchugh/Analysis/Whroo/Data/Flux_and_met/synthetic_partitioning.df')
    
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
            


