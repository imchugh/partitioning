# -*- coding: utf-8 -*-
"""
Created on Wed May 27 09:28:52 2015

@author: imchugh
"""

import numpy as np
from scipy.optimize import curve_fit
import pdb

#------------------------------------------------------------------------------

def make_data(df,Eo,rb,alpha,Aopt,k):

    Aopt_VPD = Aopt * np.exp(-k * (df['VPD'] - 1))
    index = np.where(df['VPD'] <= 1)[0]
    Aopt_VPD[index] = Aopt
    GPP = (alpha * df['PAR']) / (1 - (df['PAR'] / 2000) + (alpha * df['PAR'] / Aopt_VPD))
    Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (df['Ta'] + 46.02)))
    df['Fc_syn'] = np.where(df['PAR'] > 20, GPP + Reco, Reco)

def random_error(data_df, stats_df):
    
    day_slope = stats_df.loc['day', 'slope']
    day_int = stats_df.loc['day', 'intcpt']
    noct_slope = stats_df.loc['noct', 'slope']
    noct_int = stats_df.loc['noct', 'intcpt']
    
    sigma_delta = np.where(data_df['Fc_syn'] > 0, 
                           data_df['Fc_syn'] * noct_slope + noct_int,
                           data_df['Fc_syn'] * day_slope + day_int)     

    error = np.random.laplace(0, sigma_delta) / 2

    data_df['Fc_syn_error'] = data_df['Fc_syn'] + error
    data_df['Fc_syn_error'] = data_df['Fc_syn_error'].astype(np.float32)
    
#------------------------------------------------------------------------------

def TRF(df,Eo,rb):

    return rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (df['Ta'] + 46.02)))

def LRF(df,rb,alpha,Aopt,k):
    
    Aopt_VPD = Aopt * np.exp(-k * (df['VPD'] - 1))
    index = np.where(df['VPD'] <= 1)[0]
    Aopt_VPD[index] = Aopt
    GPP = (alpha * df['PAR']) / (1 - (df['PAR'] / 2000) + (alpha * df['PAR'] / Aopt_VPD))
    Reco = rb * np.exp(200 * (1 / (10 + 46.02) - 1 / (df['Ta'] + 46.02)))
    return GPP + Reco
    
#------------------------------------------------------------------------------

def optimise_Eo(df, var_to_fit):
    
    df = df[df.PAR < 20]
    
    params = curve_fit(TRF, df, df[var_to_fit], p0 = [200, 1])
        
    return params[0]
    
def optimise_light(df, var_to_fit):
    
    df = df[df.PAR > 20]
    print len(df)
    params = curve_fit(LRF, df, df[var_to_fit], p0 = [1, -0.1, -10, 0.5])
    
    return params[0]
    
#------------------------------------------------------------------------------