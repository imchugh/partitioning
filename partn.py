import Tkinter, tkFileDialog
from configobj import ConfigObj
import matplotlib.pyplot as plt
import netCDF4
import xlrd
import ast
import os
import pandas as pd
import numpy as np
import datetime as dt
from scipy.optimize import curve_fit
import pdb



#------------------------------------------------------------------------------
# Response functions

def TRF_Eo(local_df,rb,Eo):
    #f_VWC=1/(1+np.exp(theta_1-theta_2*local_df[VWCName]))
    Reco=rb*np.exp(Eo*(1/(Tref+46.02)-1/(local_df[tempName]+46.02)))
    return Reco

def TRF_rb(local_df,rb):
    #f_VWC=1/(1+np.exp(theta_1-theta_2*local_df[VWCName]))
    Reco=rb*np.exp(Eo*(1/(Tref+46.02)-1/(local_df[tempName]+46.02)))
    return Reco

def LRF(local_df,alpha,Aopt,rb,k):
    Aopt_VPD=Aopt*np.exp(-k*(local_df[VPDName]-D0))
    index=np.where(local_df[VPDName]<=D0)[0]
    Aopt_VPD[index]=Aopt
    GPP=(alpha*local_df[radName])/(1-(local_df[radName]/2000)+(alpha*local_df[radName]/Aopt_VPD))
    Reco=rb*np.exp(Eo*(1/(Tref+46.02)-1/(local_df[tempName]+46.02)))
    #f_VWC=1/(1+np.exp(theta_1-theta_2*local_df[VWCName]))
    return GPP+Reco

def GPP_func(local_df,alpha,Aopt):
    return (alpha*local_df)/(1-(local_df/2000)+(alpha*local_df/Aopt))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------    
# Loop through available data
def get_data(df,date_list,nocturnal_fit,window,min_n,temp_spread):

    # Set annual values for respiration parameters to global        
    global Eo
    
    # Trim for moving window (for even-numbered averaging windows)
    time_del=12 if window%2==0 else 0
        
    for i in date_list:
                
        # Slice data from the complete dataframe
        sub_df=df.ix[i-dt.timedelta(days=window/2)+dt.timedelta(hours=time_del):
                     i+dt.timedelta(days=window/2)-dt.timedelta(hours=time_del)].dropna(axis=0,how='any')
        
        # If either too few data or temperature range is less than 5C, abort optimisation
        if len(sub_df)>=min_n and sub_df[tempName].max()-sub_df[tempName].min()>=temp_spread:
            
            # Set parameters
            yr=sub_df.index[0].year
            Eo=years_df['Eo'].ix[yr]
            
            # Try optimisation - if causes error return nan array    
            if nocturnal_fit==True:
                try:
                    params_df['rb_noct'][i]=curve_fit(TRF_rb,sub_df[[tempName,VWCName]],sub_df[CfluxName],p0=1)[0]
                except RuntimeError:
                    params_df['rb_noct'][i]=np.nan
            else:
                try:
                    a=curve_fit(LRF,sub_df[[radName,tempName,VWCName,VPDName]],sub_df[CfluxName],p0=[-0.1,-10,1,1])[0]
                except RuntimeError:
                    a=[np.nan,np.nan,np.nan]
                params_df['alpha'][i]=a[0]
                params_df['Aopt'][i]=a[1]
                params_df['rb_day'][i]=a[2]                        
                params_df['k'][i]=a[3]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Program control
def main():

    # Globals
    global Tref,D0
    global params_df, years_df
    global CfluxName,tempName,VWCName,radName,VPDName

    # Get data and user-specified options and configurations
    df,d_paths,d_variables,d_options=run()    

    # Set global variables
    CfluxName=d_variables['carbon_flux']
    tempName=d_variables['temperature']
    VWCName=d_variables['soil_water']
    radName=d_variables['solar_radiation']
    VPDName=d_variables['vapour_pressure_deficit']
    Tref=d_options['reference_temperature']
    D0=d_options['VPD_threshold']
    
    # Create datetime indices and dataframes for subsequent analysis
    # Handle case of empty return!!!
    (d_dates,
     d_dfs,
     years_df,
     params_df,
     output_df)=prep(df,d_options,d_variables)
  
    # Optimisation 
    optimise(d_dfs,d_dates,d_options) 
    
    # Clean up
    QC()
    
    # Calculate Re and GPP    
    resp(df)    
    
    # Plot fits
    plot_windows(d_dfs,d_options,d_paths)

    # Plot interpolated parameter time series
    plot_paramater_ts(params_df)
    
    # Output parameters file
    params_df.to_csv(os.path.join(d_paths['results_output_path'],'Fit_parameters.csv'))
    print 'Analysis complete'    
    
    return params_df
#------------------------------------------------------------------------------
    
#------------------------------------------------------------------------------                   
# Coordinate the optimisation
def optimise(d_dfs,d_dates,d_opts):
    
    # Set nocturnal processing flag high
    nocturnal_fit=True
    
    # Retrieve objects from dictionaries
    noct_df=d_dfs['noct']
    noct_dates=d_dates['noct']
    noct_win=d_opts['window_size_night']
    day_df=d_dfs['day']
    day_dates=d_dates['day']
    day_win=d_opts['window_size_day']
    min_n=d_opts['minimum_num_window_data']
    temp_spread=d_opts['minimum_temperature_spread']
    
    # Do nocturnal optimisation to get annual Lloyd and Taylor parameters 
    print 'Calculating long term (annual) temperature sensitivity (Eo) for nocturnal data...'
    for i in years_df[years_df['process']].index:
        rb,Eo=curve_fit(TRF_Eo,noct_df[[tempName,VWCName]].ix[str(i)],noct_df[CfluxName].ix[str(i)],p0=[10,200])[0]
        years_df['Eo'].ix[i]=Eo
    mean=years_df['Eo'][years_df['process']].mean()
    years_df['Eo']=np.where(years_df['process'],years_df['Eo'],mean)
        
    # Broadcast the annual parameter values to the parameters df
    for i in years_df.index:
        params_df['Eo'].ix[str(i)]=years_df['Eo'].ix[i]
    
    # Do nocturnal optimisation using specified window to capture seasonally varying rb
    print 'Calculating temperature response function parameter rb for each window'
    get_data(noct_df,noct_dates,nocturnal_fit,noct_win,min_n,temp_spread)
    
    # Set daytime processing flag high
    nocturnal_fit=False
    
    # Fill the parameters dataframe alpha, Aopt and rb_day parameters ######
    print 'Calculating combined light/temperature response function parameters alpha, Aopt, rb and k for each window'
    get_data(day_df,day_dates,nocturnal_fit,day_win,min_n,temp_spread)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Plot daytime and nocturnal fits for each window
def plot_windows(d_dfs,d_opts,d_paths):
    
    # Set parameters from dicts
    window=d_opts['window_size_night']
    time_del=12 if window%2==0 else 0     
    path=d_paths['plot_output_path']
    noct_df=d_dfs['noct']
    day_df=d_dfs['day']
    
    # Do plotting for day and night fits
    print 'Plotting nocturnal and daytime response functions'        
    for key in d_dfs:
        df=d_dfs[key]
        daynight_ind='Status_noct' if key=='noct' else 'Status_day'
        for i in params_df[params_df[daynight_ind]=='Fitted_data'].index:
           
            # Slice data from the complete dataframe
            sub_df=df.ix[i-dt.timedelta(days=window/2)+dt.timedelta(hours=time_del):
                         i+dt.timedelta(days=window/2)-dt.timedelta(hours=time_del)].dropna(axis=0,how='any')
     
            # Specify value of global variable Eo and calculate response function
            Eo=params_df['Eo'].ix[i]
            if key=='noct':            
                sub_df['NEE_est']=TRF_rb(sub_df,params_df['rb_noct_interp'].ix[i])
                x_var=tempName
                x_lab=r'Temperature ($^{o}C$)'
            else:
                sub_df['NEE_est']=LRF(sub_df,params_df['alpha_interp'].ix[i],params_df['Aopt_interp'].ix[i],
                                      params_df['rb_day_interp'].ix[i],params_df['k'].ix[i])
                x_var=radName
                x_lab=r'PAR ($\mu mol\/photons\/m^{-2}s^{-1}$)'
            
            # Plot 
            date_str=dt.datetime.strftime(i,'%Y-%m-%d')
            fig=plt.figure(figsize=(12,8))
            fig.patch.set_facecolor('white')
            plt.plot(sub_df[x_var],sub_df[CfluxName],'bo')
            plt.plot(sub_df[x_var],sub_df['NEE_est'],'ro')
            plt.title('Fit for '+str(window)+' day window centred on '+date_str+'\n',fontsize=22)
            plt.xlabel(x_lab,fontsize=16)
            plt.ylabel(r'Fc ($\mu mol C\/m^{-2} s^{-1}$)',fontsize=16)
            plt.axhline(y=0,color='black')
            plot_out_name=key+'_'+date_str+'.jpg'
            fig.savefig(os.path.join(path,plot_out_name))
            plt.close(fig)
    
    return

#------------------------------------------------------------------------------
# Prepare dataframes and dictionaries
def prep(df,d_opts,d_vars):
    
    # Retrieve objects from dictionaries
    night_win=d_opts['window_size_night']
    night_step=d_opts['time_step_night']
    day_win=d_opts['window_size_day']
    day_step=d_opts['time_step_day']
    flux_freq=d_opts['flux_frequency']
    min_data=d_opts['minimum_pct_annual_data']
    light_threshold=d_opts['radiation_threshold']
    
    # Create datetime objects for valid dates and place into dictionary
    dates=pd.to_datetime([dt.datetime.strftime(i,'%Y-%m-%d') for i in df.asfreq('D').index])
    d_dates={'noct':dates[night_win/2+1:len(dates)-(night_win/2+1):night_step],
             'day':dates[day_win/2+1:len(dates)-(day_win/2+1):day_step]}
    
    # Yearly frequency dataframe with number of obs, number of cases and percentage of available data for each year
    years_df=pd.DataFrame({'N_obs':df[CfluxName].groupby([lambda x: x.year]).count()})
    years_df['N_recs']=[366 if i%4==0 else 365 for i in years_df.index]
    years_df['N_recs']=1440/flux_freq*years_df['N_recs']
    years_df['Data_avail_pct']=np.int8(np.float64(years_df['N_obs'])/years_df['N_recs']*100)
    years_df['process']=years_df['Data_avail_pct']>min_data
    if not np.any(years_df['process']):
        print 'Come back when you have more data... returning'
        return
    years_params=['Eo','rb']
    for i in years_params:
        years_df[i]=np.nan
   
    # Separate day and night data and place into dict
    d_dfs={}
    df[radName]=df[radName]*0.46*4.6 #Convert to PAR
    d_dfs['noct']=df[[tempName,CfluxName,VWCName]][df[radName]<light_threshold].dropna(axis=0,how='any')
    d_dfs['day']=df[[tempName,CfluxName,radName,VWCName,VPDName]][df[radName]>light_threshold].dropna(axis=0,how='any')
    
    
    # Create dataframes for fit parameters and results
    params_df=pd.DataFrame({'rb_noct':np.nan,'rb_day':np.nan,'alpha':np.nan,'Aopt':np.nan,'Eo':np.nan,
                            'k':np.nan,'Status_noct':np.nan,'Status_day':np.nan},index=dates)
    output_df=pd.DataFrame(index=df.index)
    
    return d_dates,d_dfs,years_df,params_df,output_df
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Clean up the data - this needs improvement (specifically fit statistics)
def QC():
    
    var_list=['rb_noct','rb_day','Aopt','alpha','k']
    for i in var_list:
        IQR=np.percentile(params_df[i].dropna(),75)-np.percentile(params_df[i].dropna(),25)
        len_before=len(params_df[i].dropna())
        params_df[i]=np.where((params_df[i]<np.percentile(params_df[i].dropna(),25)-2*IQR)|
                              (params_df[i]>np.percentile(params_df[i].dropna(),75)+2*IQR),
                              np.nan,
                              params_df[i])
        len_after=len(params_df[i].dropna())
        print str(len_before-len_after)+' values removed from variable '+i+' during QC'
        params_df[i+'_interp']=params_df[i].interpolate()
        params_df[i+'_interp'].fillna(method='ffill').fillna(method='bfill')
    params_df['Status_day']=np.where((np.isnan(params_df['rb_day']))|(np.isnan(params_df['Aopt']))|
                                     (np.isnan(params_df['alpha'])),
                                     'Interpolated',
                                     'Fitted_data')
    params_df['Status_noct']=np.where(np.isnan(params_df['rb_noct']),'Interpolated','Fitted_data')

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Calculate respiration using derived parameters and observed drivers
def resp(df):
    
    temp_df=pd.DataFrame()    
    
        # Calculate Re for time series for day and night estimates
    temp_df['Fre_noct']=pd.concat([TRF_Eo(df[[tempName,VWCName]].ix[i:i+dt.timedelta(hours=23.5)],
                                          params_df['rb_noct_interp'][i],Eo) for i in params_df.index])
    temp_df['Fre_day']=pd.concat([TRF_Eo(df[[tempName,VWCName]].ix[i:i+dt.timedelta(hours=23.5)],
                                         params_df['rb_day_interp'][i],Eo) for i in params_df.index])

    # Calculate parameter-based GPP for daytime                                                                            
    temp_df['GPP_Lasslop']=pd.concat([GPP_func(df[radName].ix[i:i+dt.timedelta(hours=23.5)],
                                               params_df['alpha_interp'][i],params_df['Aopt_interp'][i]) for i in params_df.index])
    
    return temp_df
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------                                          
# Fetch the data and prepare it for analysis
def run():
    
    # Prompt user for configuration file and get it
    root = Tkinter.Tk(); root.withdraw()
    cfName = tkFileDialog.askopenfilename(initialdir='')
    root.destroy()
    cf=ConfigObj(cfName)
    
    # Build dictionaries of config file contents
    d_paths={}        
    d_variables=dict(cf['variables']['data'])
    d_options=dict(cf['options'])
    
    # Set input file and output path and create directories for plots and results
    file_in=os.path.join(cf['files']['input_path'],cf['files']['input_file'])

    if cf['options']['output_results']:
        results_output_path=os.path.join(cf['files']['output_path'],'Results')
        d_paths['results_output_path']=results_output_path
        if not os.path.isdir(results_output_path): os.makedirs(results_output_path)    
    if cf['options']['output_plots']:
        plot_output_path=os.path.join(cf['files']['output_path'],'Plots')
        d_paths['plot_output_path']=plot_output_path
        if not os.path.isdir(plot_output_path): os.makedirs(plot_output_path)
        
    # Get user-set variable names from config file
    vars_data=cf['variables']['data'].values()
    vars_QC=cf['variables']['QC'].values()
    vars_all=vars_data+vars_QC
    
    # Read .nc file (write flux frequency to user options dictionary)
    nc_obj=netCDF4.Dataset(file_in)
    flux_frequency=int(nc_obj.time_step)
    dates_list=[dt.datetime(*xlrd.xldate_as_tuple(elem,0)) for elem in nc_obj.variables['xlDateTime']]
    d={}
    for i in vars_all:
        ndims=len(nc_obj.variables[i].shape)
        if ndims==3:
            d[i]=nc_obj.variables[i][:,0,0]
        elif ndims==1:    
            d[i]=nc_obj.variables[i][:]
    nc_obj.close()
    df=pd.DataFrame(d,index=dates_list)    
    
    # Replace configured error values with NaNs and remove data with unacceptable QC codes, then drop flags
    df.replace(float(d_options['nan_value']),np.nan,inplace=True)
    if 'QC_accept_codes' in d_options:    
        QC_accept_codes=ast.literal_eval(d_options['QC_accept_codes'])
        d_options.pop('QC_accept_codes')
        eval_string='|'.join(['(df[vars_QC[i]]=='+str(i)+')' for i in QC_accept_codes])
        for i in xrange(4):
            df[vars_data[i]]=np.where(eval(eval_string),df[vars_data[i]],np.nan)
    df=df[vars_data]
    
    # Prepare dictionary of user settings - drop strings or change to int / float
    for i in ['nan_value','output_results','output_plots']:
        d_options.pop(i,None)
    for key in d_options:
        if d_options[key].isdigit():
            d_options[key]=int(d_options[key])
        else:
            d_options[key]=float(d_options[key])
    d_options['flux_frequency']=flux_frequency
    
    return df,d_paths,d_variables,d_options
#------------------------------------------------------------------------------