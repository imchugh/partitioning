import Tkinter, tkFileDialog
from configobj import ConfigObj
import netCDF4
import xlrd
import ast
import os
import pandas as pd
import numpy as np
import datetime as dt
from scipy.optimize import curve_fit
import pdb


###### Functions ######

def TRF_Eo(local_df,rb,Eo,theta_1,theta_2):
    f_VWC=1#1/(1+np.exp(theta_1-theta_2*local_df[VWCName]))
    Reco=rb*np.exp(Eo*(1/(Tref+46.02)-1/(local_df[tempName]+46.02)))*f_VWC
    return Reco

def TRF_rb(local_df,rb):
    f_VWC=1#1/(1+np.exp(theta_1-theta_2*local_df[VWCName]))
    Reco=rb*np.exp(Eo*(1/(Tref+46.02)-1/(local_df[tempName]+46.02)))*f_VWC
    return Reco

def LRF(local_df,alpha,Aopt,rb):
    #Aopt=Aopt_0*np.exp(-k*(local_df[VPDName]-D0))
    #index=np.where(local_df[VPDName]<=D0)[0]
    #Aopt[index]=Aopt_0
    GPP=(alpha*local_df[radName])/(1-(local_df[radName]/2000)+(alpha*local_df[radName]/Aopt))
    f_VWC=1#1/(1+np.exp(theta_1-theta_2*local_df[VWCName]))
    Reco=rb*np.exp(Eo*(1/(Tref+46.02)-1/(local_df[tempName]+46.02)))*f_VWC
    return GPP+Reco

def GPP_func(local_df,alpha,Aopt):
    return (alpha*local_df)/(1-(local_df/2000)+(alpha*local_df/Aopt))    
    
# Do the optimisation (handle errors) then return parameters
def get_data(local_df,date_list,nocturnal_fit,time_del):
            
    for i in date_list:
                
        # Slice data from the complete dataframe
        sub_df=local_df.ix[i-dt.timedelta(days=avg_window_noct/2)+dt.timedelta(hours=time_del):
                           i+dt.timedelta(days=avg_window_noct/2)-dt.timedelta(hours=time_del)].dropna(axis=0,how='any')
        yr=sub_df.index[0].year
        global Eo, theta_1, theta_2
        Eo=years_df['Eo'].ix[yr]
        theta_1=years_df['theta_1'].ix[yr]
        theta_2=years_df['theta_2'].ix[yr]
        
        # If either too few data or temperature range is less than 5C, abort optimisation
        if len(sub_df)>=min_records and sub_df[tempName].max()-sub_df[tempName].min()>=temp_spread:
            
            # Try optimisation - if causes error return nan array    
            if nocturnal_fit==True:
                try:
                    params_df['rb_noct'][i]=curve_fit(TRF_rb,sub_df[[tempName,VWCName]],sub_df[CfluxName],p0=1)[0]
                except RuntimeError:
                    params_df['rb_noct'][i]=np.nan
            else:
                try:
                    a=curve_fit(LRF,sub_df[[radName,tempName,VPDName,VWCName]],sub_df[CfluxName],p0=[-0.1,-10,1])[0]
                except RuntimeError:
                    a=[np.nan,np.nan,np.nan]
                params_df['alpha'][i]=a[0]
                params_df['Aopt'][i]=a[1]
                #params_df['k'][i]=a[2]
                params_df['rb_day'][i]=a[2]                        


           
#------------------------#                                    
###### Main program ######

def main():

    # Globals
    global params_df, years_df

    # Get data and user-specified options and configurations
    df,d_paths,d_variables,d_options=run()    
   
    # Trim for moving window (for even-numbered averaging windows)
    time_del=12 if d_options['window_size_night']%2==0 else 0
    
    # Create datetime indices and dataframes for subsequent analysis
    (date_array,
     noct_dates,
     day_dates,
     years_df,
     noct_df,
     day_df,
     params_df,
     output_df)=prep(df,d_options,d_variables)
    
    ###### Nocturnal optimisation (rb, Eo, theta_1, theta_2) ######
    
    nocturnal_fit=True
    
    print 'Calculating temperature sensitivity (Eo) and soil moisture-dependent sigmoid modification parameters (theta_1 and theta_2) across all nocturnal data...'
    for i in years_df[years_df['process']].index:
        rb,Eo,theta_1,theta_2=curve_fit(TRF_Eo,noct_df[[tempName,VWCName]].ix[str(i)],noct_df[CfluxName].ix[str(i)],p0=[10,200,1,10])[0]
        years_df['rb'].ix[i]=rb        
        years_df['Eo'].ix[i]=Eo
        years_df['theta_1'].ix[i]=theta_1
        years_df['theta_2'].ix[i]=theta_2        
    for i in ['rb','Eo','theta_1','theta_2']:
        mean=years_df[i][years_df['process']].mean()
        years_df[i]=np.where(years_df['process'],years_df[i],mean)
   
    print 'Calculating temperature response function parameter rb for each window'
    get_data(noct_df,valid_noct_dates,nocturnal_fit,time_del)
    
    ####### Daytime optimisation (Aopt,alpha,rb,k) ######
    
    nocturnal_fit=False
    
    # Fill the parameters dataframe alpha, Aopt and rb_day parameters ######
    print 'Calculating combined light/temperature response function parameters alpha, Aopt and rb for each window'
    get_data(day_df,valid_day_dates,nocturnal_fit)

    ###### Clean up data ######
    for i in params_df.columns:
        IQR=np.percentile(params_df[i].dropna(),75)-np.percentile(params_df[i].dropna(),25)
        len_before=len(params_df[i].dropna())
        params_df[i]=np.where((params_df[i]<np.percentile(params_df[i].dropna(),25)-2*IQR)|
                              (params_df[i]>np.percentile(params_df[i].dropna(),75)+2*IQR),
                              np.nan,
                              params_df[i])
        len_after=len(params_df[i].dropna())
        print str(len_before-len_after)+' values removed from variable '+i+' during QC'
        params_df[i+'_interp']=params_df[i].interpolate()
    
                    
    ###### Calculate and output Re and GPP for time series ######
    
    # Calculate Re for time series for day and night estimates
    Fre_df['Fre_noct']=pd.concat([TRF_Eo(df[[tempName,VWCName]].ix[i:i+dt.timedelta(hours=23.5)],
                                         params_df['rb_noct_interp'][i],Eo,theta_1,theta_2) for i in params_df.index])
    Fre_df['Fre_day']=pd.concat([TRF_Eo(df[[tempName,VWCName]].ix[i:i+dt.timedelta(hours=23.5)],
                                        params_df['rb_day_interp'][i],Eo,theta_1,theta_2) for i in params_df.index])

    # Calculate parameter-based GPP for daytime                                                                            
    Fre_df['GPP_Lasslop']=pd.concat([GPP_func(df[radName].ix[i:i+dt.timedelta(hours=23.5)],
                                              params_df['alpha_interp'][i],params_df['Aopt_interp'][i]) for i in params_df.index])

    return params_df

#    # Force nocturnal GPP to 0
#    Fre_df['GPP_Lasslop'][df[radName]<light_threshold]=0     
#            
#            
#    ###### Clean up and export ######
#    
#    return params_df
#    
#    # Pull together        
#    ALL_combined=df.join(Fre_df,how='left')   
#        
#    # Reinstate original variable values if changed
#    if radType=='S' or radUnits=='Wm-2':
#        df[radName]=tempRad_S
#    
#    if CfluxUnits=='mgCO2 m-2 s-1':
#        df[CfluxName]=tempFc_S												
#    									
#    return ALL_combined

def prep(df,d_opts,d_vars):
    
    night_win=d_opts['window_size_night'], night_step=d_opts['time_step_night']
    day_win=d_opts['window_size_day'], day_step=d_opts['time_step_day']
    flux_freq=d_opts['flux_frequency'], min_data=d_opts['minimum_pct_annual_data'], light_threshold=d_opts['radiation_threshold']
    CfluxName=d_vars['carbon_flux'], tempName=d_vars['temperature'], VWCName=d_vars['soil_water']
    radName=d_vars['solar_radiation'], VPDName=d_vars['vapour_pressure_deficit']
    
    # Create datetime objects for valid dates
    date_array=np.array([dt.datetime.strftime(i,'%Y-%m-%d') for i in df.asfreq('D').index])
    noct_dates=pd.to_datetime(date_array[night_win/2+1:len(date_array)-(night_win/2+1):night_step])
    day_dates=pd.to_datetime(date_array[day_win/2+1:len(date_array)-(day_win/2+1):day_step])
    #noct_dates=pd.to_datetime(pd.Series(date_array[night_win/2+1:len(date_array)-(night_win/2+1):night_step]))
    #day_dates=pd.to_datetime(pd.Series(date_array[day_win/2+1:len(date_array)-(day_win/2+1):day_step]))
   
    # Yearly frequency dataframe with number of obs, number of cases and percentage of available data for each year
    years_df=pd.DataFrame({'N_obs':df[CfluxName].groupby([lambda x: x.year]).count()})
    years_df['N_recs']=[366 if i%4==0 else 365 for i in years_df.index]
    years_df['N_recs']=1440/flux_freq*years_df['N_recs']
    years_df['Data_avail_pct']=np.int8(np.float64(years_df['N_obs'])/years_df['N_recs']*100)
    years_df['process']=years_df['Data_avail_pct']>min_data
    if not np.any(years_df['process']):
        print 'Come back when you have more data... returning'
        return
    years_params=['Eo','rb','theta_1','theta_2']
    for i in years_params:
        years_df[i]=np.nan
   
    # Separate day and night data
    noct_df=df[[tempName,CfluxName,VWCName]][df[radName]<light_threshold].dropna(axis=0,how='any')
    day_df=df[[tempName,CfluxName,radName,VPDName,VWCName]][df[radName]>light_threshold].dropna(axis=0,how='any')
    
    # Create dataframes for fit parameters and results
    params_df=pd.DataFrame({'rb_noct':np.nan,'rb_day':np.nan,'alpha':np.nan,'Aopt':np.nan,},index=pd.to_datetime(date_array))
    output_df=pd.DataFrame(index=df.index)

    pdb.set_trace() 
   
    return date_array,noct_dates,day_dates,years_df,noct_df,day_df,params_df,output_df

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