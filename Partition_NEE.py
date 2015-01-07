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
def make_TRF_Eo(Tref):
    def TRF_Eo(Ta_S,rb,Eo):
        Reco=rb*np.exp(Eo*(1/(Tref+46.02)-1/(Ta_S+46.02)))
        return Reco
    return TRF_Eo

def make_TRF_rb(Tref,Eo):
    def TRF_rb(Ta_S,rb):
        Reco=rb*np.exp(Eo*(1/(Tref+46.02)-1/(Ta_S+46.02)))
        return Reco
    return TRF_rb

def make_LRF(Tref,Eo,D0):
    def LRF(local_df,alpha,Aopt,rb,k):
        Aopt_VPD=Aopt*np.exp(-k*(local_df[VPDName]-D0))
        index=np.where(local_df[VPDName]<=D0)[0]
        Aopt_VPD[index]=Aopt
        GPP=(alpha*local_df[radName])/(1-(local_df[radName]/2000)+(alpha*local_df[radName]/Aopt_VPD))
        Reco=rb*np.exp(Eo*(1/(Tref+46.02)-1/(local_df[tempName]+46.02)))
        return GPP+Reco
    return LRF

def est_Reco(local_df,params):
    return params['rb']*np.exp(params['Eo']*(1/(params['Tref']+46.02)-1/(local_df[tempName]+46.02)))

def est_GPP(local_df,params):
    Aopt_VPD=params['Aopt']*np.exp(-params['k']*(local_df[VPDName]-params['D0']))
    index=np.where(local_df[VPDName]<=params['D0'])[0]
    Aopt_VPD[index]=params['Aopt']
    GPP=(params['alpha']*local_df[radName])/(1-(local_df[radName]/2000)+(params['alpha']*local_df[radName]/Aopt_VPD))
    return GPP
#------------------------------------------------------------------------------
    
#------------------------------------------------------------------------------                                          
# Fetch the data and prepare it for analysis
def get_configs():
    
    # Prompt user for configuration file and get it
    root = Tkinter.Tk(); root.withdraw()
    cfName = tkFileDialog.askopenfilename(initialdir='')
    root.destroy()
    cf=ConfigObj(cfName)
    
    # Build dictionaries of config file contents
    d_paths={}        
    d_variables=dict(cf['variables'])
    d_options=dict(cf['options'])
    
    # Set input file and output path and create directories for plots and results
    d_paths['file_in']=os.path.join(cf['files']['input_path'],cf['files']['input_file'])

    if cf['options']['output_results']:
        results_output_path=os.path.join(cf['files']['output_path'],'Results')
        d_paths['results_output_path']=results_output_path
        if not os.path.isdir(results_output_path): os.makedirs(results_output_path)    
    if cf['options']['output_plots']:
        plot_output_path=os.path.join(cf['files']['output_path'],'Plots')
        d_paths['plot_output_path']=plot_output_path
        if not os.path.isdir(plot_output_path): os.makedirs(plot_output_path)
    
    # Prepare dictionary of user settings - drop strings or change to int / float
    for key in d_options:
        if d_options[key].isdigit():
            d_options[key]=int(d_options[key])
        else:
            try:
                d_options[key]=float(d_options[key])
            except ValueError:
                continue
    
    return d_paths,d_variables,d_options
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Open file and retrieve data
def get_data(file_in,d_variables,d_options):
    
    if d_options['infile_type']=='OzFluxQC':
        
        # Get user-set variable names from config file
        vars_data=d_variables['data'].values()
        vars_QC=d_variables['QC'].values()
        vars_all=vars_data+vars_QC
        
        # Read .nc file (write flux frequency to user options dictionary)
        nc_obj=netCDF4.Dataset(file_in)
        flux_freq=int(nc_obj.time_step)
        if flux_freq<>d_options['flux_frequency']:
            print 'Specified flux frequency does not match reported value in netCDF globals'
            return
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
        
        # Replace configured error values with NaNs and remove data with unacceptable QC codes, then drop QC flag variables
        df.replace(float(d_options['nan_value']),np.nan,inplace=True)
        if 'QC_accept_codes' in d_options:    
            QC_accept_codes=ast.literal_eval(d_options['QC_accept_codes'])
            d_options.pop('QC_accept_codes')
            eval_string='|'.join(['(df[vars_QC[i]]=='+str(i)+')' for i in QC_accept_codes])
            for i in xrange(4):
                df[vars_data[i]]=np.where(eval(eval_string),df[vars_data[i]],np.nan)
        df=df[vars_data]
    
    elif d_options['infile_type']=='Other':
        
        df=pd.read_pickle(file_in)
    
    return df
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Create a date index matching the specified step
def get_dates(df,d_options):
   
    night_win=d_options['window_size_night']
    night_step=d_options['time_step_night']
    day_win=d_options['window_size_day']
    day_step=d_options['time_step_day']    
    
    # Create datetime objects for valid dates and place into dictionary
    day_dates=pd.to_datetime([dt.datetime.strftime(i,'%Y-%m-%d') for i in df.asfreq('D').index])
    step_dates={'noct':day_dates[night_win/2+1:len(day_dates)-(night_win/2+1):night_step],
             'day':day_dates[day_win/2+1:len(day_dates)-(day_win/2+1):day_step]}   
    
    return day_dates,step_dates
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Program control
def main():

    # Globals
    global CfluxName,tempName,VWCName,radName,VPDName

    # Get data and user-specified options and configurations
    d_paths,d_variables,d_options=get_configs()
    
    df=get_data(d_paths['file_in'],d_variables,d_options)
    
    if not isinstance(df,pd.DataFrame):
        print 'Could not import data with current configuration settings... exiting'
        return

    # QC of raw data has been done, so drop the variables names for QC from the dictionary 
    d_variables=d_variables['data']

    # Set global variables
    CfluxName=d_variables['carbon_flux']
    tempName=d_variables['temperature']
    radName=d_variables['solar_radiation']
    VPDName=d_variables['vapour_pressure_deficit']
    
    # Create all day and stepped datetime indices
    all_dates,step_dates=get_dates(df,d_options)

    # Create dataframes for subsequent analysis    
    noct_df,day_df,years_df,params_df=make_dfs(df,d_options,d_variables,all_dates)

    # Optimisation 
    optimise(noct_df,years_df,params_df,step_dates['noct'],d_options,d_paths,True)
    optimise(day_df,years_df,params_df,step_dates['day'],d_options,d_paths,False)
    
    pdb.set_trace()    
    
    # Basic QC and interpolation
    QC_and_fill(params_df,['rb_noct','rb_day','alpha','Aopt','k'])
   
    # Calculate Re and GPP    
    output_df=resp(df, params_df)    
    
    # Output parameters file
    params_df.to_csv(os.path.join(d_paths['results_output_path'],'Fit_parameters.csv'))
    print 'Analysis complete'    
    
    return output_df,params_df
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Prepare dataframes
def make_dfs(df,d_options,d_variables,day_dates):
    
    # Retrieve objects from dictionaries
    flux_freq=d_options['flux_frequency']
    min_data=d_options['minimum_pct_annual_data']
    light_threshold=d_options['radiation_threshold']
    
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
    if d_options['convert_to_PAR']=='True':
        df[radName]=df[radName]*0.46*4.6 #Convert to PAR
    if d_options['convert_to_umol']=='True':
        df[CfluxName]=df[CfluxName]*1000/44.0 #Convert to umol
    noct_df=df[[tempName,CfluxName]][df[radName]<light_threshold].dropna(axis=0,how='any')
    day_df=df[[tempName,CfluxName,radName,VPDName]][df[radName]>light_threshold].dropna(axis=0,how='any')
     
    # Create dataframes for fit parameters and results
    params_df=pd.DataFrame({'rb_noct':np.nan,'rb_day':np.nan,'alpha':np.nan,'Aopt':np.nan,'Eo':np.nan,
                            'k':np.nan,'Status_noct':'Interpolated','Status_day':'Interpolated',
                            'Tref':d_options['reference_temperature'],'D0':d_options['VPD_threshold']},index=day_dates)
    
    return noct_df,day_df,years_df,params_df
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Coordinate the daytime optimisation
def optimise(df,years_df,params_df,dates,d_options,d_paths,noct_flag):

    # Retrieve objects from dictionaries
    window=d_options['window_size_day']
    min_n=d_options['minimum_num_window_data']
    temp_spread=d_options['minimum_temperature_spread']
    output_plots=d_options['output_plots']=='True'
    Tref=d_options['reference_temperature']
    
    # Trim for moving window (for even-numbered averaging windows)
    time_del=12 if window%2==0 else 0
    
    # Get number of windows that can potentially contain data
    num_avail_fits=len(dates)   

    # Do nocturnal optimisation to get annual Eo parameter    
    if noct_flag:
        
        print 'Calculating long term (annual) temperature sensitivity (Eo) for nocturnal data...'
        for i in years_df[years_df['process']].index:
            rb,Eo=curve_fit(make_TRF_Eo(Tref),
                            df[tempName].ix[str(i)],
                            df[CfluxName].ix[str(i)],
                            p0=[10,200])[0]
            years_df['Eo'].ix[i]=Eo
        mean=years_df['Eo'][years_df['process']].mean()
        years_df['Eo']=np.where(years_df['process'],years_df['Eo'],mean)
        print 'Done!'        
    
        # Broadcast the annual parameter value to the parameters df
        for i in years_df.index:
            params_df['Eo'].ix[str(i)]=years_df['Eo'].ix[i]    
    
    if noct_flag:
        print 'Calculating temperature response function parameter rb for each nocturnal data window...'
    else:
        print 'Calculating light (alpha, Aopt), VPD (k) and temperature response (rb) parameters for each daytime data window...'   
    
    if output_plots:    
        print 'Plotting valid fits...'
    else:
        print 'Plotting currently switched off... skipping plots'        

    # Do daytime optimisation to get light, VPD and temperature response parameters
    for i in dates:
                
        # Slice data from the complete dataframe
        sub_df=df.ix[i-dt.timedelta(days=window/2)+dt.timedelta(hours=time_del):
                     i+dt.timedelta(days=window/2)-dt.timedelta(hours=time_del)].dropna(axis=0,how='any')
             
        if noct_flag:            
            # If either too few data or temperature range is less than 5C, abort optimisation
            if len(sub_df)>=min_n and sub_df[tempName].max()-sub_df[tempName].min()>=temp_spread:

                 # Try nocturnal step optimisation - if causes error return nan array    
                try:
                    params_df['rb_noct'].ix[i]=curve_fit(make_TRF_rb(params_df['Tref'].ix[i],params_df['Eo'].ix[i]),
                                                         sub_df[tempName],
                                                         sub_df[CfluxName],
                                                         p0=1)[0]
                    params_df['Status_noct'].ix[i]='Valid fit'
                    plot_flag=True if output_plots else False
                except RuntimeError:
                    params_df['rb_noct'].ix[i]=np.nan
                    params_df['Status_noct'].ix[i]='Invalid fit'
                    plot_flag=False
                    
                # Plot if required otherwise skip
                if plot_flag:
                    # Create dictionary with file parameters
                    params_d=dict(params_df[['rb_noct','Eo','Tref']].ix[i])
                    params_d['rb']=params_d.pop('rb_noct')
                    plot_windows(sub_df,d_paths,params_d,window,i)

            else:
        
                if len(sub_df)<min_n:
                    params_df['Status_noct'][i]='Insufficient data'
                else:
                    params_df['Status_noct'][i]='Insufficient temperature spread'   
        
        else:
            # If either too few data or temperature range is less than 5C, abort optimisation
            if len(sub_df)>=min_n and sub_df[tempName].max()-sub_df[tempName].min()>=temp_spread:                
            
                # Try daytime step optimisation - if causes error return nan array    
                try:
                    params=curve_fit(make_LRF(params_df['Tref'].ix[i],params_df['Eo'].ix[i],params_df['D0'].ix[i]),
                                              sub_df[[radName,tempName,VPDName]],
                                              sub_df[CfluxName],
                                              p0=[-0.1,-10,1,1])[0]
                    params_df['Status_day'][i]='Valid fit'
                    plot_flag=True if output_plots else False
                except RuntimeError:
                    params=[np.nan,np.nan,np.nan,np.nan]
                    params_df['Status_day'].ix[i]='Invalid fit'
                    plot_flag=False
                params_df['alpha'][i]=params[0]
                params_df['Aopt'][i]=params[1]
                params_df['rb_day'][i]=params[2]
                params_df['k'][i]=params[3]
            
                # Plot if required otherwise skip
                if plot_flag:
                    # Create dictionary with file parameters
                    params_d={'Reco':dict(params_df[['rb_day','Eo','Tref']].ix[i]),
                              'GPP':dict(params_df[['Aopt','alpha','k','D0']].ix[i])}
                    params_d['Reco']['rb']=params_d['Reco'].pop('rb_day')
                    plot_windows(sub_df,d_paths,params_d,window,i)

            else:
        
                if len(sub_df)<min_n:
                    params_df['Status_day'][i]='Insufficient data'
                else:
                    params_df['Status_day'][i]='Insufficient temperature spread' 

    if noct_flag:
        num_actual_fits=len(params_df[params_df['Status_noct']=='Valid fit'])
    else:
        num_actual_fits=len(params_df[params_df['Status_day']=='Valid fit'])
    print 'Done!\nParameter optimisation was successful for '+str(num_actual_fits)+' of '+str(num_avail_fits)+' available fit windows'        
#------------------------------------------------------------------------------   

#------------------------------------------------------------------------------
# Plot daytime and nocturnal fits for each window
def plot_windows(df,d_paths,params_d,window,date):

    # Set parameters from dicts
    path=d_paths['plot_output_path']
    
    # Check whether plotting day or night and configure appropriately below
    daynight_ind='day' if 'GPP' in params_d.keys() else 'noct'
        
    # Calculate model estimated value
    if daynight_ind=='noct':
        df['NEE_est']=est_Reco(df,params_d)
        x_var=tempName
        x_lab=r'Temperature ($^{o}C$)'
    else:
        df['NEE_est']=est_Reco(df,params_d['Reco'])+est_GPP(df,params_d['GPP'])
        x_var=radName
        x_lab=r'PAR ($\mu mol\/photons\/m^{-2}s^{-1}$)'
    
    # Plot 
    date_str=dt.datetime.strftime(date,'%Y-%m-%d')
    fig=plt.figure(figsize=(12,8))
    fig.patch.set_facecolor('white')
    plt.plot(df[x_var],df[CfluxName],'bo')
    plt.plot(df[x_var],df['NEE_est'],'ro')
    plt.title('Fit for '+str(window)+' day window centred on '+date_str+'\n',fontsize=22)
    plt.xlabel(x_lab,fontsize=16)
    plt.ylabel(r'Fc ($\mu mol C\/m^{-2} s^{-1}$)',fontsize=16)
    plt.axhline(y=0,color='black')
    plot_out_name=daynight_ind+'_'+date_str+'.jpg'
    fig.savefig(os.path.join(path,plot_out_name))
    plt.close(fig)
    
    return
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Clean up the data - this needs improvement (specifically fit statistics)
def QC_and_fill(params_df,var_list):
    
    print 'Basic QC:'
    for i in var_list:
        IQR=np.percentile(params_df[i].dropna(),75)-np.percentile(params_df[i].dropna(),25)
        len_before=len(params_df[i].dropna())
        params_df[i]=np.where((params_df[i]<np.percentile(params_df[i].dropna(),25)-2*IQR)|
                              (params_df[i]>np.percentile(params_df[i].dropna(),75)+2*IQR),
                              np.nan,
                              params_df[i])
        len_after=len(params_df[i].dropna())
        print '          '+str(len_before-len_after)+' values removed from variable '+i
        params_df[i]=params_df[i].interpolate()
        params_df[i].fillna(method='ffill',inplace=True)
        params_df[i].fillna(method='bfill',inplace=True)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Calculate respiration using derived parameters and observed drivers
def resp(df,params_df):
    
    # Create output dataframe
    output_df=pd.DataFrame()
    output_df['NEE_measured']=df[CfluxName]    
    
    # Calculate Reco for time series using nocturnal rb estimates
    params_df.rename(columns={'rb_noct':'rb'},inplace=True)
    output_df['Reco_noct']=pd.concat([est_Reco(df.ix[i:i+dt.timedelta(hours=23.5)],
                                              params_df[['rb','Eo','Tref']].ix[i]) 
                                     for i in params_df.index])
    params_df.rename(columns={'rb':'rb_noct'},inplace=True)
    
    # Calculate GPP by difference using Reco_noct
    output_df['GPP (NEE-Reco_noct)']=output_df['NEE_measured']-output_df['Reco_noct']    
    
    # Calculate Re for time series for night estimates
    params_df.rename(columns={'rb_day':'rb'},inplace=True)
    output_df['Reco_day']=pd.concat([est_Reco(df.ix[i:i+dt.timedelta(hours=23.5)],
                                              params_df[['rb','Eo','Tref']].ix[i]) 
                                    for i in params_df.index])
    params_df.rename(columns={'rb':'rb_day'},inplace=True)
    
    # Calculate GPP by difference using Reco_noct
    output_df['GPP (NEE-Reco_day)']=output_df['NEE_measured']-output_df['Reco_day']    
    
    # Calculate parameter-based GPP
    output_df['GPP_parameter_based']=pd.concat([est_GPP(df.ix[i:i+dt.timedelta(hours=23.5)],
                                                        params_df[['Aopt','alpha','k','D0']].ix[i])
                                                for i in params_df.index])
    
    return output_df
#------------------------------------------------------------------------------