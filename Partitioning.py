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


###### Constants ######

time_step_noct=5            #Number of steps in days between successive windows for nocturnal fitting
avg_window_noct=10          #Width of window for nocturnal fitting
time_step_day=5             #Number of steps in days between successive windows for daytime fitting
avg_window_day=15           #Width of window for daytime fitting
light_threshold=5           #Minimum light level for fitting
VPD_threshold=1             #Estimated (literature) threshold for VPD (kPa) stomatal response
Tref=10                     #Temperature reference for basal respiration (celsius)
D0=1                        #Threshold for stomatal response to VPD (kPa)
data_avail_threshold=50     #Percentage of data required in each year for fitting of Eo
data_period='30Min'         #Set the measurement period
temp_spread=5               #Minimum acceptable temperature (C) spread across sample for fitting
min_records=50              #Minimum number of acceptable records for fitting
min_pct_yr=30

###### Variable names ######

radName='Fsd_Con'           #Input name for solar radiation
radType='S'                 #Input type of measurement ('S' for full short-wave solar radiation specturm; 'PAR' for 400-700nm)
radUnits='Wm-2'             #Input units for incoming radiation ('Wm-2' or 'umol m-2 s-1')
tempName='Ta_Con'           #Input name for temperature
VWCName='Sws_Con'           #Input name for volumetric water content
VPDName='VPD_Con'           #Input name for vapour pressure deficit
VPDUnits='kPa'              #Input units for VPD
CfluxName='Fc'              #Input name for carbon flux data
CfluxUnits='umolCO2 m-2 s-1'#Input units for carbon flux data ('mgCO2 m-2 s-1' or 'umolCO2 m-2 s-1')
filterName='Fc_Con_QCFlag'  #Input name for filter variable in dataframe




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
def get_data(local_df,date_list,nocturnal_fit):
            
    for i in date_list:
                
        # Slice data from the complete dataframe
        sub_df=local_df.ix[i-dt.timedelta(days=avg_window_noct/2)+dt.timedelta(hours=t_del):
                           i+dt.timedelta(days=avg_window_noct/2)-dt.timedelta(hours=t_del)].dropna(axis=0,how='any')
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

    df,d=run()    
    
    pdb.set_trace()    
    
    ###### Housekeeping ######    
            
    # Globals
    global t_del, params_df, years_df
    
    # Trim for moving window (for even-numbered averaging windows)
    if avg_window_noct%2==False:
        t_del=12
    else:
        t_del=0
   
        
    ###### Create working arrays / dataframes ######
    
    # Date arrays for all days in time series and for centre points of moving windows
    date_array=np.array([dt.datetime.strftime(i,'%Y-%m-%d') for i in df.asfreq('D').index]) # Str array of all day dates in df
    valid_noct_dates=pd.to_datetime(pd.Series(date_array[avg_window_noct/2+1:len(date_array)-(avg_window_noct/2+1):time_step_noct])) # Datetime index of daytime dates to be sampled
    valid_day_dates=pd.to_datetime(pd.Series(date_array[avg_window_day/2+1:len(date_array)-(avg_window_day/2+1):time_step_day])) # Datetime index of nocturnal dates to be sampled
    
    # Yearly frequency dataframe with number of obs, number of cases and percentage of available data for each year
    years_df=pd.DataFrame({'N_obs':df[CfluxName].groupby([lambda x: x.year]).count()})
    years_df['N_recs']=[366 if i%4==0 else 365 for i in years_df.index]
    years_df['N_recs']=1440/int(data_period[:2])*years_df['N_recs']
    years_df['Data_avail_pct']=np.int8(np.float64(years_df['N_obs'])/years_df['N_recs']*100)
    years_df['process']=years_df['Data_avail_pct']>min_pct_yr
    if not np.any(years_df['process']):
        print 'Come back when you have more data... returning'
        return
    years_params=['Eo','rb','theta_1','theta_2']
    for i in years_params:
        years_df[i]=np.nan
    
    # Separate day and night data
    noct_df=df[[tempName,CfluxName,VWCName]][df[radName]<light_threshold].dropna(axis=0,how='any')
    day_df=df[[tempName,CfluxName,radName,VPDName,VWCName]][df[radName]>light_threshold].dropna(axis=0,how='any')
    
    # Parameters dataframe to contain fit parameters of temperature and light response functions
    params_df=pd.DataFrame({'rb_noct':np.nan,'rb_day':np.nan,'alpha':np.nan,'Aopt':np.nan,},index=pd.to_datetime(date_array))
    
    # Output dataframe to contain estimates of Re
    Fre_df=pd.DataFrame(index=df.index)
    
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
    get_data(noct_df,valid_noct_dates,nocturnal_fit)
    
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

# Fetch the data and prepare it for analysis
def run():
    
    # Prompt user for configuration file and get it
    root = Tkinter.Tk(); root.withdraw()
    cfName = tkFileDialog.askopenfilename(initialdir='')
    root.destroy()
    cf=ConfigObj(cfName)
    
    # Set input file and output path and create directories for plots and results
    file_in=os.path.join(cf['files']['input_path'],cf['files']['input_file'])
    path_out=cf['files']['output_path']
    plot_path_out=os.path.join(path_out,'Plots')
    if not os.path.isdir(plot_path_out): os.makedirs(os.path.join(path_out,'Plots'))
    results_path_out=os.path.join(path_out,'Results')
    if not os.path.isdir(results_path_out): os.makedirs(os.path.join(path_out,'Results'))    
    
    # Get user-set variable names from config file
    vars_data=[cf['variables']['data'][i] for i in cf['variables']['data']]
    vars_QC=[cf['variables']['QC'][i] for i in cf['variables']['QC']]
    vars_all=vars_data+vars_QC
       
    # Read .nc file
    nc_obj=netCDF4.Dataset(file_in)
    flux_frequency=int(nc_obj.time_step)
    dates_list=[dt.datetime(*xlrd.xldate_as_tuple(elem,0)) for elem in nc_obj.variables['xlDateTime']]
    d={}
    for i in vars_all:
        d[i]=nc_obj.variables[i][:]
    nc_obj.close()
    df=pd.DataFrame(d,index=dates_list)    
        
    # Build dictionary of additional configs
    d={}
    d['radiation_threshold']=int(cf['options']['radiation_threshold'])
    d['flux_frequency']=flux_frequency
    d['min_T_spread']=int(cf['options']['minimum_temperature_spread'])
    d['time_step_night']=int(cf['options']['time_step_night'])
    d['time_step_day']=int(cf['options']['time_step_day'])
    d['window_size_night']=int(cf['options']['window_size_night'])
    d['window_size_day']=int(cf['options']['window_size_day'])
    d['ref_T']=int(cf['options']['reference_temperature'])
    d['min_pct_data']=int(cf['options']['minimum_pct_annual_data'])
    d['min_num_data']=int(cf['options']['minimum_num_window_data'])
    d['ustar_filter']=float(cf['options']['ustar_filter'])
    
    pdb.set_trace()    
    
    if cf['options']['output_plots']=='True':
        d['plot_output_path']=plot_path_out
    if cf['options']['output_results']=='True':
        d['results_output_path']=results_path_out
        
    # Replace configured error values with NaNs and remove data with unacceptable QC codes, then drop flags
    df.replace(int(cf['options']['nan_value']),np.nan)
    if 'QC_accept_codes' in cf['options']:    
        QC_accept_codes=ast.literal_eval(cf['options']['QC_accept_codes'])
        eval_string='|'.join(['(df[vars_QC[i]]=='+str(i)+')' for i in QC_accept_codes])
        for i in xrange(4):
            df[vars_data[i]]=np.where(eval(eval_string),df[vars_data[i]],np.nan)
    df=df[vars_data]
    
    return df,d