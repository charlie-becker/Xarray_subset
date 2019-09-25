# AUTHOR: CHALRIE BECKER
# DATE: 06/13/2019

import xarray as xr
import numpy as np
from varlist import var_list
import time
import glob
import matplotlib.pyplot as plt

t0 = time.time()

# Calculate GDD (Growing Degree Days)
def calc_GDD(tmin, tmax):
    
    adjusted_tmin = np.where(tmin <= 283.15, 283.15, tmin)
    adjusted_tmax = np.where(tmax >= 303.15, 303.15, tmax)
    adjusted_tmax = np.where(adjusted_tmax < 283.15, 283.15, adjusted_tmax)
    GDD = (adjusted_tmax+adjusted_tmin)/2 - 283.15
    return GDD

# Decummulate Precipitation 
def calc_precip(cum_precip, bucket_precip):
    
    total_precip = cum_precip + bucket_precip * 100.0
    PRCP = np.zeros(total_precip.shape)
    
    for i in np.arange(1,PRCP.shape[0]):
        
        PRCP[i,:,:] = total_precip[i,:,:].values - total_precip[i-1,:,:].values
        
    return PRCP

# Calculate Frost Hours
def calc_frost_hours(temp):
    
    frost_hours = (temp < 273.15).resample(XTIME = '24H').sum(axis = 0)
    
    return frost_hours 

######### MAIN #################

for year in np.arange(1988, 1990):
    
    # Select all files except last (which is a repeat of the first hour of the following water year) for each WY
    files = sorted(glob.glob('/mnt/wrf_history/vol??/wrf_out/wy_' + str(year) + '/d02/wrfout_d02_*'))[:-1]
    
    # Open multi-file data
    d = xr.open_mfdataset(files, drop_variables = var_list, concat_dim = 'Time', parallel = True)
    
    # Swap Time dimenstions to use 'resample'
    d = d.swap_dims({'Time':'XTIME'})

    # Daily Aggregations (hourly to daily)
    d['PRCP'] = d['RAINNC']
    d['PRCP'].values = calc_precip(d['RAINNC'],d['I_RAINNC'])
    new_array = d[['SWDOWN','SWNORM','Q2','T2']].resample(XTIME = '24H').mean(dim = 'XTIME') # create daily means of few variables
    new_array['TMIN'] = d['T2'].resample(XTIME = '24H').min(dim = 'XTIME') # create daily minimum temperature
    new_array['TMAX'] = d['T2'].resample(XTIME = '24H').max(dim = 'XTIME')  # create daily maximum temperature
    new_array['PRCP'] = d['PRCP'].resample(XTIME = '24H').sum(dim = 'XTIME')
    new_array['GDD'] = new_array['TMIN']
    new_array['GDD'].values = calc_GDD(new_array['TMIN'],new_array['TMAX'])
    new_array['FROSTH'] = calc_frost_hours(d['T2'])

    # rename T2 as TMEAN
    new_array = new_array.rename({'T2' : 'TMEAN'}) 

    # Adjust some meta data
    new_array['TMEAN'].attrs = [('description','DAILY MEAN GRID SCALE TEMPERATUTE'), ('units','K')]
    new_array['TMIN'].attrs = [('description','DAILY MINIMUM GRID SCALE TEMPERATURE'), ('units','K')]
    new_array['TMAX'].attrs = [('description','DAILY MAXIMUM GRID SCALE TEMPERATURE'), ('units','K')]
    new_array['Q2'].attrs = [('description','DAILY MEAN GRID SCALE SPECIFIC HUMIDITY'), ('units','')]
    new_array['SWDOWN'].attrs = [('description','DAILY MEAN DOWNWARD SHORT WAVE FLUX AT GROUND SURFACE'), ('units','W m^2')]
    new_array['SWNORM'].attrs = [('description','DAILY MEAN NORMAL SHORT WAVE FLUX AT GROUND SURFACE (SLOPE-DEPENDENT)'), ('units','W m^2')]
    new_array['GDD'].attrs = [('description','DAILY GROWING DEGREE DAYS'), ('units','')]
    new_array['PRCP'].attrs = [('description','DAILY ACCUMULATED GRID SCALE PRECIPITATION'), ('units','mm')]
    new_array['FROSTH'].attrs = [('description','DAILY HOURS BELOW FREEZING (< 273.15 K)'), ('units','h')]

    # Write new netcdf file
    new_array.to_netcdf("/mnt/selway/data/data_02/charlie/subsets/WRF30YR/d02/DAILY_30YR_D02_WY" + str(year) + ".nc")
    
    print("Successfully wrote: /mnt/selway/data/data_02/charlie/subsets/WRF30YR/d02/DAILY_30YR_D02_WY" + str(year) + ".nc")
    del d, new_array
    
t1 = time.time()
print("Total time to create this subset was:", t1 - t0, "seconds.")