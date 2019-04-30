
import xarray as xr
import numpy as np
from varlist import var_list
import time

t0 = time.time()

for year in np.arange(1987, 2018):
    
    # open multi-file dataset (this function accepts unix wildcards)
    d = xr.open_mfdataset('/mnt/wrf_history/vol??/wrf_out/wy_' + str(year) + '/d02/wrfout_d02_*', drop_variables=var_list, concat_dim='Time')
    
    # Get mean/min/max by day of year for desired variables 
    new_array = d[['T2','Q2','SWDOWN','SWNORM']].groupby('XTIME.dayofyear').mean(dim='Time') # create daily means of few variables
    new_array['TMIN'] = d['T2'].groupby('XTIME.dayofyear').min(dim='Time') # create daily minimum temperature
    new_array['TMAX'] = d['T2'].groupby('XTIME.dayofyear').max(dim='Time') # create daily maximum temperature
    new_array = new_array.rename({'T2' : 'TMEAN'}) # rename T2 as TMEAN
    
    # Adjust some meta data
    new_array['TMEAN'].attrs = [('description','DAILY MEAN GRID SCALE TEMPERATUTE'), ('units','K')]
    new_array['TMIN'].attrs = [('description','DAILY MINIMUM GRID SCALE TEMPERATURE'), ('units','K')]
    new_array['TMAX'].attrs = [('description','DAILY MAXIMUM GRID SCALE TEMPERATURE'), ('units','K')]
    new_array['Q2'].attrs = [('description','DAILY MEAN GRID SCALE SPECIFIC HUMIDITY'), ('units','')]
    new_array['SWDOWN'].attrs = [('description','DAILY MEAN DOWNWARD SHORT WAVE FLUX AT GROUND SURFACE'), ('units','W m^2')]
    new_array['SWNORM'].attrs = [('description','DAILY MEAN NORMAL SHORT WAVE FLUX AT GROUND SURFACE (SLOPE-DEPENDENT)'), ('units','W m^2')]

    # Write new netcdf file
    new_array.to_netcdf("/mnt/selway/data/data_02/charlie/subsets/test/forLejo/Biome-BGC-WY-" + str(year) + ".nc")
	
    del d, new_array	

t1 = time.time()
print("Total time to create this subset was:", t1 - t0, "seconds.")

