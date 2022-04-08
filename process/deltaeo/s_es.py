# -*- coding: utf-8 -*-
"""S_es.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1kbMoZ9uCT7hXYa3LdXjdffl8FxlZR_AY
"""

from google.colab import drive
drive.mount('/content/drive/')

# Commented out IPython magic to ensure Python compatibility.
# %%capture
# !pip install netcdf4
# !pip install geopandas

import xarray as xr
import numpy as np

"""
# Equation
def D_es(DELTA_i, GAMMA_i, ws_i, tmax_i, tmin_i):
  d_es_d_eo = (900*GAMMA_i*ws_i)/((DELTA_i + GAMMA_i*(0.34*ws_i + 1))*((tmax_i/2) + (tmin_i/2) + 273))
  # response = d_tmax_d_eo * ((tmax_i - tmax_min)/eo_i)
  return np.round(d_es_d_eo, 2)
"""

def S_es(DELTA_i, GAMMA_i, ws_i, tmax_i, tmin_i, es_i, es_min, eo_i):
  d_es_d_eo = (900*GAMMA_i*ws_i)/((DELTA_i + GAMMA_i*(0.34*ws_i + 1))*((tmax_i/2) + (tmin_i/2) + 273))
  response = d_es_d_eo * ((es_i - es_min)/eo_i)
  return np.round(response, 2)

path_netcd_01 = "/content/drive/MyDrive/DATA_FT/version_rc/"

for year_i in range(1989,1990):
  for month_i in range(1, 13):

    delta_daily = xr.open_dataset(path_netcd_01+"delta/daily/delta_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc")
    gamma_daily = xr.open_dataset(path_netcd_01+"gamma/gamma.nc")
    ws_daily = xr.open_dataset(path_netcd_01+"/ws/normal/ws_mean.nc").isel(time=month_i-1)
    #ea_daily = xr.open_dataset(path_netcd_01+"ea/daily/ea_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc")
    #es_daily = xr.open_dataset(path_netcd_01+"es/daily/es_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc")
    tmax_daily = xr.open_dataset(path_netcd_01+"tmax/daily/tmax_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc")
    #tmax_min = xr.open_dataset("/content/drive/MyDrive/variables_min_clim/var_min_clim.nc")
    tmin_daily = xr.open_dataset(path_netcd_01+"tmin/daily/tmin_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc")
    es_daily = xr.open_dataset(path_netcd_01+"es/daily/es_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc")     
    es_min = xr.open_dataset("/content/drive/MyDrive/variables_min_clim/var_min_clim.nc")
    eo_daily = xr.open_dataset(path_netcd_01+"eo/daily/eo_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc")

    Ses = xr.apply_ufunc(S_es,
                           delta_daily.delta,
                           gamma_daily.gamma,
                           ws_daily.ws,
                           tmax_daily.tmax,
                           tmin_daily.tmin,
                         es_daily.es,
                         es_min.es,
                         eo_daily.eo,
                           input_core_dims=[['time'], [], [], ["time"], ["time"], ["time"], [], ["time"]],
                           output_core_dims=[["time"]],
                           vectorize=True,
                           output_dtypes=['float32'])
    
    encoding = {v: {'zlib': True, 'complevel': 9} for v in ["Ses"]}
    Ses.to_dataset(name="Ses").to_netcdf(path_netcd_01+"S/Ses/Ses_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc", encoding=encoding, engine='netcdf4')

Ses

"""
for year_i in range(1981,1982):
  for month_i in range(1, 2):

    delta_daily = xr.open_dataset(path_netcd_01+"delta/daily/delta_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc")
    gamma_daily = xr.open_dataset(path_netcd_01+"gamma/gamma.nc")
    ws_daily = xr.open_dataset(path_netcd_01+"/ws/normal/ws_mean.nc").isel(time=month_i-1)
    #ea_daily = xr.open_dataset(path_netcd_01+"ea/daily/ea_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc")
    #es_daily = xr.open_dataset(path_netcd_01+"es/daily/es_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc")
    tmax_daily = xr.open_dataset(path_netcd_01+"tmax/daily/tmax_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc")
    #tmax_min = xr.open_dataset("/content/drive/MyDrive/variables_min_clim/var_min_clim.nc")
    tmin_daily = xr.open_dataset(path_netcd_01+"tmin/daily/tmin_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc")    
    #eo_daily = xr.open_dataset(path_netcd_01+"eo/daily/eo_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc")

    Des = xr.apply_ufunc(D_es,
                           delta_daily.delta,
                           gamma_daily.gamma,
                           ws_daily.ws,
                           tmax_daily.tmax,
                           tmin_daily.tmin,
                           input_core_dims=[['time'], [], [], ["time"], ["time"]],
                           output_core_dims=[["time"]],
                           vectorize=True,
                           output_dtypes=['float32'])
    
    encoding = {v: {'zlib': True, 'complevel': 9} for v in ["Des"]}
    Des.to_dataset(name="Des").to_netcdf(path_netcd_01+"D/Des/Des_daily"+"_"+str(year_i)+ "_"+str(month_i).zfill(2)+".nc", encoding=encoding, engine='netcdf4')
"""

