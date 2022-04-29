# -*- coding: utf-8 -*-
"""ERA5land_2m_temperature.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1ZqNgYH0iE_ZIzafH2D-f-7mKkpspvWL4
"""

from google.colab import drive
drive.mount('/content/drive')

url = 'url: https://cds.climate.copernicus.eu/api/v2'
key = 'key: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

with open('/root/.cdsapirc', 'w') as f:
    f.write('\n'.join([url, key]))

with open('/root/.cdsapirc') as f:
    print(f.read())

!pip install cdsapi

import cdsapi
import glob
c = cdsapi.Client()

for year in range(1981, 2021):
    c.retrieve("reanalysis-era5-land",
               {
                   "variable": "2m_temperature",
                   "product_type": "reanalysis",
                   "year": str(year),
                   'month': ['01'],
                   'day': ['01', '02', '03',
                           '04', '05', '06',
                           '07', '08', '09',
                           '10', '11', '12',
                           '13', '14', '15',
                           '16', '17', '18',
                           '19', '20', '21',
                           '22', '23', '24',
                           '25', '26', '27',
                           '28', '29', '30',
                           '31'],
                   'time': ['00:00', '01:00', '02:00',
                            '03:00', '04:00', '05:00',
                            '06:00', '07:00', '08:00',
                            '09:00', '10:00', '11:00',
                            '12:00', '13:00', '14:00',
                            '15:00', '16:00', '17:00',
                            '18:00', '19:00', '20:00',
                            '21:00', '22:00', '23:00'],
                   "area": "2/-82/-19/-66",  # North, West, South, East
                   "format": "netcdf"
               },
               "/content/drive/MyDrive/Google_Colab_temp/ERA5land/2m_temperature/" + str(year) + ".nc")

import xarray as xr
import numpy as np
import pandas as pd
import glob
!pip install netcdf4
!pip install rasterio

for year_i in range(1981,2021):
  print(year_i)
  netcdf_year_i = xr.open_dataset("/content/drive/MyDrive/Google_Colab_temp/ERA5land/2m_temperature/" + str(year_i) + ".nc")
  netcdf_year_i_plus_1 = xr.open_dataset("/content/drive/MyDrive/Google_Colab_temp/ERA5land/2m_temperature/" + str(year_i+1) + ".nc")
  netcdf_year_i_plus_1 = netcdf_year_i_plus_1.sel(time=slice(str(year_i+1) + "-01-01T00:00:00", str(year_i+1) + "-01-01T04:00:00"))

  netcdf_year_i = xr.concat([netcdf_year_i, netcdf_year_i_plus_1], dim="time")
  netcdf_year_i["time"] = pd.date_range(str(year_i) + '-01-01', freq="H", periods=len(netcdf_year_i.time)) - pd.Timedelta(hours=5)
  netcdf_year_i = netcdf_year_i.where((-82 < netcdf_year_i.longitude) & (netcdf_year_i.longitude < -67) & (-19 < netcdf_year_i.latitude) & (netcdf_year_i.latitude < 2), drop=True)
  netcdf_year_i = netcdf_year_i-273.15

  # variables
  var1 = netcdf_year_i.resample(time="1D").max(dim="time").sel(time=slice(str(year_i) + "-01-01", str(year_i) + "-12-31")).t2m
  encoding_tmax = {v: {'zlib': True, 'complevel': 5} for v in ["tmax"]}
  (var1.to_dataset(name="tmax")).to_netcdf("/content/drive/MyDrive/Google_Colab_temp/ERA5land/for_eo_computation/tmax/tmax_" + str(year_i) + ".nc", encoding=encoding_tmax, engine='netcdf4')
  
  var2 = netcdf_year_i.resample(time="1D").min(dim="time").sel(time=slice(str(year_i) + "-01-01", str(year_i) + "-12-31")).t2m
  encoding_tmin = {v: {'zlib': True, 'complevel': 5} for v in ["tmin"]}
  (var2.to_dataset(name="tmin")).to_netcdf("/content/drive/MyDrive/Google_Colab_temp/ERA5land/for_eo_computation/tmin/tmin_" + str(year_i) + ".nc", encoding=encoding_tmin, engine='netcdf4')

ERA5_Land_tmax = sorted(glob.glob("/content/drive/MyDrive/Google_Colab_temp/ERA5land/for_eo_computation/tmax/*.nc"))
ERA5_Land_tmax = [xr.open_dataset(grid) for grid in ERA5_Land_tmax]
ERA5_Land_tmax = xr.concat(ERA5_Land_tmax, dim="time")
ERA5_Land_tmax.to_netcdf("/content/drive/MyDrive/repos/PISCOt_v1-2/data/gridded/co_variables/ERA5land_tmax_1981_2020.nc")
ERA5_Land_tmax

ERA5_Land_tmin = sorted(glob.glob("/content/drive/MyDrive/Google_Colab_temp/ERA5land/for_eo_computation/tmin/*.nc"))
ERA5_Land_tmin = [xr.open_dataset(grid) for grid in ERA5_Land_tmin]
ERA5_Land_tmin = xr.concat(ERA5_Land_tmin, dim="time")
ERA5_Land_tmin.to_netcdf("/content/drive/MyDrive/repos/PISCOt_v1-2/data/gridded/co_variables/ERA5land_tmin_1981_2020.nc")
ERA5_Land_tmin

