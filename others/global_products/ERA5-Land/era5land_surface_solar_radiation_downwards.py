# -*- coding: utf-8 -*-
"""ERA5land_surface_solar_radiation_downwards.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1enmP03gEm6HOAYHkXh_WlfVmpMEWLHKO
"""

from google.colab import drive
drive.mount('/content/drive')

# url = 'url: https://cds.climate.copernicus.eu/api/v2'
# key = 'key: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

# with open('/root/.cdsapirc', 'w') as f:
#     f.write('\n'.join([url, key]))

# with open('/root/.cdsapirc') as f:
#     print(f.read())

# !pip install cdsapi

# import cdsapi
# import glob
# c = cdsapi.Client()

# for year in range(1981, 2021):
#     c.retrieve("reanalysis-era5-land",
#                {
#                    "variable": "surface_solar_radiation_downwards",
#                    "product_type": "reanalysis",
#                    "year": str(year),
#                    'month': ['01', '02', '03',
#                              '04', '05', '06',
#                              '07', '08', '09',
#                              '10', '11', '12'],
#                    'day': ['01', '02', '03',
#                            '04', '05', '06',
#                            '07', '08', '09',
#                            '10', '11', '12',
#                            '13', '14', '15',
#                            '16', '17', '18',
#                            '19', '20', '21',
#                            '22', '23', '24',
#                            '25', '26', '27',
#                            '28', '29', '30',
#                            '31'],
#                    'time': ['00:00', '01:00', '02:00',
#                             '03:00', '04:00', '05:00',
#                             '06:00', '07:00', '08:00',
#                             '09:00', '10:00', '11:00',
#                             '12:00', '13:00', '14:00',
#                             '15:00', '16:00', '17:00',
#                             '18:00', '19:00', '20:00',
#                             '21:00', '22:00', '23:00'],
#                    "area": "2/-82/-19/-66",  # North, West, South, East
#                    "format": "netcdf"
#                },
#                "/content/drive/MyDrive/Google_Colab_temp/ERA5land/surface_solar_radiation_downwards/" + str(year) + ".nc")

import xarray as xr
import numpy as np
import pandas as pd
import glob
!pip install netcdf4
!pip install rasterio

# https://confluence.ecmwf.int/display/CUSF/ERA5+Land+Hourly+data+-+Surface++Net+Solar+Radiation
# https://confluence.ecmwf.int/pages/viewpage.action?pageId=197702790
# https://edis.ifas.ufl.edu/pdf%5CAE%5CAE45900.pdf

data_pseudo_yearly = [] 

for year_i in range(1981,2021):
  print(year_i)
  netcdf_year_i = xr.open_dataset("/content/drive/MyDrive/Google_Colab_temp/ERA5land/surface_solar_radiation_downwards/" + str(year_i) + ".nc")

  # following https://confluence.ecmwf.int/display/CUSF/ERA5+Land+Hourly+data+-+Surface++Net+Solar+Radiation
  # as ssrd is an accumulated variable, I just compute the daily from the UTC. No UTC to local time conversion is done
  
  netcdf_year_i = netcdf_year_i.sel(time=netcdf_year_i.time.dt.hour == 00)
  netcdf_year_i = (netcdf_year_i/86400)*0.0864 # to MJ/m??day??
  time_step_corrected = pd.date_range(pd.to_datetime(str(netcdf_year_i.isel(time=0).time.values)).strftime('%Y-%m-%d'), freq="D", periods=len(netcdf_year_i.time)) - pd.Timedelta(days=1)
  netcdf_year_i["time"] = time_step_corrected
  netcdf_year_i = netcdf_year_i.where((-82 < netcdf_year_i.longitude) & (netcdf_year_i.longitude < -67) & (-19 < netcdf_year_i.latitude) & (netcdf_year_i.latitude < 2), drop=True)
  
  # variables
  data_pseudo_yearly.append(netcdf_year_i.rename({"ssrd":"rs"}))

data_pseudo = xr.concat(data_pseudo_yearly, dim="time")

for year_i in range(1981,2020):
  print(year_i)
  var_year = data_pseudo.sel(time=data_pseudo.time.dt.year == year_i)
  var_year = var_year.astype("float32")
  encoding = {v: {'zlib': True, 'complevel': 5} for v in ["rs"]}
  var_year.to_netcdf("/content/drive/MyDrive/Google_Colab_temp/ERA5land/for_eo_computation/rs/rs_" + str(year_i) + ".nc", encoding=encoding, engine='netcdf4')

