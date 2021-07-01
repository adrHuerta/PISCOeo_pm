rm(list = ls())

source('./src/checking_files.R')

# checking files
checking_files(var_grid = "tmax", time_range = c("1981-01-01", "2016-12-31"))
checking_files(var_grid = "tmin", time_range = c("1981-01-01", "2016-12-31"))
checking_files(var_grid = "td", time_range = c("1981-01-01", "2016-12-31"))
checking_files(var_grid = "sd", time_range = c("1981-01-01", "2016-12-31"))

# checking no NA files
summary(get_pixel_ts_from_grids(var_grid = "tmax"))
summary(get_pixel_ts_from_grids(var_grid = "tmin"))
summary(get_pixel_ts_from_grids(var_grid = "td"))
summary(get_pixel_ts_from_grids(var_grid = "sd"))

# python code
reticulate::use_virtualenv("/home/adrian/Documents/Repos/Evapotranspiration/venv", required = TRUE)
reticulate::repl_python()

# packages
import urllib.request
import pandas as pd
import numpy as np
import datetime
import xarray as xr
from joblib import Parallel, delayed

# (FAO56) Penman-Monteith (PM) equation for FT
PM_equation = urllib.request.urlopen("https://git.io/J3Iy5")
exec(PM_equation.read())

# range time for computation
range_time = pd.date_range("1981-01-01", "2016-12-31", freq="d")

# path files
file_path_grids = "data/processed/gridded/"

# PM function
def apply_PM_function(time_step):
  
  Time_i = int(time_step.strftime('%j'))
  Tmax_i = xr.open_dataset(file_path_grids + "sub_variables/values/tmax/tmax_" + time_step.strftime('%Y-%m-%d') + ".nc").layer
  Tmin_i = xr.open_dataset(file_path_grids + "sub_variables/values/tmin/tmin_" + time_step.strftime('%Y-%m-%d') + ".nc").layer
  Sd_i = xr.open_dataset(file_path_grids + "sub_variables/values/sd/sd_" + time_step.strftime('%Y-%m-%d') + ".nc").layer
  Td_i = xr.open_dataset(file_path_grids + "sub_variables/values/td/td_" + time_step.strftime('%Y-%m-%d') + ".nc").layer
  Ws_m = xr.open_dataset(file_path_grids + "sub_variables/normals/ws/ws_" + time_step.strftime('%m') + ".nc").layer
  LAT = xr.open_dataset(file_path_grids + "co_variables/Y.nc").Y
  Z = xr.open_dataset(file_path_grids + "co_variables/DEM.nc").DEM
  
  pm_eo_evap = penman_monteith_FAO56_FT(time_i = Time_i,
                                        tmax_i = Tmax_i.values,
                                        tmin_i = Tmin_i.values,
                                        sd_i = Sd_i.values,
                                        td_i = Td_i.values,
                                        u2_i = Ws_m.values,
                                        lat_i = LAT.values,
                                        z_i = Z.values)
  
  to_save = xr.open_dataset(file_path_grids + "co_variables/DEM.nc")
  to_save["eo"] = (('latitude', 'longitude'), pm_eo_evap)
  to_save["eo"].to_netcdf(file_path_grids + "eo/eo_" + time_step.strftime('%Y-%m-%d') + ".nc")


# applying PM function
Parallel(n_jobs=2, verbose=50)(
  delayed(apply_PM_function)(i) for i in range_time
  )
