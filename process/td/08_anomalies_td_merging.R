rm(list = ls())

library(xts)
library(raster)
library(gstat)
"%>%" = magrittr::`%>%`

source('./src/from_PISCOt/Merging/MG_make_covariables.R')
source('./src/from_PISCOt/Merging/MG_RK.R')

output_anomalies <- "./data/processed/gridded/sub_variables/values"

# obs
qc_data <- readRDS("./data/processed/obs/td/Anomalies_OBS_td.RDS")

# gridded
LST_mean <- raster::brick("./data/processed/gridded/co_variables/LST_mean.nc")
DEM <- raster::raster("./data/processed/gridded/co_variables/DEM.nc")/1000
X <- raster::raster("./data/processed/gridded/co_variables/X.nc")
Y <- raster::raster("./data/processed/gridded/co_variables/Y.nc")
tdi_grided <- raster::raster("./data/processed/gridded/co_variables/TDI.nc")

# gridded
td_normals <- file.path("./data/processed/gridded/sub_variables/normals",
                        sprintf("%s/td_%02d.nc", "td", 1:12)) %>%
  lapply(function(x) raster::raster(x)) %>%
  raster::brick()

# making list of covs
covs_list_td <- list(dynamic = list(LST = LST_mean),
                     static = list(DEM = DEM, X = X, Y = Y, TDI = tdi_grided))


parallel::mclapply(seq_along(time(qc_data$values$td)), 
                   function(i){
                     
                     date_i <- time(qc_data$values$td)[i]
                     month_value_i <- as.numeric(format(as.Date(date_i), "%m"))
                     
                     
                     td_i <- make_Anomaly_coVariables(day_date = date_i,
                                                      var = "td",
                                                      covs_list = covs_list_td,
                                                      obs = qc_data)
                     
                     grid_i  <- (RK(obs_cov_data = td_i, resFitting = 10) + td_normals[[month_value_i]])
                     grid_i <- round(grid_i, 2)
                     
                     raster::writeRaster(x = grid_i, 
                                         filename = file.path(output_anomalies, 
                                                              sprintf("%s/td_%s.nc", "td",  date_i)),
                                         datatype = 'FLT4S', force_v4 = TRUE, compression = 7)
                     
                   }, mc.cores = 10)
