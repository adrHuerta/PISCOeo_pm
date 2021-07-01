rm(list = ls())

library(xts)
library(raster)
library(gstat)
"%>%" = magrittr::`%>%`

source('./src/from_PISCOt/Merging/MG_make_covariables.R')
source('./src/from_PISCOt/Merging/MG_RK.R')
source('./src/physical_conditions_for_sub_variables.R')

output_anomalies <- "./data/processed/gridded/sub_variables/values"

# obs
qc_data <- readRDS("./data/processed/obs/sd/Anomalies_OBS_sd.RDS")

# gridded
CC <- raster::brick("./data/processed/gridded/co_variables/CC.nc")
DEM <- raster::raster("./data/processed/gridded/co_variables/DEM.nc")/1000
X <- raster::raster("./data/processed/gridded/co_variables/X.nc")
Y <- raster::raster("./data/processed/gridded/co_variables/Y.nc")
tdi_grided <- raster::raster("./data/processed/gridded/co_variables/TDI.nc")

# making list of covs
covs_list_sd <- list(dynamic = list(CC = CC),
                     static = list(DEM = DEM, X = X, Y = Y, TDI = tdi_grided))

# gridded
sd_normals <- file.path("./data/processed/gridded/sub_variables/normals",
                        sprintf("%s/sd_%02d.nc", "sd", 1:12)) %>%
  lapply(function(x) raster::raster(x)) %>%
  raster::brick()

# for Nmax
mask_for_Nmax <- CC[[1]]
mask_for_Nmax[mask_for_Nmax > 0] <- 1


parallel::mclapply(4001:13149, 
                   function(i){
                     
                     date_i <- time(qc_data$values$sd)[i]
                     month_value_i <- as.numeric(format(as.Date(date_i), "%m"))
                     
                     sd_i <- make_Anomaly_coVariables(day_date = date_i,
                                                      var = "sd",
                                                      covs_list = covs_list_sd,
                                                      obs = qc_data)
                     
                     grid_i  <- (RK(obs_cov_data = sd_i, resFitting = 10) + sd_normals[[month_value_i]])

                     # physical condition
                     max_sd <- maximum_lenght_sd(jday_i = as.numeric(format(date_i, "%j")), lat_i = Y) * mask_for_Nmax
                     grid_i[grid_i < 0] <- 0
                     grid_i@data@values[grid_i@data@values > max_sd@data@values] <- NA
                     grid_i <- raster::overlay(grid_i, max_sd, fun=function(x, y) ifelse(is.na(x), y, x))
                     
                     grid_i <- round(grid_i, 2)
                     
                     raster::writeRaster(x = grid_i, 
                                         filename = file.path(output_anomalies, 
                                                              sprintf("%s/sd_%s.nc", "sd",  date_i)),
                                         datatype = 'FLT4S', force_v4 = TRUE, compression = 7)
                     
                   }, mc.cores = 11)