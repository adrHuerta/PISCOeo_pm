rm(list = ls())

library(raster)
library(gstat)
"%>%" = magrittr::`%>%`

source('./src/from_PISCOt/Merging/MG_make_covariables.R')
source('./src/from_PISCOt/Merging/MG_RK.R')
source('./src/physical_conditions_for_sub_variables.R')

output_normals <- "./data/processed/gridded/sub_variables/normals"

# obs
qc_data <- readRDS("./data/processed/obs/sd/Normals_OBS_sd.RDS")

# gridded
CC <- raster::brick("./data/processed/gridded/co_variables/CC.nc")
DEM <- raster::raster("./data/processed/gridded/co_variables/DEM.nc")/1000
X <- raster::raster("./data/processed/gridded/co_variables/X.nc")
Y <- raster::raster("./data/processed/gridded/co_variables/Y.nc")

# making list of covs

covs_list_sd <- list(dynamic = list(CC = CC),
                     static = list(DEM = DEM, X = X, Y = Y))

# for Nmax
mask_for_Nmax <- CC[[1]]
mask_for_Nmax[mask_for_Nmax > 0] <- 1

seq(as.Date("1981-01-01"), as.Date("1981-12-31"), by = "day") %>%
  .[format(., "%d") == "15"] %>% format("%j") %>% as.numeric() -> j_day_clim

for(i in 1:12){
  
  cc_i <- make_Normal_coVariables(month_value = i,
                                  var = "sd",
                                  covs_list = covs_list_sd,
                                  obs = qc_data)
  
  grid_i <- RK(obs_cov_data = cc_i, resFitting = 10)
  
  # physical condition
  max_sd <- maximum_lenght_sd(jday_i = j_day_clim[i], lat_i = Y) * mask_for_Nmax
  grid_i[grid_i < 0] <- 0
  grid_i@data@values[grid_i@data@values > max_sd@data@values] <- NA
  grid_i <- raster::overlay(grid_i, max_sd, fun=function(x, y) ifelse(is.na(x), y, x))
  
  grid_i <- round(grid_i, 2)

  raster::writeRaster(x = grid_i,
                      filename = file.path(output_normals,
                                           sprintf("%s/sd_%02d.nc", "sd",  i)),
                      datatype = 'FLT4S', force_v4 = TRUE, compression = 7)
  
  }

