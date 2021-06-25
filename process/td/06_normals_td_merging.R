rm(list = ls())

library(raster)
library(gstat)
"%>%" = magrittr::`%>%`

source('./src/from_PISCOt/Merging/MG_make_covariables.R')
source('./src/from_PISCOt/Merging/MG_RK.R')

output_normals <- "./data/processed/gridded/sub_variables/normals/td"

# obs
qc_data <- readRDS("./data/processed/obs/td/Normals_OBS_td.RDS")

# gridded
LST_mean <- raster::brick("./data/processed/gridded/co_variables/LST_mean.nc")
DEM <- raster::raster("./data/processed/gridded/co_variables/DEM.nc")
X <- raster::raster("./data/processed/gridded/co_variables/X.nc")
Y <- raster::raster("./data/processed/gridded/co_variables/Y.nc")

# making list of covs
covs_list_td <- list(dynamic = list(LST_mean = LST_mean),
                     static = list(DEM = DEM, X = X, Y = Y))

for(i in 1:12){
  
  cc_i <- make_Normal_coVariables(month_value = i,
                                  var = "td",
                                  covs_list = covs_list_td,
                                  obs = qc_data)
  
  grid_i <- RK(obs_cov_data = cc_i, resFitting = 10)

  
  raster::writeRaster(x = grid_i, 
                      filename = file.path(output_normals, 
                                           sprintf("Normals_%s/td_%02d.nc", "td",  i)),
                      datatype = 'FLT4S', force_v4 = TRUE, compression = 7)
  
  }
