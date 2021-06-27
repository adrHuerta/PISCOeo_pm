rm(list = ls())

library(raster)
library(gstat)
"%>%" = magrittr::`%>%`

source('./src/from_PISCOt/Merging/MG_make_covariables.R')
source('./src/from_PISCOt/Merging/MG_RK.R')

output_normals <- "./data/processed/gridded/sub_variables/normals"

# obs
qc_data <- readRDS("./data/processed/obs/ws/Normals_OBS_ws.RDS")
qc_data$values$ws <- sqrt(qc_data$values$ws)

# gridded
WS_worldclim <- raster::brick("./data/processed/gridded/co_variables/WS_worldclim.nc")

# making list of covs
covs_list_ws <- list(dynamic = list(WS_worldclim = WS_worldclim),
                     static = NA)

for(i in 1:12){
  
  cc_i <- make_Normal_coVariables(month_value = i,
                                  var = "ws",
                                  covs_list = covs_list_ws,
                                  obs = qc_data)
  
  grid_i <- RK(obs_cov_data = cc_i, resFitting = 10)
  grid_i <- round(grid_i^2, 2)
  
  raster::writeRaster(x = grid_i, 
                      filename = file.path(output_normals, 
                                           sprintf("%s/ws_%02d.nc", "ws",  i)),
                      datatype = 'FLT4S', force_v4 = TRUE, compression = 7)
  
}
