rm(list = ls())

library(xts)
library(raster)
"%>%" = magrittr::`%>%`

source('./src/from_PISCOt/Merging/MG_make_single_point.R')
source('./src/from_PISCOt/Merging/MG_normal_anomaly_values.R')

# data
qc_data <- readRDS("./data/processed/obs/td/qc_gf_hmg_td_obs.RDS")

# grid data
gridded_data <- raster::brick("./data/processed/gridded/co_variables/DEM.nc")[[1]]
gridded_data_a <- raster::aggregate(gridded_data, 5)

# data to sp
xyz_sp <- sp::SpatialPointsDataFrame(coords = qc_data$xyz[, c("LON", "LAT")],
                                     data = qc_data$xyz,
                                     proj4string = sp::CRS(raster::projection(gridded_data)))

# point that are in the same
points_to_be_merged <- make_single_point(pts = xyz_sp,
                                         rgrid = gridded_data_a)

# getting normal values
normal_td <- lapply(qc_data$values, function(x) get_monthly_normals(daily_time_serie = x)) %>%
  do.call("cbind", .)

# save data
saveRDS(object = list(values = list(td = normal_td),
                      xyz = xyz_sp),
        file = "./data/processed/obs/td/Normals_OBS_td.RDS")

saveRDS(object = list(values = list(td = qc_data$values),
                      xyz = xyz_sp),
        file = "./data/processed/obs/td/OBS_td.RDS")