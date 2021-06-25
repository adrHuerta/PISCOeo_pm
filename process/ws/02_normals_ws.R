rm(list = ls())

library(xts)
library(raster)
"%>%" = magrittr::`%>%`

source('./src/from_PISCOt/Merging/MG_make_single_point.R')
source('./src/from_PISCOt/Merging/MG_normal_anomaly_values.R')

# data
qc_data <- readRDS("./data/processed/obs/ws/qc_data_ws_obs.RDS")

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
# filtering data that has at least 5-years of data

n_months_by_stations <- qc_data$values %>% 
  lapply(function(xts_obj){
    apply.monthly(xts_obj, function(x) ifelse(sum(!is.na(x)) > 20, mean(x, na.rm = TRUE), NA) ) %>%
      apply.yearly(function(x) ifelse(sum(!is.na(x)) >= 12, mean(x, na.rm = TRUE), NA) )
  }) %>% do.call(cbind, .) %>%
  apply(2, function(x) sum(!is.na(x))) %>%
  .[. >= 5]

normal_ws <- lapply(qc_data$values[, names(n_months_by_stations)], function(x) get_monthly_normals(daily_time_serie = x)) %>%
  do.call("cbind", .)

# save data
xyz_sp <- xyz_sp[match(names(n_months_by_stations), xyz_sp@data$ID), ]
rownames(xyz_sp@data) <- NULL

saveRDS(object = list(values = list(ws = normal_ws),
                      xyz = xyz_sp),
        file = "./data/processed/obs/ws/Normals_OBS_ws.RDS")

saveRDS(object = list(values = list(ws = qc_data$values),
                      xyz = xyz_sp),
        file = "./data/processed/obs/ws/OBS_ws.RDS")
