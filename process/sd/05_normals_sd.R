rm(list = ls())

library(xts)
library(raster)
"%>%" = magrittr::`%>%`

source('./src/from_PISCOt/Merging/MG_make_single_point.R')
source('./src/from_PISCOt/Merging/MG_normal_anomaly_values.R')

# data
qc_data <- readRDS("./data/processed/obs/sd/qc_gf_hmg_sd_obs.RDS")

# for cv
# xyz_cv <- read.csv("./data/processed/obs/PISCOeo_pm_xyz_for_cv.csv")
# qc_data$values <- qc_data$values[, match(as.character(xyz_cv$ID), colnames(qc_data$values))]
# qc_data$xyz <- qc_data$xyz[match(as.character(xyz_cv$ID), as.character(qc_data$xyz$ID)), ]

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
normal_sd <- lapply(qc_data$values, function(x) get_monthly_normals(daily_time_serie = x)) %>%
  do.call("cbind", .)


# save data
saveRDS(object = list(values = list(sd = normal_sd),
                      xyz = xyz_sp),
        file = "./data/processed/obs/sd/Normals_OBS_sd.RDS")

saveRDS(object = list(values = list(sd = qc_data$values),
                      xyz = xyz_sp),
        file = "./data/processed/obs/sd/OBS_sd.RDS")
