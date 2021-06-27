rm(list = ls())

library(xts)
library(raster)
"%>%" = magrittr::`%>%`

source('./src/from_PISCOt/Merging/MG_normal_anomaly_values.R')

OBS_data <- readRDS("./data/processed/obs/sd/OBS_sd.RDS")

# normals values based on grid
# to approximate the obs value as much as possible

sd_normals <- file.path("./data/processed/gridded/sub_variables/normals",
                        sprintf("%s/sd_%02d.nc", "sd", 1:12)) %>%
  lapply(function(x) raster::raster(x)) %>%
  raster::brick() %>%
  raster::extract(x = ., y = OBS_data$xyz) %>%
  t()

# to anomaly values
anomaly_sd <- lapply(seq_len(ncol(OBS_data$values$sd)),
                     function(x){
                       get_anomaly_values2(daily_time_serie = OBS_data$values$sd[, x],
                                           monthly_values = sd_normals[, x])
                     }) %>% 
  do.call("cbind", .)

# save data
saveRDS(object = list(values = list(sd = anomaly_sd),
                      xyz = OBS_data$xyz),
        file = "./data/processed/obs/sd/Anomalies_OBS_sd.RDS")
