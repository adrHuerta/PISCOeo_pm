rm(list = ls())

library(xts)
library(raster)
"%>%" = magrittr::`%>%`

source('./src/from_PISCOt/Merging/MG_normal_anomaly_values.R')

OBS_data <- readRDS("./data/processed/obs/td/OBS_td.RDS")

# normals values based on grid
# to approximate the obs value as much as possible

td_normals <- file.path("./data/processed/gridded/sub_variables/normals",
                        sprintf("%s/td_%02d.nc", "td", 1:12)) %>%
  lapply(function(x) raster::raster(x)) %>%
  raster::brick() %>%
  raster::extract(x = ., y = OBS_data$xyz) %>%
  t()

# to anomaly values
anomaly_td <- lapply(seq_len(ncol(OBS_data$values$td)),
                     function(x){
                       get_anomaly_values2(daily_time_serie = OBS_data$values$td[, x],
                                           monthly_values = td_normals[, x])
                     }) %>% 
  do.call("cbind", .)

# save data
saveRDS(object = list(values = list(td = anomaly_td),
                      xyz = OBS_data$xyz),
        file = "./data/processed/obs/td/Anomalies_OBS_td.RDS")
