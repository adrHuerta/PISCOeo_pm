rm(list = ls())

library(lubridate)
library(xts)
library(magrittr)

"%>%" = magrittr::`%>%`

source('src/reanalysis_bias.R')
source("src/qc_functions.R")
source("src/from_PISCOt/GapFilling/GF_build_neigh_matrix.R")
source("src/from_PISCOt/GapFilling/GF_std_dep_imputation.R")
source("src/from_PISCOt/GapFilling/GF_daily_climatology_filling.R")
source('src/qc_spatial_neighbors.R')

# read data and xyz from sd obs & sd ERA5
sd_cws <- readRDS('data/processed/obs/sd/qc_sd_plus_era_data.RDS')
sd_cws_data <- sd_cws$values
sd_cws_xyz <- sd_cws$xyz

# Select stations to complete:
sd_cws_xyz_to_be_used_plusERA5 <- sd_cws_xyz # stations + ERA5
sd_cws_xyz <- sd_cws_xyz[sd_cws_xyz$SRC == 'OBS' & sd_cws_xyz$QC == 1,] # stations
qc_data_values_sd_ERA5_filled <- sd_cws_data

# Search neighbor stations to complete by parameters:
# i) distance 70km, difference de 1000m, limit 8 stations
# ii)  distance 100km, difference de 1000m, limit 8 stations
# iii) distance 150km, difference de 5500m, limit 8 stations
param_spt <- list(list(lmt_dist = 70,
                       lmt_elv = 1000,
                       lmt_n = 8),
                  list(lmt_dist = 100,
                       lmt_elv = 1000,
                       lmt_n = 8),
                  list(lmt_dist = 150,
                       lmt_elv = 5500,
                       lmt_n = 8))

# gap-filling (2) 
# gap-filling using daily climatology
# (only in stations where there are at least 1% (or less) of NAs)
for(xi in seq_along(param_spt)){
  
  sd_era_to_be_fill <- parallel::mclapply(seq_along(sd_cws_xyz$ID),
                                          function(station_j){
                                            
                                            spt_neighrs(id_station = sd_cws_xyz$ID[station_j],
                                                        stations_database = sd_cws_xyz_to_be_used_plusERA5,
                                                        lmt_dist = param_spt[[xi]]$lmt_dist,
                                                        lmt_elv = param_spt[[xi]]$lmt_elv,
                                                        lmt_n = param_spt[[xi]]$lmt_n) -> step1
                                            
                                            build_matrix(id_stations = step1,
                                                         time_series_database = qc_data_values_sd_ERA5_filled,
                                                         list(r_cor = .6, n_daily_cycle = 5)) -> step2
                                            
                                            std_dep_imputation(stat_data = step2) -> step3
                                            step3$filled
                                            
                                          }, mc.cores = 4)

  do.call("cbind", sd_era_to_be_fill) -> sd_ERA5_to_be_filled
  setNames(sd_ERA5_to_be_filled, sd_cws_xyz$ID) -> sd_ERA5_to_be_filled
  
  qc_data_values_sd_ERA5_filled[, colnames(sd_ERA5_to_be_filled)] <- sd_ERA5_to_be_filled
}

qc_data_sd_ERA5_filled <- qc_data_values_sd_ERA5_filled[, sd_cws_xyz$ID]
qc_data_sd_ERA5_filled[qc_data_sd_ERA5_filled<0] <- 0

qc_data_sd_ERA5_filled_xts <- xts(qc_data_sd_ERA5_filled, order.by = index(sd_cws_data))
qc_data_sd_ERA5_filled_rds <- list(values=qc_data_sd_ERA5_filled_xts, xyz = sd_cws$xyz)

# gap-filling (2)
# gap-filling using daily climatology
# (only in stations where there are at least 1% (or less) of NAs)

for(station_j in colnames(qc_data_sd_ERA5_filled)){
  
  sample_station_j <- qc_data_sd_ERA5_filled[, station_j]
  size_percent_na <- sum(is.na(sample_station_j))*100/length(sample_station_j)
  
  if(size_percent_na <= 1 & size_percent_na > 0) {
    
    qc_data_sd_ERA5_filled[, station_j] <- daily_climatology_filling(ts_data = sample_station_j)
    
  } else {
    
    next
  }
  
}

# selection without NA stations
sd_ERA5_filled_sel <- insel_na_3(qc_data_sd_ERA5_filled)
sd_cws_xyz <- sd_cws$xyz

sd_xyz <- data.frame()
for (est in names(sd_ERA5_filled_sel)) {
  xyz_sel <- sd_cws_xyz[sd_cws_xyz$ID==est,]
  sd_xyz <- rbind(sd_xyz, xyz_sel)
}

sd_ERA5_filled_sel_xts <- xts(sd_ERA5_filled_sel, order.by = index(sd_cws$values))
row.names(sd_xyz) <- NULL
sd_cws_gf <- list(values = sd_ERA5_filled_sel_xts, xyz=sd_xyz)

saveRDS(sd_cws_gf, 'data/processed/obs/sd/qc_gf_sd_obs.RDS')
