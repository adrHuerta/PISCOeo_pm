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

# read data and xyz from td obs & td ERA5
td_cws <- readRDS('data/processed/obs/td/qc_td_plus_era_data.RDS')
td_cws_data <- td_cws$values
td_cws_xyz <- td_cws$xyz

# Select stations to complete:
td_cws_xyz_to_be_used_plusERA5 <- td_cws_xyz # stations + ERA5
td_cws_xyz <- td_cws_xyz[td_cws_xyz$SRC == 'OBS' & td_cws_xyz$QC == 1,] # stations
qc_data_values_td_ERA5_filled <- td_cws_data

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
  
  td_era_to_be_fill <- parallel::mclapply(seq_along(td_cws_xyz$ID),
                                          function(station_j){
                                            
                                            spt_neighrs(id_station = td_cws_xyz$ID[station_j],
                                                         stations_database = td_cws_xyz_to_be_used_plusERA5,
                                                         lmt_dist = param_spt[[xi]]$lmt_dist,
                                                         lmt_elv = param_spt[[xi]]$lmt_elv,
                                                         lmt_n = param_spt[[xi]]$lmt_n) -> step1
                                            
                                            build_matrix(id_stations = step1,
                                                         time_series_database = qc_data_values_td_ERA5_filled,
                                                         param_neigh = list(r_cor = .6, n_daily_cycle = 5)) -> step2
                                            
                                            std_dep_imputation(stat_data = step2) -> step3
                                            step3$filled
                                             
                                           }, mc.cores = 10)
  
  do.call("cbind", td_era_to_be_fill) -> td_ERA5_to_be_filled
  setNames(td_ERA5_to_be_filled, td_cws_xyz$ID) -> td_ERA5_to_be_filled
  
  qc_data_values_td_ERA5_filled[, colnames(td_ERA5_to_be_filled)] <- td_ERA5_to_be_filled
}

qc_data_td_ERA5_filled <- qc_data_values_td_ERA5_filled[, td_cws_xyz$ID]

qc_data_td_ERA5_filled_xts <- xts(qc_data_td_ERA5_filled, order.by = index(td_cws_data))
qc_data_td_ERA5_filled_rds <- list(values=qc_data_td_ERA5_filled_xts, xyz = td_cws$xyz)

# gap-filling (2)
# gap-filling using daily climatology
# (only in stations where there are at least 1% (or less) of NAs)

for(station_j in colnames(qc_data_td_ERA5_filled)){
  
  sample_station_j <- qc_data_td_ERA5_filled[, station_j]
  size_percent_na <- sum(is.na(sample_station_j))*100/length(sample_station_j)
  
  if(size_percent_na <= 1 & size_percent_na > 0) {
    
    qc_data_td_ERA5_filled[, station_j] <- daily_climatology_filling(ts_data = sample_station_j)
    
  } else {
    
    next
  }
  
}

# selection without NA stations
td_ERA5_filled_sel <- insel_na_3(qc_data_td_ERA5_filled)
td_cws_xyz <- td_cws$xyz

td_xyz <- data.frame()
for (est in names(td_ERA5_filled_sel)) {
  xyz_sel <- td_cws_xyz[td_cws_xyz$ID==est,]
  td_xyz <- rbind(td_xyz, xyz_sel)
}

td_ERA5_filled_sel_xts <- xts(td_ERA5_filled_sel, order.by = index(td_cws$values))
row.names(td_xyz) <- NULL
td_cws_gf <- list(values = td_ERA5_filled_sel_xts, xyz=td_xyz)

saveRDS(td_cws_gf, 'data/processed/obs/td/qc_gf_td_obs.RDS')
