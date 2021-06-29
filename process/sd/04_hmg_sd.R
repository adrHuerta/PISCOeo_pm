library(lubridate)
library(xts)

"%>%" = magrittr::`%>%`

source("src/qc_spatial_neighbors.R")
source("src/qc_functions.R")
source("src/from_PISCOt/GapFilling/GF_build_neigh_matrix.R")
source("src/from_PISCOt/GapFilling/GF_daily_climatology_filling.R")
source('src/from_PISCOt/Homogenization/HG_pairwise_snht.R')
source('src/from_PISCOt/Homogenization/HG_simple_snht.R')
source('src/from_PISCOt/Homogenization/HG_hmgFactor2daily.R')

# E1: import data
sd_gf   <- readRDS("data/processed/obs/sd/qc_gf_sd_obs.RDS")
sd_gf_data <- sd_gf$values

# E2: daily to monthly
month_sd <- do.call("cbind", lapply(sd_gf_data, function(x) round(xts::apply.monthly(x, mean), 2)))
qc_monthly_values_rs_ERA5_hmg <- month_sd

# E3: Search neighbor stations to complete by parameters:
param_spt <- list(list(lmt_dist = 1000,
                       lmt_elv = 1000,
                       lmt_n = 8),
                  list(lmt_dist = 1000,
                       lmt_elv = 1000,
                       lmt_n = 8),
                  list(lmt_dist = 1000,
                       lmt_elv = 1000,
                       lmt_n = 8))

# E4: monthly homogenization 
for(xi in seq_along(param_spt)){
  
  sd_to_be_hmg <- list()
  
  for (station_j in seq_along(sd_gf$xyz$ID)){
    ID_stat_s  <- spt_neighrs(id_station = sd_gf$xyz$ID[station_j],
                              stations_database = sd_gf$xyz,
                              lmt_dist = param_spt[[xi]]$lmt_dist,
                              lmt_elv = param_spt[[xi]]$lmt_elv,
                              lmt_n = param_spt[[xi]]$lmt_n)
    build_matrix(id_stations = ID_stat_s,
                 time_series_database = qc_monthly_values_sd_ERA5_hmg,
                 neigh_r = .55) -> ID_stat_s
    
    if(ncol(ID_stat_s) >= 4){
      
      response <- pha_hmg(ts_data = ID_stat_s)$hmg
      
    } else {
      
      response <- snht_hmg(ts_data = ID_stat_s[, 1])$hmg
      
    }
    
    sd_to_be_hmg[[station_j]] <- response
  }
  
  qc_monthly_values_sd_ERA5_hmg <- do.call("cbind", sd_to_be_hmg)
  colnames(qc_monthly_values_sd_ERA5_hmg) <- sd_gf$xyz$ID
}

# E5: daily correction
qc_daily_values_sd_ERA5_hmg <- sd_gf_data

for(station_j in sd_gf$xyz$ID){
  qc_daily_values_sd_ERA5_hmg[, station_j] <- hmgFactor2daily(monthly_ts = month_sd[, station_j],
                                                                monthly_ts_hmg = qc_monthly_values_sd_ERA5_hmg[, station_j],
                                                                daily_ts = sd_gf_data[, station_j])
}

qc_daily_values_sd_ERA5_hmg[qc_daily_values_sd_ERA5_hmg<0] <- 0

sd_hmg <- list(values = xts(qc_daily_values_sd_ERA5_hmg, order.by = index(sd_gf$values)), 
               xyz=sd_gf$xyz)
saveRDS(sd_hmg, 'data/processed/obs/sd/qc_gf_hmg_sd_obs.RDS')
