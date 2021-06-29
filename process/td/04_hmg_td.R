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
td_gf   <- readRDS("data/processed/obs/td/qc_gf_td_obs.RDS")
td_gf_data <- td_gf$values

# E2: daily to monthly
month_td <- do.call("cbind", lapply(td_gf_data, function(x) round(xts::apply.monthly(x, mean), 2)))
qc_monthly_values_td_ERA5_hmg <- month_td

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
  
  td_to_be_hmg <- list()
  
  for (station_j in seq_along(td_gf$xyz$ID)){
    ID_stat_s  <- spt_neighrs(id_station = td_gf$xyz$ID[station_j],
                              stations_database = td_gf$xyz,
                              lmt_dist = param_spt[[xi]]$lmt_dist,
                              lmt_elv = param_spt[[xi]]$lmt_elv,
                              lmt_n = param_spt[[xi]]$lmt_n)
    build_matrix(id_stations = ID_stat_s,
                 time_series_database = qc_monthly_values_td_ERA5_hmg,
                 neigh_r = .6) -> ID_stat_s
    
    if(ncol(ID_stat_s) >= 4){
      
      response <- pha_hmg(ts_data = ID_stat_s)$hmg
      
    } else {
      
      response <- snht_hmg(ts_data = ID_stat_s[, 1])$hmg
      
    }
    
    td_to_be_hmg[[station_j]] <- response
  }
  
  qc_monthly_values_td_ERA5_hmg <- do.call("cbind", td_to_be_hmg)
  colnames(qc_monthly_values_td_ERA5_hmg) <- td_gf$xyz$ID
}

# E5: daily correction
qc_daily_values_td_ERA5_hmg <- td_gf_data

for(station_j in td_gf$xyz$ID){
  qc_daily_values_td_ERA5_hmg[, station_j] <- hmgFactor2daily(monthly_ts = month_td[, station_j],
                                                              monthly_ts_hmg = qc_monthly_values_td_ERA5_hmg[, station_j],
                                                              daily_ts = td_gf_data[, station_j])
}

qc_daily_values_td_ERA5_hmg[qc_daily_values_td_ERA5_hmg<0] <- 0

td_hmg <- list(values = xts(qc_daily_values_td_ERA5_hmg, order.by = index(td_gf$values)), 
               xyz=td_gf$xyz)
saveRDS(td_hmg, 'data/processed/obs/td/qc_gf_hmg_td_obs.RDS')
