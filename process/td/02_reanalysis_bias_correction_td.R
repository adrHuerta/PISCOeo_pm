rm(list = ls())

library(raster)
library(ncdf4)
library(lubridate)
library(xts)
library(lattice)
library(qmap)

source('src/qc_functions.R')
source('src/from_PISCOt/GapFilling/GF_dreqm.R')
source('src/from_PISCOt/GapFilling/GF_daily_climatology_filling.R')
source('src/reanalysis_bias.R')

# E1: import data reanalysis NetCDF
td_era5 <- raster::brick('data/processed/gridded/co_variables/ERA5land_td_1981_2019.nc') + 0
td_list <- readRDS('data/processed/obs/td/qc_td_obs.RDS')
td_xyz  <- td_list$xyz
td_xyz  <- td_xyz[td_xyz$QC==1,]
row.names(td_xyz) <- NULL

points_xy <- raster::extract(td_era5[[1]], td_xyz[,2:3], method = 'simple', cellnumbers = TRUE)[,1]
td_era <- as.data.frame(t(td_era5[points_xy]))
names(td_era) <- td_xyz$ID

td_era_xts <- xts(td_era, order.by = index(td_list$values))
td_era <- list(values=td_era_xts, xyz=td_xyz)

# E2: apply correction
td_era_xts <- td_era$values
td_era_id <- td_era$xyz$ID
td_list <- readRDS('data/processed/obs/td/qc_td_obs.RDS')
td_list_df <-  as.data.frame(td_list$values)
td_list_df <- td_list_df[,td_era_id]
td_data_xts <- xts(td_list_df, order.by = index(td_list$values))

td_erac <- parallel::mclapply(1:ncol(td_data_xts),
                              function(i){
                                
                                daily_varying_anom_qmap(ts_obs = td_data_xts[,i], 
                                                        ts_model = td_era_xts[,i])
                                
                              }, mc.cores = 5)

td_erac <- do.call("cbind", td_erac)
names(td_erac) <- names(td_data_xts)

td_erac_xts <- xts(td_erac, order.by = index(td_data_xts))
td_erac_list <- list(values=td_erac_xts, xyz=td_era$xyz)

# E3: ERA5 series extreme data control
td_erac <- as.data.frame(td_erac_list$values)
td_eracq <- qc3_func_1(td_erac, 3.5)
td_eracq[td_eracq>0] <- NA
td_eracq[td_eracq==0] <- 1
td_erac_c <- td_eracq[,1:ncol(td_erac)]*td_erac

td_eracc_xts <- xts(td_erac_c, order.by = as.Date(row.names(td_erac_c)))
td_eracc_list <- list(values=td_eracc_xts, xyz=td_erac_list$xyz)

# E4: keep if r >= 0.6 (ERA5 vs obs)
td_list <- readRDS('data/processed/obs/td/qc_td_obs.RDS')

td_era_xts <- td_eracc_list$values
td_era_id <- td_eracc_list$xyz$ID
td_data <-  as.data.frame(td_list$values)
td_data <- td_data[,td_era_id]
td_data_xts <- xts(td_data, order.by = index(td_list$values))

td_erarc <- compar_cc(as.data.frame(td_eracc_list$values),
                      as.data.frame(td_data_xts), 0.55)
td_erarc_xts <- xts(td_erarc, order.by = as.Date(row.names(td_erarc)))
td_erarc_list <- list(values=td_erarc_xts, xyz=td_eracc_list$xyz)

# E5: complete with daily climatologies
td_erarc_xts <- td_erarc_list$values

td_eracc <- as.data.frame((matrix(1:nrow(td_erarc_xts), ncol = 1)))

for (i in 1:ncol(td_erarc_xts)) {
  ERA5_c <- daily_climatology_filling(ts_data = td_erarc_xts[,i])
  td_eracc[,i] <- ERA5_c
}

td_eracc[sapply(td_eracc, is.nan)] <- NA
names(td_eracc) <- names(td_erarc_xts)

td_eracc_xts <- xts(td_eracc, order.by = index(td_erarc_xts))
td_eracc_list <- list(values=td_eracc_xts, xyz=td_erarc_list$xyz)

# E6: merge obs + ERA5c
td_era_cd <- as.data.frame(td_eracc_list$values)
td_data_df <- as.data.frame(readRDS('data/processed/obs/td/qc_td_obs.RDS')$values)

names(td_era_cd) <- paste0('ERA5_',names(td_era_cd))
td_era_cds <- insel_na_1(td_era_cd)

td_db <- cbind(td_data_df, td_era_cds)
td_db <- xts(td_db, order.by = as.Date(row.names(td_data_df)))

# E7: add metadata by ID
td_xyz <- readRDS('data/processed/obs/td/qc_td_obs.RDS')$xyz 

td_era_xyz <- insel_na_2(as.data.frame(td_eracc_list$values), td_eracc_list$xyz)
td_xyz$SRC <- 'OBS'
td_era_xyz$SRC <- 'ERA5'

td_xyz_db <- rbind(td_xyz, td_era_xyz)
td_xyz_db$ID <- names(td_db)
row.names(td_xyz_db) <- NULL

td_db_xts <- xts(td_db, order.by = index(td_db))
td_db_list <- list(values=td_db_xts, xyz = td_xyz_db)
saveRDS(td_db_list, 'data/processed/obs/td/qc_td_plus_era_data.RDS')
