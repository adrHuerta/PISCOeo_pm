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
rs_era5 <- raster::brick('data/processed/gridded/co_variables/ERA5land_rs_1981_2019.nc') + 0
sd_list <- readRDS('data/processed/obs/sd/qc_sd_obs.RDS')
sd_xyz  <- sd_list$xyz
sd_xyz  <- sd_xyz[sd_xyz$QC==1,]
row.names(sd_xyz) <- NULL

points_xy <- raster::extract(rs_era5[[1]], sd_xyz[,2:3], method = 'simple', cellnumbers = TRUE)[,1]
rs_era <- as.data.frame(t(rs_era5[points_xy]))
names(rs_era) <- sd_xyz$ID

rs_era_xts <- xts(rs_era, order.by = index(sd_list$values))
rs_era <- list(values=rs_era_xts, xyz=sd_xyz)

# E2: convert radiation (ERA5) to sunshine duration
sd_era <- conv_rs_hs(rs_era$xyz, as.data.frame(rs_era$values))
sd_era[sd_era<0] <- 0
sd_era_xts <- xts(sd_era, order.by = index(rs_era$values))
names(sd_era_xts) <- paste0('ERA5_',names(rs_era$values))
sd_era <- list(values=sd_era_xts, xyz=rs_era$xyz)

# E3: apply correction
sd_era_xts <- sd_era$values
sd_era_id <- sd_era$xyz$ID
sd_list <- readRDS('data/processed/obs/sd/qc_sd_obs.RDS')
sd_data_df <-  as.data.frame(sd_list$values)
sd_data_df <- sd_data_df[,sd_era_id]
sd_data_xts <- xts(sd_data_df, order.by = index(sd_list$values))

sd_erac <- parallel::mclapply(1:ncol(sd_data_xts),
                              function(i){
                                
                                daily_varying_anom_qmap(ts_obs = sd_data_xts[,i], 
                                                        ts_model = sd_era_xts[,i])
                                
                              }, mc.cores = 5)

sd_erac <- do.call("cbind", sd_erac)
names(sd_erac) <- names(sd_data_xts)

sd_erac_xts <- xts(sd_erac, order.by = index(sd_data_xts))
sd_erac_list <- list(values=sd_erac_xts, xyz=sd_era$xyz)

# E4: ERA5 series extreme data control
sd_erac <- as.data.frame(sd_erac_list$values)
sd_eracq <- qc3_func_1(sd_erac, 3.5)
sd_eracq[sd_eracq>0] <- NA
sd_eracq[sd_eracq==0] <- 1
sd_erac_c <- sd_eracq[,1:ncol(sd_erac)]*sd_erac

sd_eracc_xts <- xts(sd_erac_c, order.by = as.Date(row.names(sd_erac_c)))
sd_eracc_list <- list(values=sd_eracc_xts, xyz=sd_erac_list$xyz)

# E5: keep if r >= 0.6 (ERA5 vs obs)
sd_list <- readRDS('data/processed/obs/sd/qc_sd_obs.RDS')

sd_era_xts <- sd_eracc_list$values
sd_era_id <- sd_eracc_list$xyz$ID
sd_data <-  as.data.frame(sd_list$values)
sd_data <- sd_data[,sd_era_id]
sd_data_xts <- xts(sd_data, order.by = index(sd_list$values))

sd_erarc <- compar_cc(as.data.frame(sd_eracc_list$values),
                      as.data.frame(sd_data_xts), 0.55)
sd_erarc_xts <- xts(sd_erarc, order.by = as.Date(row.names(sd_erarc)))
sd_erarc_list <- list(values=sd_erarc_xts, xyz=sd_eracc_list$xyz)

# E6: complete with daily climatologies
sd_erarc_xts <- sd_erarc_list$values

sd_eracc <- as.data.frame((matrix(1:nrow(sd_erarc_xts), ncol = 1)))

for (i in 1:ncol(sd_erarc_xts)) {
  ERA5_c <- daily_climatology_filling(ts_data = sd_erarc_xts[,i])
  sd_eracc[,i] <- ERA5_c
}

sd_eracc[sapply(sd_eracc, is.nan)] <- NA
names(sd_eracc) <- names(sd_erarc_xts)

sd_eracc_xts <- xts(sd_eracc, order.by = index(sd_erarc_xts))
sd_eracc_list <- list(values=sd_eracc_xts, xyz=sd_erarc_list$xyz)

# E7: merge obs + ERA5c
sd_era_cd <- as.data.frame(sd_eracc_list$values)
sd_data_df <- as.data.frame(readRDS('data/processed/obs/sd/qc_sd_obs.RDS')$values)

names(sd_era_cd) <- paste0('ERA5_',names(sd_era_cd))
sd_era_cds <- insel_na_1(sd_era_cd)

sd_db <- cbind(sd_data_df, sd_era_cds)
sd_db <- xts(sd_db, order.by = as.Date(row.names(sd_data_df)))

# E8: add metadata by ID
sd_xyz <- readRDS('data/processed/obs/sd/qc_sd_obs.RDS')$xyz 

sd_era_xyz <- insel_na_2(as.data.frame(sd_eracc_list$values), sd_eracc_list$xyz)
sd_xyz$SRC <- 'OBS'
sd_era_xyz$SRC <- 'ERA5'

sd_xyz_db <- rbind(sd_xyz, sd_era_xyz)
sd_xyz_db$ID <- names(sd_db)
row.names(sd_xyz_db) <- NULL

sd_db_xts <- xts(sd_db, order.by = index(sd_db))
sd_db_list <- list(values=sd_db_xts, xyz = sd_xyz_db)
saveRDS(sd_db_list, 'data/processed/obs/sd/qc_sd_plus_era_data.RDS')
