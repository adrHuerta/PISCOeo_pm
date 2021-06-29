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
rs_ERA5 <- raster::stack('data/raw/ERA5land_rs_1981_2019.nc')
sd_list <- readRDS('data/processed/sd_wqc.RDS')
sd_xyz <- sd_list$xyz
sd_xyz <- sd_xyz[sd_xyz$QC==1,]
row.names(sd_xyz) <- NULL

rs_ERA <- as.data.frame(t(raster::extract(rs_ERA5, sd_xyz[,2:3], method='simple')))
names(rs_ERA) <- sd_xyz$ID
rs_ERA <- round(((rs_ERA/3600)/(1000000/86000)),2) #to MJ/m2 day

rs_ERA_xts <- xts(rs_ERA, order.by = index(sd_list$values))
rs_ERA <- list(values=rs_ERA_xts, xyz=sd_xyz)

# E2: convert radiation (ERA5) to sunshine duration
sd_ERA <- conv_rs_hs(rs_ERA$xyz, as.data.frame(rs_ERA$values))
sd_ERA[sd_ERA<0] <- 0
sd_ERA_xts <- xts(sd_ERA, order.by = index(rs_ERA$values))
names(sd_ERA_xts) <- paste0('ERA5_',names(rs_ERA$values))
sd_ERAm <- list(values=sd_ERA_xts, xyz=rs_ERA$xyz)

# E3: apply correction
sd_ERAm_xts <- sd_ERAm$values
sd_ERAm_id <- sd_ERAm$xyz$ID
sd_list <- readRDS('data/processed/sd_wqc.RDS')
sd_data_df <-  as.data.frame(sd_list$values)
sd_data_df <- sd_data_df[,sd_ERAm_id]
sd_data_xts <- xts(sd_data_df, order.by = index(sd_list$values))

sd_ERAc <- as.data.frame((matrix(1:nrow(sd_ERAm_xts), ncol = 1)))

for (i in 1:ncol(sd_data_xts)) {
  ERA5_c <- daily_varying_anom_qmap(ts_obs = sd_data_xts[,i], ts_model = sd_ERAm_xts[,i])
  sd_ERAc[,i] <- ERA5_c
}

names(sd_ERAc) <- names(sd_data_xts)

sd_ERAc_xts <- xts(sd_ERAc, order.by = index(sd_data_xts))
sd_ERAc_list <- list(values=sd_ERAc_xts, xyz=sd_ERAm$xyz)

# E4: ERA5 series extreme data control
sd_ERAc <- as.data.frame(sd_ERAc_list$values)
sd_ERAcq <- qc3_func_1(sd_ERAc, 3.5)
sd_ERAcq[sd_ERAcq>0] <- NA
sd_ERAcq[sd_ERAcq==0] <- 1
sd_ERAc_c <- sd_ERAcq[,1:ncol(sd_ERAc)]*sd_ERAc

sd_ERAcc_xts <- xts(sd_ERAc_c, order.by = as.Date(row.names(sd_ERAc_c)))
sd_ERAcc_list <- list(values=sd_ERAcc_xts, xyz=sd_ERAc_list$xyz)

# E5: keep if r >= 0.6 (ERA5 vs obs)
sd_list <- readRDS('data/processed/sd_wqc.RDS')

sd_ERA_xts <- sd_ERAcc_list$values
sd_ERA_id <- sd_ERAcc_list$xyz$ID
sd_data <-  as.data.frame(sd_list$values)
sd_data <- sd_data[,sd_ERA_id]
sd_data_xts <- xts(sd_data, order.by = index(sd_list$values))

sd_ERArc <- COMPAR_REO(as.data.frame(sd_ERAcc_list$values), 
                         as.data.frame(sd_data_xts), 0.55)
sd_ERArc_xts <- xts(sd_ERArc, order.by = as.Date(row.names(sd_ERA_rvs)))
sd_ERArc_list <- list(values=sd_ERArc_xts, xyz=sd_ERAcc_list$xyz)

# E6: complete with daily climatologies
sd_ERArc_xts <- sd_ERArc_list$values

sd_ERAcc <- as.data.frame((matrix(1:nrow(sd_ERArc_xts), ncol = 1)))

for (i in 1:ncol(sd_ERArc_xts)) {
  ERA5_c <- daily_climatology_filling(ts_data = sd_ERArc_xts[,i])
  sd_ERAcc[,i] <- ERA5_c
}

sd_ERAcc[sapply(sd_ERAcc, is.nan)] <- NA
names(sd_ERAcc) <- names(sd_ERArc_xts)

sd_ERAcc_xts <- xts(sd_ERAcc, order.by = index(sd_ERArc_xts))
sd_ERAcc_list <- list(values=sd_ERA_CC_xts, xyz=sd_ERArc_list$xyz)

# E7: merge obs + ERA5c
sd_ERA_cd <- as.data.frame(sd_ERAcc_list$values)
sd_data_df <- as.data.frame(readRDS('data/processed/sd_wqc.RDS')$values)

names(sd_ERA_cd) <- paste0('ERA5_',names(sd_ERA_cd))
sd_ERA_cds <- insel_na_1(sd_ERA_cd)

sd_db <- cbind(sd_data_df, sd_ERA_cds)
sd_db <- xts(sd_db, order.by = as.Date(row.names(sd_data_df)))

# E8: add metadata by ID
sd_xyz <- readRDS('data/processed/sd_wqc.RDS')$xyz 

sd_ERA_xyz <- insel_na_2(as.data.frame(sd_ERAcc_list$values), sd_ERAcc_list$xyz)
sd_xyz$SRC <- 'OBS'
sd_ERA_xyz$SRC <- 'ERA5'

sd_xyz_db <- rbind(sd_xyz, sd_ERA_xyz)
sd_xyz_db$ID <- names(sd_db)
row.names(sd_xyz_db) <- NULL

sd_db_xts <- xts(sd_db, order.by = index(sd_db))
sd_db_list <- list(values=sd_db_xts, xyz = sd_xyz_db)
saveRDS(sd_db_list, 'data/processed/sd_data_era.RDS')
