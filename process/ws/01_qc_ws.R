# ws_data_cws:  Wind speed data from conventional weather station
# qc:           quality control
# ws_xyz_cws:   Wind speed location from conventional weather station
rm(list = ls())

library(tidyverse)
library(xts)
library(raster)
library(lubridate)
library(ggpubr)
library(geosphere)

source('src/qc_functions.R')
source('src/qc_visual.R')
source('src/qc_spatial_neighbors.R')

# qc1: coding errors
ws_data <- read.csv('data/raw/obs/ws/ws_data_cws.csv', row.names = 1)
ws_qc1 <- qc1_func(ws_data, 15)

# qc2: maximum thresholds and physical limits
ws_data <- read.csv('data/raw/obs/ws/ws_data_cws.csv', row.names = 1)
ws_data[ws_data>40] <- 1000
ws_qc2 <- qc2_func(ws_data)

# qc3: percentile thresholds
ws_data <- read.csv('data/raw/obs/ws/ws_data_cws.csv', row.names = 1)
ws_qc3 <- (qc3_func_1(ws_data, 4.5))[,1:ncol(ws_data)]

# qc4: neighbors compare
ws_data <- read.csv('data/raw/obs/ws/ws_data_cws.csv', row.names = 1)
ws_xyz <- read.csv('data/raw/obs/ws/ws_xyz_cws.csv', stringsAsFactors = F)
mins_ws <- min_stt(ws_xyz)
ws_qc4 <- qc4_func(values = ws_data, xyz = ws_xyz, min_st = mins_ws, 100)

# grouped automatic qc
ws_data <- read.csv('data/raw/obs/ws/ws_data_cws.csv', row.names = 1)
ws_qcg <- ws_qc1+ws_qc2+ws_qc3+ws_qc4
ws_qcg[ws_qcg == 0] <- -1
ws_qcg[ws_qcg > 0] <- NA
ws_qcg[ws_qcg == -1] <- 1
ws_qca <- ws_data*ws_qcg

# qc5: visual control
graph_hmg(ws_qca,'data/raw/obs/ws/ws_hmg/', 'ws')

# delete year with errors
# stations that should be deleted: c("X117043", "X116060", "X116011")
ws_qcv <- delet_year(ws_qca, ws_visual_qc(ws_qca))

# minimum n years with 365 days
ws_qcmd <- qc_ny365(clima_data = ws_qcv, nyear = 3)

# minimum information 20D & 12M & n year (FLAG):"1" nyear |"0" 
ws_qcma <- qc_min_info(clima_data = ws_qcmd, nyear = 5)

# ws with qc in list: data & xyz
ws_xyz <- data.frame(ws_xyz, QC = ws_qcma$SYT)
ws_data_qcf <- sel_data_wqc(ws_qcmd)
ws_xyz_qcf <- sel_xyz_wqc(ws_qcmd, ws_xyz)
row.names(ws_xyz_qcf) <- NULL

ws_data_xts <- xts(ws_data_qcf, order.by = as.Date(row.names(ws_data_qcf)))
ws_climatology <- data.frame(t(hydroTSM::monthlyfunction(ws_data_xts, FUN=mean)))

ws_wqc <- list(values = ws_data_xts, xyz = ws_xyz_qcf)
saveRDS(ws_wqc, 'data/processed/obs/ws/qc_ws_obs.RDS')
