# td_data_cws:  Sunshine duration data from conventional weather station
# qc:           quality control
# td_xyz_cws:   Sunshine duration location from conventional weather station

library(tidyverse)
library(xts)
library(raster)
library(lubridate)
library(ggpubr)
library(geosphere)

source('src/qc_functions.R')
source('src/qc_visual.R')
source('src/qc_spatial_neighbors.R')

# qc1: coding errors (7)
td_data <- read.csv('data/raw/obs/td/td_data_cws.csv', row.names = 1)
td_qc1    <- qc1_func(td_data, 7)

# qc2: maximum thresholds and physical limits
td_data <- read.csv('data/raw/obs/td/td_data_cws.csv', row.names = 1)
td_data[td_data > 50 | td_data < -35] <- 1000
td_qc2  <- qc2_func(td_data)

# qc3: percentile thresholds
td_data <- read.csv('data/raw/obs/td/td_data_cws.csv', row.names = 1)
td_qc3  <- (qc3_func_1(td_data, 3.5))[,1:ncol(td_data)]

# qc4: neighbors compare
td_data <- read.csv('data/raw/obs/td/td_data_cws.csv', row.names = 1)
td_xyz  <- read.csv('data/raw/obs/td/td_xyz_cws.csv', stringsAsFactors = F)
mins_sd <- min_stt(td_xyz)
td_qc4  <- qc4_func(values = td_data, xyz = td_xyz, min_st = mins_sd, 80)

# grouped automatic qc
td_data <- read.csv('data/raw/obs/td/td_data_cws.csv', row.names = 1)
td_qcg  <- td_qc1+td_qc2+td_qc3+td_qc4
td_qcg[td_qcg == 0]   <- -1
td_qcg[td_qcg > 0]    <- NA
td_qcg[td_qcg == -1]  <- 1
td_qca <- td_data*td_qcg

# qc5: visual control
graph_hmg(td_qca,'data/raw/obs/td/td_hmg/', 'td')

# delete year with errors
td_qcv <- delet_year(td_qca, td_visual_qc(td_qca))

# minimum n years with 365 days
td_qcmd <- qc_ny365(clima_data = td_qcv, nyear = 5)

# minimum information 20D & 12M & n year (FLAG):"1" nyear |"0" 
td_qcma <- qc_min_info(clima_data = td_qcmd, nyear = 7)

# sd with qc in list: data & xyz
td_xyz      <- data.frame(td_xyz, QC = td_qcma$SYT)
td_data_qcf <- sel_data_wqc(td_qcmd)
td_xyz_qcf  <- sel_xyz_wqc(td_qcmd, td_xyz)
row.names(td_xyz_qcf) <- NULL

td_data_xts <- xts(td_data_qcf, order.by = as.Date(row.names(td_data_qcf)))
td_climatology <- data.frame(t(hydroTSM::monthlyfunction(td_data_xts, FUN=mean)))

td_wqc <- list(values=td_data_xts, xyz=td_xyz_qcf, climat = td_climatology)

saveRDS(td_wqc, 'data/processed/obs/td/qc_td_obs.RDS')
