# sd_data_cws:  Sunshine duration data from conventional weather station
# qc:           quality control
# sd_xyz_cws:   Sunshine duration location from conventional weather station

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
sd_data <- read.csv('data/raw/sd_data_cws.csv', row.names = 1)
#0 (15)
sd_sqc0 <- sd_data
sd_sqc0[sd_sqc0>0] <- NA
sd_0    <- qc1_func(sd_sqc0, 15)

#>0 (10)
sd_sqc0m <- sd_data
sd_sqc0m[sd_sqc0m==0] <- NA
sd_0m   <- qc1_func(sd_sqc0m, 10)

sd_qc1  <- sd_0+sd_0m

# qc2: maximum thresholds and physical limits
sd_data <- read.csv('data/raw/sd_data_cws.csv', row.names = 1)
sd_data[sd_data>13 | sd_data<0] <- 1000
sd_qc2  <- qc2_func(sd_data)

# qc3: percentile thresholds
sd_data <- read.csv('data/raw/sd_data_cws.csv', row.names = 1)
sd_qc3  <- (qc3_func_1(sd_data, 3.5))[,1:ncol(sd_data)]

# qc4: neighbors compare
sd_data <- read.csv('data/raw/sd_data_cws.csv', row.names = 1)
sd_xyz  <- read.csv('data/raw/sd_xyz_cws.csv', stringsAsFactors = F)
mins_sd <- min_stt(sd_xyz)
sd_qc4  <- qc4_func(values = sd_data, xyz = sd_xyz, min_st = mins_sd, 80)
names(sd_qc4) <- names(sd_data)

# grouped automatic qc
sd_data <- read.csv('data/raw/sd_data_cws.csv', row.names = 1)
sd_qcg  <- sd_qc1+sd_qc2+sd_qc3+sd_qc4
sd_qcg[sd_qcg == 0] <- -1
sd_qcg[sd_qcg > 0] <- NA
sd_qcg[sd_qcg == -1] <- 1
sd_qca <- sd_data*sd_qcg

# qc5: visual control
graph_hmg(sd_qca,'data/raw/graphics/sd_hmg/', 'sd')

# delete year with errors
sd_qcv <- delet_year(sd_qca, sd_visual_qc(sd_qca))

# minimum n years with 365 days
sd_qcmd <- qc_ny365(clima_data = sd_qcv, nyear = 5)

# minimum information 20D & 12M & n year (FLAG):"1" nyear |"0" 
sd_qcma <- qc_min_info(clima_data = sd_qcmd, nyear = 7)

# sd with qc in list: data & xyz
sd_xyz <- data.frame(sd_xyz, QC = sd_qcma$SYT)
sd_data_qcf <- sel_data_wqc(sd_qcmd)
sd_xyz_qcf <- sel_xyz_wqc(sd_qcmd, sd_xyz)
row.names(sd_xyz_qcf) <- NULL

sd_data_xts <- xts(sd_data_qcf, order.by = as.Date(row.names(sd_data_qcf)))
sd_climatology <- data.frame(t(hydroTSM::monthlyfunction(sd_data_xts, FUN=mean)))

sd_wqc <- list(values=sd_data_xts, xyz=sd_xyz_qcf, climat = sd_climatology)

saveRDS(sd_wqc, 'data/processed/sd_wqc.RDS')
