yeardays <- function(time){
  time <- as.POSIXlt(time)
  time$mon[] <- time$mday[] <- time$sec[] <- time$min <- time$hour <- 0
  time$year <- time$year + 1
  return(as.POSIXlt(as.POSIXct(time))$yday + 1)
}

year2day <- function(clim_year){
  fyears <- row.names(clim_year)
  yrs <- as.data.frame(t(matrix(1:ncol(clim_year), ncol = 1)))
  colnames(yrs) <- colnames(clim_year)
  for (i in 1:length(fyears)) {
    pc1 <- do.call("rbind", replicate(yeardays(fyears[i]), clim_year[i,], simplify = FALSE))
    yrs <- rbind(yrs, pc1)
  }
  yrs_f <- yrs[2:nrow(yrs),]
  return(yrs_f)
}

qc1_func <- function(val_clim, umbral){
  recop_temp <- val_clim
  for (i in 1:ncol(val_clim)) {
    repx <- rle(val_clim[,i])[[1]]
    rep_cs <- cumsum(repx)
    rep_df <- data.frame(rlev= repx, acum = rep_cs)
    rep_df$rlev[rep_df$rlev <= umbral] <- NA 
    rep_df <- na.omit(rep_df)
    for (y in 1:nrow(rep_df)) {
      ifelse(nrow(rep_df)==0, rep_df<-rep_df,recop_temp[((rep_df[y,2]-rep_df[y,1]+1):rep_df[y,2]), i] <- 1000)
    }
  }
  recop_temp[recop_temp<1000] <- 0
  recop_temp[is.na(recop_temp)] <- 0
  recop_temp[recop_temp==1000] <- 1
  return(recop_temp)
}

qc2_func <- function(val_clim){
  val_clim[val_clim != 1000] <- 0
  val_clim[val_clim == 1000] <- 1
  val_clim[is.na(val_clim)] <- 0
  return(val_clim)
}

qc3_func_2 <- function(clim_eval, val_clim, mlt){
  for (i in 1:ncol(clim_eval)) {
    q25p <- as.numeric(quantile(val_clim[,i],0.25, na.rm = T)-as.numeric(mlt)*IQR(val_clim[,i], na.rm = T))
    q75p <- as.numeric(quantile(val_clim[,i],0.75, na.rm = T)+as.numeric(mlt)*IQR(val_clim[,i], na.rm = T))
    clim_eval[,i][clim_eval[,i] > q75p  | clim_eval[,i] < q25p] <- 1000
  }
  clim_eval[clim_eval != 1000] <- 0
  clim_eval[clim_eval == 1000] <- 1
  clim_eval[is.na(clim_eval)] <- 0
  return(clim_eval)
}

qc3_func_1 <- function(val_clim, mlt){
  clim_sq_ts <- dplyr::mutate(val_clim, MES=as.numeric(format(as.Date(row.names(val_clim),format="%Y-%m-%d"),format="%m")))
  q3_temp_fn <- data.frame()
  for (yr in 1:12) {
    clim_01 <- subset(clim_sq_ts, MES == yr)
    clim_eval <- clim_01
    q3_temp_01 <- qc3_func_2(clim_eval, clim_01, mlt)
    q3_temp_fn <- rbind(q3_temp_fn, q3_temp_01)
  }
  q3_temp_fn2 <- xts(x = q3_temp_fn, order.by = as.Date(row.names(q3_temp_fn)))
  q3_temp_fn2 <- as.data.frame(q3_temp_fn2)
  return(q3_temp_fn2)
}

min_stt <- function(xyz_data){
  minim <- data.frame(1:5)
  for (i in 1:nrow(xyz_data)) {
    est_min2 <- data.frame(EST = spt_neighrs(id_station = xyz_data[i,1],
                                             stations_database = xyz_data,
                                             lmt_dist = 70,
                                             lmt_elv = 1000,
                                             lmt_n = 4)) 
    nr <- nrow(est_min2)
    est_min2 <- data.frame(ifelse(nrow(est_min2)!=5, rbind(est_min2, data.frame(EST = rep(xyz_data[i,1], (5-nr)))), est_min2))
    names(est_min2) <- 'EST'
    minim[,i] <- est_min2[,1]
  }
  colnames(minim) <- paste0('EST_', 1:nrow(xyz_data))
  
  level_vec <- xyz_data$ID
  replacement_vec <- paste(1:nrow(xyz_data))
  minim[] <- lapply(minim, function(x) forcats::lvls_revalue(factor(x, levels = level_vec),replacement_vec))
  
  min_station <- minim
  for (i in 1:nrow(xyz_data)) {
    min_station[ , i] <- as.numeric(as.character(min_station[ , i]))
  }
  
  min_station <- min_station[2:5,]
  return(min_station)
}

percentils <- function(vect){
  pt1 <- quantile(vect, probs = seq(0, 1, by = 0.01), type = 7, na.rm=T)
  pt2 <- unique(as.data.frame(pt1), fromLast = TRUE)
  pt3 <- rownames(pt2)
  pt4 <- as.integer(strsplit(pt3, "%"))
  dato <- pt4[as.integer(cut(vect, c(-Inf, pt2$pt1), labels = 1:length(pt3)))]
  return(dato)
}

mean_cn <- function(a,n){
  count <- sum(!is.na(a))
  if (count >= n) {
    tot <- mean(a, na.rm=T)
  } else tot <- NA
  return (tot)
}

mean_n <- function(a,n){
  count <- sum(is.na(a))
  if (count <= n) {
    tot <- mean(a, na.rm=T)
  } else tot <- NA
  return (tot)
}

qc4_func <- function(values = ws_data, xyz = ws_xyz, min_st = mins_ws, umb = 100){
  est_wq4 <- data.frame(1:nrow(values))
  for (w in 1:nrow(xyz)) {
    minim_4e <- min_st[,w]
    vect_ob <- ifelse(class(values[,w])=='numeric', percentils(values[,w]), NA)
    vect_vc <- data.frame(p1 = ifelse(class(values[,minim_4e[1]])=='numeric', percentils(values[,minim_4e[1]]), NA),
                          p2 = ifelse(class(values[,minim_4e[2]])=='numeric', percentils(values[,minim_4e[2]]), NA),
                          p3 = ifelse(class(values[,minim_4e[3]])=='numeric', percentils(values[,minim_4e[3]]), NA),
                          p4 = ifelse(class(values[,minim_4e[4]])=='numeric', percentils(values[,minim_4e[4]]), NA))
    vect_vc <- dplyr::mutate(vect_vc, PM = ifelse(is.nan(rowMeans(dplyr::select(vect_vc, p1,p2,p3,p4), na.rm = TRUE)), NA, 
                                                  rowMeans(dplyr::select(vect_vc, p1,p2,p3,p4), na.rm = TRUE)))
    comp_vec <- data.frame(OB=vect_ob, VC=vect_vc[,5])
    comp_vec <- dplyr::mutate(comp_vec, DIF = comp_vec[,1]-comp_vec[,2])
    est_wq4[,w] <- comp_vec[,3]
  }
  est_wq4[est_wq4 == 'NaN'] <- NA
  est_wq4[is.na(est_wq4)] <- 0
  est_wq4_abs <- abs(est_wq4)
  est_wq4_abs[est_wq4_abs<umb] <- 0
  est_wq4_abs[est_wq4_abs>0] <- 1
  row.names(est_wq4_abs) <- row.names(values)
  return(est_wq4_abs)
}

graph_hmg <- function(clim_data, path, tipdata){
  for (i in 1:ncol(clim_data)) {
    if(identical(clim_data[,i], rep(NA, length(clim_data[,i])))) {
      
    } else {
      clim_dc <- round(abs(round(clim_data[,i],1) - trunc(round(clim_data[,i],1))), 1)
      clim_df <- data.frame(prc = clim_dc, dates = rownames(clim_data),
                            year = format(as.Date(row.names(clim_data),format="%Y-%m-%d"),format="%Y"))
      clim_df$prc[clim_df$prc==0] <- 1
      breaks <- c(seq(0.05, 0.95, by = 0.1),1.1)
      clim_df$rnk <- cut(clim_df$prc,breaks = breaks, right = FALSE,
                         labels=c(seq(0.1, 0.9, by = 0.1),1))
      clim_df$rnk <- factor(clim_df$rnk, levels = c('1','0.1','0.2','0.3','0.4',
                                                    '0.5','0.6','0.7','0.8','0.9'),
                            labels = c('0', '0.1','0.2','0.3','0.4',
                                       '0.5','0.6','0.7','0.8','0.9'))
      
      graf <- ggplot(clim_df, aes(x=year, y=prc))+
        geom_bar(aes(fill = factor(rnk, levels=rev(levels(rnk)))), stat="identity")+
        scale_fill_manual("legend", values = c("0" = "black", "0.1" = "#7C766F", 
                                               "0.2" = "darkred", "0.3" = "red",
                                               "0.4" = "orange", "0.5" = "gold",
                                               "0.6" = "yellow", "0.7" = "greenyellow",
                                               "0.8" = "green", "0.9" = "forestgreen")) +
        guides(fill = guide_legend(title.position = "right",
                                   reverse=F,title = 'Decimals',
                                   title.vjust = 0.5))+
        theme_bw()+
        theme(panel.background = element_rect(fill = "transparent"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_line(color = 'gray40', linetype = 'dashed', size = 0.1),
              axis.text.y = element_text(color = 'black',size=10, hjust=1),
              axis.text.x = element_text(color = 'black',size=10, angle = 90, vjust=0.5),
              axis.title = element_blank(),
              legend.position="bottom")

      
      time_series <- ggplot(clim_data, aes(x=as.Date(row.names(clim_data)), y=clim_data[,i]))+
        scale_x_date(date_breaks = "1 year", date_labels = "%Y",
                     limit = c(as.Date(row.names(clim_data)[1]), as.Date(row.names(clim_data)[nrow(clim_data)])))+
        geom_point(color='darkblue', cex=0.2)+
        ggtitle(paste0(tipdata,'_station_',colnames(clim_data)[i]))+
        theme_bw() +
        theme(panel.background = element_rect(fill = "transparent"),
              panel.grid.minor = element_line(linetype = "dotted"),
              panel.grid.major = element_line(color = 'gray40', linetype = 'dashed', size = 0.1),
              axis.text.y = element_text(color = 'black',size=10, hjust=1),
              axis.text.x = element_text(color = 'black',size=10, angle = 90, vjust=0.5),
              axis.title = element_blank(),
              legend.position="none")
      
      clim_st <- xts(clim_data, order.by=as.Date(row.names(clim_data)))
      clim_st <- apply.monthly(clim_st, MARGIN=2, FUN=apply, mean_cn, n=20)
      clim_st <- as.data.frame(clim_st)
      
      time_series_m <- ggplot(clim_st, aes(x=as.Date(row.names(clim_st)), y=clim_st[,i]))+
        scale_x_date(date_breaks = "1 year", date_labels = "%Y",
                     limit = c(as.Date(row.names(clim_st)[1]), as.Date(row.names(clim_st)[nrow(clim_st)])))+
        geom_line(color='darkblue', cex=0.5)+
        ggtitle(paste0(tipdata,'_station_',colnames(clim_st)[i]))+
        theme_bw() +
        theme(panel.background = element_rect(fill = "transparent"),
              panel.grid.minor = element_line(linetype = "dotted"),
              panel.grid.major = element_line(color = 'gray40', linetype = 'dashed', size = 0.1),
              axis.text.y = element_text(color = 'black',size=10, hjust=1),
              axis.text.x = element_text(color = 'black',size=10, angle = 90, vjust=0.5),
              axis.title = element_blank(),
              legend.position="none")
      
      plotcomb <- ggarrange(time_series, time_series_m, graf,
                            labels = c("a)", "b)", 'c)'),vjust = 4,
                            ncol = 1, nrow = 3, font.label=list(size=12,face="bold", color="black")) 
      
      ggsave(plot=plotcomb, paste0(path, tipdata,'_',i,'_',colnames(clim_data)[i],'.png'), 
             units = "mm", width = 200, height = 200, dpi = 300)
    }
  }
}

delet_year <- function(clima_data, dlt_station){
  names(dlt_station) <- names(clima_data)
  row.names(dlt_station) <- seq(ymd(as.Date('1981-01-01')), ymd(as.Date('2019-12-31')), by='year')
  clima_dts <- year2day(dlt_station)
  clima_fmi <- clima_data*clima_dts
  return(clima_fmi)
}

qc_ny365 <- function(clima_data = ws_data, nyear=5){
  clima_ts <- xts(x = clima_data, order.by = as.Date(row.names(clima_data)))
  clima_yts <- as.data.frame(apply.yearly(clima_ts, FUN=apply, MARGIN=2, mean_cn, n=365))
  clima_yts[clima_yts < 1000] <- 1
  clima_yts_rs <- data.frame(syt = colSums(clima_yts, na.rm = T))
  clima_yts_rs[clima_yts_rs < nyear] <- NA
  clima_yts_rs[clima_yts_rs > 0] <- 1
  for (i in 1:nrow(clima_yts_rs)) {
    if(is.na((clima_yts_rs$syt)[i])){
      clima_data[,i] <- NA
    }else{
      
    }
  }
  return(clima_data)
}

qc_min_info <- function(clima_data = ws_data, nyear = 5){
  clima_ts <- xts(x = clima_data, order.by = as.Date(row.names(clima_data)))
  clima_mts <- as.data.frame(apply.monthly(clima_ts,FUN=apply,MARGIN=2, mean_cn, n=20))
  clima_yts <- as.data.frame(apply.yearly(clima_mts,FUN=apply,MARGIN=2, mean_n, n=0))
  clima_yts[clima_yts < 1000] <- 1
  
  clima_yts_rs <- data.frame(SYT = colSums(clima_yts, na.rm = T))
  clima_yts_rs[clima_yts_rs < nyear] <- 0
  clima_yts_rs[clima_yts_rs > 0] <- 1
  return(clima_yts_rs)
}

sel_data_wqc <- function(clima_data){
  metad <- data.frame()
  for (i in 1:ncol(clima_data)) {
    data <- ifelse(identical(clima_data[,i], rep(NA, length(clima_data[,i]))), NA, 1)
    metad <- rbind(metad, data)
  }
  names(metad) <- 'NST'
  metad$ESTAC <- names(clima_data)
  clima_t <- as.data.frame(t(as.matrix(clima_data)))
  clima_f <- clima_t[!is.na(metad$NST), ]
  clima_f2 <- as.data.frame(t(as.matrix(clima_f)))
  return(clima_f2)
}

sel_xyz_wqc <- function(clima_data){
  metad <- data.frame()
  for (i in 1:ncol(clima_data)) {
    data <- ifelse(identical(clima_data[,i], rep(NA, length(clima_data[,i]))), NA, 1)
    metad <- rbind(metad, data)
  }
  names(metad) <- 'NST'
  metad$ESTAC <- names(clima_data)
  return(metad)
}

sel_xyz_wqc <- function(clima_data, clima_xyz){
  metad <- data.frame()
  for (i in 1:ncol(clima_data)) {
    data <- ifelse(identical(clima_data[,i], rep(NA, length(clima_data[,i]))), NA, 1)
    metad <- rbind(metad, data)
  }
  names(metad) <- 'NST'
  metad$ESTAC <- names(clima_data)
  xyz_cf <- clima_xyz[!is.na(metad$NST), ]
  return(xyz_cf)
}
