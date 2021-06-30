conv_rs_hs <- function(UBIC_EHS, RADS_ERA){
  HSOL_ERA <- RADS_ERA
  for (e in 1:nrow(UBIC_EHS)) {
    lat <- UBIC_EHS[e,3]
    rad_est <- RADS_ERA[,e]
    diaj <- yday(as.Date(rownames(RADS_ERA)))
    latr <- (pi/180*lat)#latitud en radianes
    drts <- (1+0.033*cos(2*pi/365*diaj))#distancia relativa inversa Tierra-Sol
    decs <- 0.409*sin((2*pi*(diaj/365))-1.39)
    arad <- acos(-tan(latr)*tan(decs))#angulo de rad a la hora puesta del sol
    ang1 <- sin(latr)*sin(decs)
    ang2 <- cos(latr)*cos(decs)
    N <- 24/pi*arad
    RADX <- ((24*60/pi)*(0.082*drts)*(arad*(ang1)+(ang2)*sin(arad)))
    HSOL_C <- round(((((rad_est/RADX)-0.25)/0.5)*N),2)
    HSOL_r <- ifelse(HSOL_C>N, N, HSOL_C)
    HSOL_ERA[,e] <- HSOL_r
  }
  return(HSOL_ERA)
}

conv_hs_rs <- function(UBIC_EHS, HSOL_DATA){
  RAD_DATA <- HSOL_DATA
  for (e in 1:nrow(UBIC_EHS)) {
    lat <- UBIC_EHS[e,3]
    rad_est <- HSOL_DATA[,e]
    diaj <- yday(as.Date(rownames(HSOL_DATA)))
    latr <- pi/180*lat #latitud en radianes
    drts <- 1+0.033*cos(2*pi/365*diaj) #distancia relativa inversa Tierra-Sol
    decs <- 0.409*sin((2*pi*(diaj/365))-1.39) #declinacion solar
    arad <- acos(-tan(latr)*tan(decs)) #angulo de rad a la hora puesta del sol
    ang1 <- sin(latr)*sin(decs)
    ang2 <- cos(latr)*cos(decs)
    N <- 24/pi*arad
    RADX <- (24*60/pi)*(0.082*drts)*(arad*(ang1)+(ang2)*sin(arad))
    
    RADS <- round(((0.25+0.5*(rad_est/N))*RADX),2)
    RAD_DATA[,e] <- RADS
  }
  return(RAD_DATA)
}

insel_na_1 <- function(HSOL){
  metad <- data.frame()
  for (i in 1:ncol(HSOL)) {
    data <- ifelse(identical(HSOL[,i], as.numeric(rep(NA, length(HSOL[,i])))), NA, 1)
    metad <- rbind(metad, data)
  }
  names(metad) <- 'NST'
  metad$ESTAC <- names(HSOL)
  HSOL_T <- as.data.frame(t(as.matrix(HSOL)))
  HSOL_F <- HSOL_T[!is.na(metad$NST), ]
  HSOL_F2 <- as.data.frame(t(as.matrix(HSOL_F)))
  return(HSOL_F2)
}

insel_na_2 <- function(RADS_ERA_CD, UBIC_RADS_OB){
  metad <- data.frame()
  for (i in 1:ncol(RADS_ERA_CD)) {
    data <- ifelse(identical(RADS_ERA_CD[,i], as.numeric(rep(NA, length(RADS_ERA_CD[,i])))), NA, 1)
    metad <- rbind(metad, data)
  }
  names(metad) <- 'NST'
  metad$ESTAC <- names(RADS_ERA_CD)
  HSOL_F <- UBIC_RADS_OB[!is.na(metad$NST), ]
  return(HSOL_F)
}

insel_na_3 <- function(RADS_ERA_CD){
  RADS_ERA_t <- as.data.frame(t(as.matrix(RADS_ERA_CD)))
  RADS_ERA_na <- na.omit(RADS_ERA_t)
  RADS_ERA_t2 <- as.data.frame(t(as.matrix(RADS_ERA_na)))
  return(RADS_ERA_t2)
}

compar_cc <- function(clima_era, clima_obs, r){
  clima_era_c <- clima_era
  cc_n <- data.frame()
  for (i in 1:ncol(clima_era)) {
    vect_ob <- cor(clima_era[,i], clima_obs[,i], use="pairwise.complete.obs")
    cc_n <- rbind(cc_n, vect_ob)
  }
  cc_n <- round(cc_n,2)
  names(cc_n) <- 'R_EO'
  cc_n[cc_n>=r] <- 1
  cc_n[cc_n<1] <- NA
  
  for (i in 1:nrow(cc_n)) {
    if(is.na((cc_n[,1])[i])){
      clima_era_c[,i] <- NA
    }
  }
  return(clima_era_c)
}
