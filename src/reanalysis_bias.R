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
