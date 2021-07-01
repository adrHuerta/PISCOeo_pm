# http://www.fao.org/3/x0490s/x0490s.pdf 
# page 48, eq. 34

maximum_lenght_sd <- function(jday_i, lat_i)
{
  lat_i <- lat_i*pi/180
  DELTA_solar_declintacion <- 0.409 * sin( (jday_i * ((2*pi)/365)) - 1.39 )
  W_s <- raster::calc(-raster::calc(lat_i, tan) * tan(DELTA_solar_declintacion), acos)
  (24/pi)*W_s
}
