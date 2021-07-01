checking_files <- function(var_grid, time_range)
  {
  
  # checking files
  path_gridded_sub_variables <- "./data/processed/gridded/sub_variables/values"
  files_to_check <- sprintf("%s_%s.nc",
                            var_grid,
                            seq(as.Date(time_range[1]), as.Date(time_range[length(time_range)]), by = "day"))
  
  are_files <- files_to_check %in% dir(file.path(path_gridded_sub_variables, var_grid))
  

  if(all(isTRUE(are_files)) == FALSE){
    
    return(list(files_var = files_to_check[are_files %in% FALSE],
                files_id = match(files_to_check[are_files %in% FALSE], files_to_check)))
    
  }
  
}


get_pixel_ts_from_grids <- function(var_grid, pixel = 5000,
                                    n_cores_ = 4)
{
  
  list_of_grids <- dir(file.path("./data/processed/gridded/sub_variables/values", var_grid), full.names = TRUE)
  ts_of_grids <- parallel::mclapply(list_of_grids,
                                    function(x) raster::raster(x)[pixel],
                                    mc.cores = n_cores_)
  unlist(ts_of_grids)
}
