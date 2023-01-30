#The aim of this script is to convert the  bird density rasters to a probability 
##of occurrence. 
##
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rcartocolor, rgdal, rgeos, future, furrr, reproducible, RColorBrewer, 
               colorspace, ggspatial, ggpubr, gridExtra, hrbrthemes, terra, stringr, glue, 
               sf, tidyverse, RStoolbox, fs, future.apply, fst, trend, crayon)


g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data ---------------------------------------------------------------
root <- './outputs'
dirs <- fs::dir_ls(root, type = 'directory')
spcs <- basename(dirs)
spcs <- spcs[1:72]

edeh_shp <- sf::st_read('./inputs/Edehzhie/Edehzhie_Boundary.shp')
dehcho_shp <- sf::st_read('./inputs/Dehcho/Dehcho.shp')

targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

# Extract by mask for the study area ---------------------------------------
plot(sf::st_geometry(edeh_shp))
#use this objects for cropping and masking
edeh <- sf::st_transform(x = edeh_shp, crs = targetCRS)
dehcho <- sf::st_transform(x = dehcho_shp, crs = targetCRS)
plot(sf::st_geometry(edeh))
plot(sf::st_geometry(dehcho))

# Function -----------------------------------------------------------------
get_probOcc <- function(spc,studyAreaName){
 
  paths <- list.files(glue('./outputs_masked/{studyAreaName}'), full.names = TRUE) 
  #spc <- spcs[38] # Run and erase
  message(crayon::green("Starting with:", spc))
  fle <- grep(spc, paths, value = TRUE)
  yrs <- parse_number(basename(fle))
  yrs <- unique(yrs)
  gcm <- str_sub(basename(fle), start = 21, end = nchar(basename(fle)) - 4)
  gcm <- unique(gcm)
  
  occurRas<- map(.x = 1:length(gcm), .f = function(k){
    message(crayon::green('Loading files for', gcm[k]))
    fl <- grep(gcm[k],fle, value = TRUE)
    fl <- as.character(fl)

    proOccRas<-  map(.x = 1:length(yrs), .f = function(yr){
      message(crayon::green('Year', yrs[yr]))
      sfl <- grep(yrs[yr], fl, value = TRUE)
      rst <- raster::raster(sfl)
      dps <- raster::calc(x = rst, fun = function(pxl){1- dpois(x = 0, lambda = pxl * pi)})
      out <- glue('./outputs/occur/{studyAreaName}')
      ifelse(!dir.exists(out), dir.create(out, recursive = TRUE), print('Folder already exist'))
      raster::writeRaster(x = dps, 
                  filename = glue('{out}/occu_{studyAreaName}_{spc}_{yrs[yr]}_{gcm[k]}.tif'),
                  overwrite = TRUE)
      cat('Done!\n')
    })
  })
}

# Apply the function ------------------------------------------------------
map(.x = spcs, 'edeh', .f = get_probOcc)
