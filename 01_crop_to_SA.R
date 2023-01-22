# Load libraries --------------------------------------------------------
require(pacman)

pacman::p_load(dplyr, glue, sf, terra, tidyverse)

g <- gc(reset = TRUE)
rm(list = ls())

# Load data ---------------------------------------------------------------
root <- './outputs'
dirs <- fs::dir_ls(root, type = 'directory')
spcs <- basename(dirs)
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

# Clip rasters to study area ----------------------------------------------

rasterToTable<- function(spc,studyArea,studyAreaName){
 # spc <- spcs[1] # for proof. Run an erase or comment out 
  message(crayon::green('Loading data for: ', spc))
  dir <- grep(spc, dirs, value = TRUE)
  fls <- fs::dir_ls(dir, regexp = '.tif$')
  fls <- grep('mean', fls, value = TRUE)
  yrs <- parse_number(basename(fls))
  yrs <- unique(yrs)
  yrs <- na.omit(yrs)
  gcm <- str_sub(basename(fls), start = 16, end = nchar(basename(fls)) - 4)
  gcm <- unique(gcm)
  
  cat('Raster to table\n')
  message(crayon::green('Cropping raster for: ', spc))
  dfm <- map(.x = 1:length(gcm), .f = function(k){
    
    cat(gcm[k], '\n')
    fl <- grep(gcm[k], fls, value = TRUE)
  
      cropRas <- map(.x = 1:length(yrs), .f = function(yr){ 
        message(crayon::green('Year', yrs[yr]))
        sfl <- grep(yrs[yr], fl, value = TRUE)
        tr <- terra::rast(sfl)
        rs <- terra::crop(tr, studyArea)
        rs <- terra::mask(rs, studyArea)
        names(rs) <- paste0('y', yrs[yr])
        out <- glue('./outputs_masked/{studyAreaName}')
        ifelse(!dir.exists(out), dir.create(out, recursive = TRUE), print('Folder already exist'))
        terra::writeRaster(rs, filename = glue('{out}/mean_mask_{spc}_{yrs[yr]}_{gcm[k]}.tif'),
                           overwrite = TRUE)
        message(crayon::green('Converting raster to table'))
        df <- terra::as.data.frame(rs, xy = TRUE, na.rm = TRUE)
        df <- as_tibble(df)
        df <- mutate(df, gc = gcm[k])
        return(df)
      })
      fdf <- Reduce(full_join, cropRas)
   })
    rsl <- bind_rows(dfm)
    rsl <- rsl %>% relocate('y2011', .after = gc)
    out <- glue('./qs_{studyAreaName}')
    ifelse(!dir.exists(out), dir.create(out, recursive = TRUE), print('Folder already exist'))
    qs::qsave(x = rsl, file = glue('{out}/tbl_yrs_{studyAreaName}_{spc}.qs'))
  
  cat('------- Done -------\n')
  return(rsl)
}
### Raster to table ---------------------------------------------------------
dfrm <- map(.x = spcs, edeh,'edeh' , .f = rasterToTable)
  
  

