# Load libraries --------------------------------------------------------
require(pacman)

pacman::p_load(raster, rgdal, rgeos, terra, stringr, glue, sf, tidyverse, RStoolbox, fs, fst, trend)

g <- gc(reset = TRUE)
rm(list = ls())

# Load data ---------------------------------------------------------------
root <- './outputs'
dirs <- fs::dir_ls(root, type = 'directory')
spcs <- basename(dirs)
spcs <- spcs[1:72]


# See the changes  --------------------------------------------------------
raster_to_table <- function(spc, studyAreaName){
  
  # Proof
  #spc <- spcs[2] # Run and comment (after)
  message(crayon::green("Starting with:", spc))
  path <- glue('./outputs/occur/{studyAreaName}')
  fls <- list.files(path, pattern = 'occu', full.names = TRUE)
  fls <- grep(spc,fls, value = TRUE)
  yrs <- parse_number(basename(fls))
  yrs <- unique(yrs)
  gcm <- str_sub(basename(fls), start = 21, end = nchar(basename(fls)) - 4) # start = 23 for dehcho area and 21 for edeh
  gcm <- unique(gcm)
  
  cat('Raster to table\n')
  dfm <- map(.x = 1:length(gcm), .f = function(k){
    
    message(crayon::green(gcm[k]))
    fl <- grep(gcm[k], fls, value = TRUE)
    rs <- terra::rast(fl)
    df <- terra::as.data.frame(rs, xy = TRUE, na.rm = TRUE)
    colnames(df)[3:8]<- glue('y{yrs}')
    df <- as_tibble(df)
    df <- df %>% mutate(gc = gcm[k],
                        specie = spc)
    return(df)
  })
  
  rsl <- bind_rows(dfm)
  out <- glue('./qs_{studyAreaName}/occur')
  ifelse(!dir.exists(out), dir.create(out, recursive = TRUE), print('Folder already exist'))
  qs::qsave(x = rsl, file = glue('{out}/occ_yrs_{spc}.qs'))
  
  cat('------- Done -------\n')
  return(rsl)
  
}
### Raster to table ---------------------------------------------------------
dfrm <- map(.x = spcs,'edeh', .f = raster_to_table)
