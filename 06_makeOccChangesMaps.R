# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(ggspatial, 
               ggpubr, gridExtra, stringr, glue, sf, tidyverse, fasterize,RColorBrewer,
               RStoolbox, fs, fst, trend, colorspace, hrbrthemes,exactextractr, purrr, future, spatialEco)

g <- gc(reset = TRUE)
rm(list = ls())

# Load data ---------------------------------------------------------------4
species <- c("ALFL", "AMCR", "AMRE", "AMRO", "ATSP", "ATTW", "BAWW", "BBWA", 
             "BBWO", "BCCH", "BHCO", "BHVI", "BLPW", "BOCH", "BRBL", "BRCR",
             "BTNW", "CAWA", "CCSP", "CHSP", "CONW", "CORA", "COYE", "DEJU", 
             "EAKI", "EVGR", "FOSP", "GCKI", "GCTH", "GRAJ", "HAFL", "HETH",
             "HOLA", "LCSP", "LEFL", "LEYE", "LISP", "MAWA", "NOWA", "OCWA", 
             "OSFL", "OVEN", "PAWA", "PHVI", "PISI", "PUFI", "RBGR", "RBNU", 
             "RCKI", "REVI", "RUBL", "RUGR", "RWBL", "SAVS", "SOSP", "SWSP", 
             "SWTH", "TEWA" ,"TRES", "VATH", "WAVI", "WCSP", "WETA", "WEWP", 
             "WIWA", "WIWR", "WTSP", "WWCR", "YBFL", "YBSA", "YEWA", "YRWA")


edeh_shp <- sf::st_read('./inputs/Edehzhie/Edehzhie_Boundary.shp')
dehcho_shp <- sf::st_read('./inputs/Dehcho/Dehcho.shp')

targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")


# Extract by mask for the study area ---------------------------------------
edeh <- sf::st_transform(x = edeh_shp, crs = targetCRS)
dehcho <- sf::st_transform(x = dehcho_shp, crs = targetCRS)
plot(sf::st_geometry(edeh))
plot(sf::st_geometry(dehcho))

# Function to use ---------------------------------------------------------
changes_rasters <- function(spc, studyArea, studyAreaName, yr1, yr2){
  #spc <- species[2] # Run and comment (after)
  message(crayon::green("Calculating changes for: ", spc))
  path <- glue('./qs_{studyAreaName}/occur')
  fls <- list.files(path, pattern =  'occ_yrs', full.names = TRUE)
  fle <- grep(spc, fls, value = TRUE)
  tbl <- qs::qread(file = glue('{path}/occ_yrs_{spc}.qs'))
  tbl <- dplyr::select(tbl, x, y, gc, everything())
  names(tbl)[1:2] <- c('lon', 'lat')
  tbl <- mutate(tbl, avg = rowMeans(tbl[,4:9]))
  tbl <- as_tibble(tbl)
  gcm <- unique(tbl$gc)
  
  message(crayon::green(glue('Estimating change from {yr1} to {yr2} year\n')))
  tbl <- mutate(tbl, change = tbl$y2091 - tbl$y2011)  #change 2100 for 2091
  tbl <- mutate(tbl, perctChange = ((tbl$y2091 - tbl$y2011) / tbl$y2011) * 100) # use if you want to calculate the percent change from baselinea
  tbl <- mutate(tbl, ratio = tbl$y2091/tbl$y2011)
  tbl <- mutate(tbl, logRatio = log2(ratio))
  tbl <- mutate(tbl, gc = as.factor(gc))
  out <- glue('./qs_{studyAreaName}/occur')
  ifelse(!dir.exists(out), dir.create(out, recursive = TRUE), print('Folder already exist'))
  
  qs::qsave(x = tbl, file = glue('{out}/occur_changes_{studyAreaName}_{yr1}-{yr2}_{spc}.qs'))
  breaks <- seq(-1, 1, by = 0.25)
  message(crayon::green('Making map for: ',spc))
  ggRatio <- ggplot() +
    geom_tile(data = tbl, aes(x = lon, y = lat, color = change)) + # use logRatio or change
    #geom_sf(data = dehcho, fill = NA, col = 'gray10') + # !!! use only when studyAreaName = 'dehcho'
    geom_sf(data = edeh, fill= NA, col = 'gray10') +
 #   scale_fill_gradientn(colours = colorspace::diverging_hcl(n = 10, palette = 'Green-Brown'),rev = TRUE,
  #                       na.value = '#999999',breaks = breaks,limits = c(-1, 1) * max(abs(tbl$change))) +
   scale_colour_gradientn(colours = brewer.pal(n = 10, name = 'BrBG'), limit = c(-1,1))+#limit = c(-1, 1) * max(abs(tbl$change)))+
    facet_wrap(. ~ gc, ncol = 3, nrow = 1) +
    ggtitle(label = glue('{spc} ({yr1}-{yr2})')) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          legend.key.width = unit(4, 'line'),
          plot.title = element_text(size = 16, face = 'bold', hjust = 0, vjust = 0.7), 
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12), 
          legend.text = element_text(size = 10), 
          legend.title = element_text(size = 12, face = 'bold'), 
          strip.text = element_text(size =12)) +
    labs(x = 'Longitude', y = 'Latitude', fill = 'Change')
  out <- glue('./maps/occur/changes/{studyAreaName}') 
  ifelse(!dir.exists(out), dir.create(out, recursive = TRUE), 'Folder already exist')
  ggsave(plot = ggRatio,filename = glue('{out}/ocurr_change_{yr1}-{yr2}_{studyAreaName}_{spc}.png'),
         units = 'in', width = 12, height = 9, dpi = 700)
}

# Apply the function -----------------------------------------------------
purrr::map(.x= species, edeh,'edeh', '2011', '2031', .f = changes_rasters)

