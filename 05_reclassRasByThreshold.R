# Load libraries --------------------------------------------------------
require(pacman)

pacman::p_load(raster, rgdal, rgeos, terra, stringr, glue, sf, tidyverse, RStoolbox, fs, fst, trend)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data ---------------------------------------------------------------
thrs <- read_csv('./inputs/prevOcc625.csv')
root <- './outputs'
dirs <- fs::dir_ls(root, type = 'directory')

# Load data ---------------------------------------------------------------4
path <- './outputs/occur'
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

# See the changes  --------------------------------------------------------
reclass_ras <- function(spc, studyAreaName){
  fls <- list.files(glue('{path}/{studyAreaName}'), full.names = TRUE) 
  # Proof
 # spc <- species[37] # Run and comment (after)
  cat('Start ', spc, '\n')
  fls <- grep(spc,fls, value = TRUE)
  fls <- as.character(fls)
  yrs <- parse_number(basename(fls))
  yrs <- unique(yrs)
  gcm <- str_sub(basename(fls), start = 21, end = nchar(basename(fls)) - 4)
  gcm <- unique(gcm)
  thr <- filter(thrs, spec == spc)
  vle <- unique(thr$pOccMean)
  
  cat('Reclassifying\n')
  dfm <- map(.x = 1:length(gcm), .f = function(k){
    
    cat(k, '\n')
    fl <- grep(gcm[k], fls, value = TRUE)
    cat(gcm[k],'\n')
    rs <- terra::rast(fl)
    rs[rs < vle] <- 0
    names(rs)<- yrs
    out <- glue ('./outputs/{spc}/occuReclass/{studyAreaName}')
    ifelse(!dir.exists(out), dir.create(out, recursive = TRUE), print('Dir already exists'))
    writeRaster(rs, glue('{out}/occuReclass_{spc}_{gcm[k]}_{names(rs)}.tif'),
                overwrite = TRUE) # write each of the layers as an individual raster 
    df <- terra::as.data.frame(x = rs, xy = TRUE, na.rm = TRUE)
    colnames(df) <- c('x', 'y', yrs)
    
    df <- as_tibble(df)
    df <- mutate(df, gc = gcm[k])
    return(df)
    
  })
  
  rsl <- bind_rows(dfm)
  out <- glue ('./qs_{studyAreaName}/occuReclass')
  ifelse(!dir.exists(out),dir.create(out, recursive = TRUE), print('Dir already exists'))
  qs::qsave(x = rsl, file = glue('{out}/occmsk_reclass_yrs_{spc}.qs'))
  
  cat('------- Done -------\n')
  return(rsl)
} 
# Apply the function ------------------------------------------------------
rslt <- map(.x = species, 'edeh', .f = reclass_ras)



# Now raster to table -----------------------------------------------------
path <- glue()
fles <- glue('./outputs/{spcs}/occur/occmsk_yrs_{spcs}_625.qs')
fles <- as.character(fles)
head(qs::qread(fles[1]))
test <- qs::qread(fles[1])
dout <- './graphs/figs/occur'

# Function table to raster ------------------------------------------------
tbl2rst <- function(fle){
  
  #fle <- fles[1] # Run and erase (after)
  cat('Start ', basename(fle), '\n')
  tbl <- qs::qread(fle)
  spc <- basename(fle) %>% str_split(., '_') %>% sapply(., `[[`, 3) %>% gsub('.qs', '', .) 
  gcm <- unique(tbl$gc)
  rsl <- map(gcm, function(i){
    cat(i, '\n')
    tb <- filter(tbl, gc == i)
    colnames(tb) <- c('x', 'y', 'y2011', 'y2031', 'y2051', 'y2071', 'y2091', 'y2100', 'gc')
    tb <- mutate(tb, change = y2091 - y2011) #ratio = (y2091/y2011), logRatio = log2(ratio))
    tb <- dplyr::select(tb, -gc)
    tr <- terra::rast(tb[1:9], type = 'xyz')
    do <- glue('{dout}/{spc}')
    ifelse(!dir.exists(do), dir_create(do), print('Directory already exists'))
    terra::writeRaster(x = tr, filename = glue('{do}/occ_change_msk_{spc}_{i}_625.tif'), overwrite = TRUE)
    cat('Done!\n')
    return(tb[,c(1, 2, 9)])
  })
  rsl <- map(.x = 1:length(gcm), .f = function(i){
    mutate(rsl[[i]], gc = gcm[i])
  })
  rsl <- bind_rows(rsl)
  #rsl <- mutate(rsl, class = ifelse(change < 0, 'Loss', ifelse(change == 0, 'No change', 'Gain')))
  #rsl <- mutate(rsl, class = factor(class, levels = c('Loss', 'No change', 'Gain')))
  pal <- wesanderson::wes_palette('Zissou1', 100, type = 'continuous')
  colours = c("#FFFFFF", "#9999FF", "#66FFFF", "#66FFCC", "#99FF33", "#FFFF33", 
              "#FFCC33", "#FF9900", "#FF6600", "#FF0000", "#CC0000")
  
  cat('To make the map\n')
  ggp <- ggplot() +
    geom_tile(data = rsl, aes(x = x, y = y, fill = change)) + 
    geom_sf(data = limt, fill = NA, col = '#c8c8c8') +
    geom_sf(data = ecrg, fill = NA, col = '#D3D3D3') +
    scale_fill_distiller(palette = 'RdBu',
                         type = 'div',
                         limits = c(-1, 1) * max(abs(rsl$change)),
                         #limits = c(min(rsl$change), max(rsl$change))* max(abs(rsl$change)),
                         direction = + 1,
                         guide= 'colourbar') +
    facet_wrap(.~gc, ncol = 3, nrow = 1) + 
    ggtitle(label = glue('{spc} (2011-2091)')) +
    labs(x = 'Longitude', y = 'Latitude') +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.text.y = element_text(angle = 90, hjust = 0.5), 
          plot.title = element_text(size = 14, face = 'bold'), 
          legend.position = 'bottom', 
          legend.key.width = unit(1.5, 'cm')) +
    coord_sf()
  
  
  dou <- glue('./graphs/figs/occur/{spc}/occ_change2_{spc}_11_91_625.png')
  ggsave(plot = ggp, filename = dou, units = 'in', width = 12, 
         height = 9, dpi = 700)
  cat('------ Done ------!\n')
}

# Apply the function ------------------------------------------------------
rslt <- map(.x = fles, .f = tbl2rst)
