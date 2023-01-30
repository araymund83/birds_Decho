library(pacman)

pacman::p_load(glue, qs, stringr, sf, tidyverse,fs)

rm(list = ls())

# Load data ---------------------------------------------------------------
root <- './outputs'
dirs <- fs::dir_ls(root, type = 'directory')
species <- basename(dirs)
species <- species[1:72]

#Create summary table summing all rows for each species ------------------
sum_table <- function(specie, studyAreaName, yr1, yr2){
  
 #specie <- species[1]
  message(crayon::blue('Starting\n', specie, '\n'))
  out <- glue('./qs_{studyAreaName}/changes/')
  table <- qs::qread(file = glue('{out}changes_{studyAreaName}_{yr1}-{yr2}_{specie}.qs'))
  sumDF <-table %>% group_by(gc) %>% summarise(across(y2011:y2100, sum))
  sumDF <- mutate(sumDF, species = specie)
  cat('Done \n')
  return(sumDF)
}  

sumTable <- map(.x = species,'edeh', '2011','2091', .f = sum_table)

totalTable <- bind_rows(sumTable)  
out <- glue('./qs_{studyAreaName}/changes/')
ifelse(!dir.exists(out), dir.create(out, recursive = TRUE), print('Folder already exist'))
qs::qsave(totalTable, glue('{out}/totalTable_{studyAreaName}_allsp.qs'))

cellArea <- 6.25 # comes from pixel resolution 250 * 250m it will be use to calculate the total abundance 

message(crayon::green(glue('Estimating change in abundance from {yr1} and {yr2}\n')))
abundTable <- totalTable %>% mutate(across((y2011:y2100), function(x) x * cellArea)) ##gets the total abundance

qs::qsave(x = abundTable, file = glue('{out}abundanceTable_{studyAreaName}_allsp.qs'))

tbl <- abundTable %>% mutate(change = y2031 - y2011,
                             ratio = (y2031/y2011),
                             logRatio = log2(ratio),
                             change2 = y2091 - y2011,
                             ratio2 = (y2091/y2011),
                             logRatio2 = log2(ratio2),
                             gc = as.factor(gc)) %>% 
  select(gc, species, logRatio, logRatio2) %>% rename('2011-2031' = 'logRatio',
                                                      '2011-2091' = 'logRatio2')

tbl<- tbl %>%  group_by(gc, species) %>% arrange(species) %>%
  pivot_wider(names_from = gc, values_from= c('2011-2031','2011-2091')) 
#save the table to open in excel
write.csv(tbl, glue('./outputs/logRatioTable_{studyAreaName}_allsp.csv'))




