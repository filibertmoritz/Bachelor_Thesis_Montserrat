###### Run MacKenzie-Bailaey Goodness-of-Fit Tests on dynamic occupancy models as a patch process ###########
###### Script written by Filibert Heim in June 2024, filibert.heim@posteo.de ############


# load packages 
library(AICcmodavg)
library(unmarked)
library(tidyverse)
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# set working directory 
getwd() 
setwd('C:/Users/filib/Documents/Studium/Bachelorarbeit/R_BachelorThesisMontserrat/')

# store species names 
SPECIES <- c('MTOR','FOTH','BRQD','TREM','ACHU','PTCA','PETH','GTCA','SBTH','SNPI','CAEL','BANA')

##### 1. Load best model data per species and store them in list #####

best_models <- vector(mode = 'list', length = length(SPECIES)) # create input vector with mode 'list'
for(s in 1:length(SPECIES)){
  best_models[[s]] <- readRDS(file = sprintf('output/data/best_models/%s_best_model.rds', SPECIES[s])) # read in file
  names(best_models)[s] <- SPECIES[s]} # name each part of the list with species name


##### 2. Run goodness of fit in a loop for all species together 

for(s in 1:length(SPECIES)){
  (best_model_gof <- mb.gof.test(best_models[[s]], nsim = 1000)) # specify number of iterations - min: 1000
  saveRDS(best_model_gof, file = sprintf('output/data/GOF/%s_gof_mb.rds', SPECIES[s]))
}

