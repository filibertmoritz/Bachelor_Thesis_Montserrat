##### making the data more pretty for presentation in tables and figures etc. #####
##### written by Filibert Heim in Feb 2024, filibert.heim@posteo.de ###############

# load required packages
library(tidyverse)
library(AICcmodavg)
library(stringr)
#install.packages("patchwork")
library(patchwork)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(grid)
#install.packages('stargazer')
#install.packages('tsutils')
library(stargazer)
library(tsutils)
library(gt) # for tables
library(webshot2) # to save tables
library('unmarked')
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# set working directory 
setwd('C:/Users/filib/Documents/Studium/Bachelorarbeit/R_BachelorThesisMontserrat/output')

# load species names (abbreviation and full names) and arrange them by alphabetical order
SPECIES <- c('MTOR','FOTH','BRQD','TREM','ACHU','PTCA','PETH','GTCA','SBTH','SNPI','CAEL','BANA')
Species_names <- c('Montserrat Oriole','Forest Thrush','Bridled Quail-Dove','Brown Trembler','Antillean Crested Hummingbird','Purple-Throated Carib','Pearly-Eyed Thrasher','Green-Throated Carib','Scaly-Breasted Thrasher','Scaly-Naped Pigeon','Caribbean Elaenia','Bananaquit')
data <- data.frame(names = Species_names, codes = SPECIES) %>% arrange(SPECIES)
Species_names_sorted <- data$names
SPECIES_sorted <- data$codes
rm('data')

#### 1. access goodness of fit tests from RDS #####

# read GOF Rdata in a list 
GOF_mb <- vector(mode = 'list', length = length(SPECIES)) # create input vector with mode 'list'
for(s in 1:length(SPECIES)){
  GOF_mb[[s]] <- readRDS(file = sprintf('data/GOF/%s_gof_mb.rds', SPECIES[s])) # read in file
  names(GOF_mb)[s] <- SPECIES[s]} # name each part of the list with species name

# access GOF data for summary table
GOF_df <- data.frame(Species = NA, c_hat = NA, p_value = NA)
for(z in 1:length(SPECIES)){
 GOF_df[z,] <- c(SPECIES[z], round(GOF_mb[[z]]$c.hat.est, digits = 2), round(GOF_mb[[z]]$p.global, digits = 3))}
GOF_df %>% arrange(Species)

#### 2. full modSel-tables ####

# read data in a list 
modsel_full <- vector(mode = 'list', length = length(SPECIES)) # create input vector with mode 'list'
for(s in 1:length(SPECIES)){
  modsel_full[[s]] <- data.frame(read.csv(file = sprintf('data/modSel/%s_modSel_full.csv', SPECIES[s]), header = T, sep = ','))
  names(modsel_full)[s] <- SPECIES[s]}

# that is another way to store all the files as objects in global environment
#for(s in 1:length(SPECIES)){
#  assign(sprintf('modSel_full_%s', SPECIES[s]),
#         read.csv2(file = sprintf('data/modSel/%s_modSel_full.csv', SPECIES[s]), header = T, sep = ','))
#} 

# create data.frame containing all species modSel tables
modsel_full_species <- data.frame(modsel_full[[1]]) # initial data.frame with first species 
for(s in 2:length(SPECIES)){
  modsel_full_species  <- rbind(modsel_full_species, modsel_full[[s]])
  modsel_full_species  <- na.omit(modsel_full_species)}

# make some adjustments (rename columns, round values to 4 digits, safe values as.character)
modsel_full_species <- modsel_full_species %>% 
  rename(Species = species, Model = model, Formula = formula, Step = step, 
         'ΔAIC' = delta, 'AIC weight' = AICwt, 'Cumulative Akaike weights' = cumltvWt)
modsel_full_species[,7:10] <- round(modsel_full_species[,7:10], digits = 4)
modsel_full_species <- modsel_full_species %>% 
  mutate(`AIC weight` = as.character(`AIC weight`), `ΔAIC` = as.character(`ΔAIC`), 
         `Cumulative Akaike weights` = as.character(`Cumulative Akaike weights`))

# count overall number of fitted models!
dim(modsel_full_species %>% group_by(Species) %>% count(Formula))[1]
modsel_full_species %>% group_by(Species) %>% count(Formula) %>% filter(Species == 'GTCA') # number of models per species 
modsel_full_species %>% group_by(Species) %>% count(Formula) %>% filter(Species == 'MTOR') 
modsel_full_species %>% ungroup()

# derive different tables from this parent table and Formula# derive different tables from this parent table and export them

###### a) modSel full species (for Appendix full species,full parameters,full steps) #####

# derive table from global modsel table by selecting needed columns and renaming columns
modsel_full_species_export <- modsel_full_species # safe new object to prevent changes in global table

# add GOF c-hat and p-value to table
minAIC <- modsel_full_species %>% filter(Step == 'g_e') %>% aggregate(AIC ~ Species, FUN = min) %>% mutate(Step = 'g_e') # finds minimal AIC for each species 
modsel_GOF <- merge(GOF_df, minAIC, by = c('Species')) # add GOF data
modsel_full_species_export <- full_join(x = modsel_full_species_export, y = modsel_GOF, by = c('Species', 'Step', 'AIC')) # join both data sets

# make the table prettier
modsel_full_species_export <- modsel_full_species_export %>%
  select(-X, -Model, -`Cumulative Akaike weights`) %>% 
  arrange(Species,  match(Step, c('g_e', 'psi', 'p')), AIC) %>% 
  mutate_all(~ifelse(is.na(.), "-", .)) %>% 
  mutate(p_value = if_else(p_value == '0', false = as.character(p_value), true = '>0.001')) %>% 
  select(Species, Step, everything()) %>%
  rename(`Formula (~Ψ~γ~ε~P)` = Formula, `p-value` = p_value, `c-hat` = c_hat)  
modsel_full_species_export$Step[modsel_full_species_export$Step == 'p'] <- 'P'
modsel_full_species_export$Step[modsel_full_species_export$Step == 'psi'] <- 'Ψ'
modsel_full_species_export$Step[modsel_full_species_export$Step == 'g_e'] <- 'ε,γ'

modsel_full_species_export %>% 
  mutate(RowNumber = row_number()) %>% 
  group_by(Species) %>% 
  filter(AIC == min(AIC)) %>% 
  select(RowNumber)

# make pretty table for modselFullSpecies using gt
gt_modsel_full_species_export = gt(modsel_full_species_export)
gt_modsel_full_species_export  <- gt_modsel_full_species_export   %>%
  tab_header(
    title = md("Model Comparison Tables for All Species and All Steps")) %>%
  tab_spanner(
    label = "GoF",
    columns = c("c-hat", "p-value")) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")),
    locations = list(cells_column_labels(
      columns = everything()), cells_title())) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")),
    locations = list(cells_body(rows = c(1,24,47,70,96,119,145,168,191,214,237,260)))) %>%
  cols_align(
    align = "right",
    columns = c("Formula (~Ψ~γ~ε~P)", "AIC", "ΔAIC", "AIC weight", "nPars", "p-value", "c-hat")) %>%
  cols_align(
    align = "left", 
    columns = c("Species")) %>% 
  cols_width(
    'Formula (~Ψ~γ~ε~P)' ~ px(595)) %>%
  tab_source_note( # add notes below the table coded as markdown by md()
    source_note = md("**Models**: Refers to the candidate models in the last modelling step (ε, γ) (more detailed explaination in the methods section).<br>
    **Model Parameters**: Ψ = Initial Occupancy; γ = Colonization Probability;  ε = Extinction Probability; P = Detection Probability.<br>
    **Model Predictors**: day = survey day; time = survey time; rain = 1 - 'rain' and 0 - 'no rain' during survey; wind = wind strength ('calm' - 1, 'light' - 2, 'moderate' - 3, 'high' - 4) during survey; activity = general bird activity; altitude = altitude in m a.s.l., ridge = 1 - 'ridge' or 0 - 'valley'/'midslope'; dbh = diameter-breast height (DBH); treeheight = tree-height around the survey point; year = survey year (see Table 1 for details of measurements).<br>
    **Assessment of Model Fit**: GoF = MacKenzie-Bailey Goodness-of-Fit Test; c-hat = Overdispersion Parameter.<br>
    **Model Comparison Parameter**: nPars = Number of Parameters; AIC = Aikake Information Criterion; ΔAIC = Delta AIC; AIC weight = Akaike weights.<br>
    **Species**: ACHU = Antillean Crested Hummingbird; BANA = Bananaquit; BRQD = Bridled Quail-Dove; 
    CAEL = Caribbean Elaenia; FOTH = Forest Thrush; GTCA = Green-Throated Carib; MTOR = Montserrat Oriole; 
    PETH = Pearly-Eyed Thrasher; PTCA = Purple-Throated Carib; SBTH = Scaly-Breasted Thrasher; SNPI = Scaly-Naped Pigeon; 
    TREM = Brown Trembler.")) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('top'), color = 'black', weight = px(2.25))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('bottom'), color = 'black', weight = px(1.5))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c("top"), color = 'gray', weight = px(1))),
    locations = list(cells_column_labels(), cells_body(rows = c(1,24,47,70,96,119,145,168,191,214,237,260)))) %>%
  tab_style(
    style = list(
      cell_borders(sides = c("bottom"), color = 'black', weight = px(2.25))),
    locations = list(cells_body(rows = dim(modsel_full_species_export)[1]))) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('bottom'), color = 'black', weight = px(1.5))),
    locations = list(cells_column_labels()))
gt_modsel_full_species_export

#save the table: create an html file first and store it using gtsave
html_file = "MODSEL_FULL.html"
gtsave(gt_modsel_full_species_export, html_file, path = paste0(getwd(), "/plot/tables/")) 


html_file = "MODSEL_FULL.html"
gtsave(gt_modsel_full_species_export, html_file, path = paste0(getwd(), "/plot/tables/")) 
webshot2::webshot(paste0(getwd(), "/plot/tables/", html_file), 
                  file = paste0(getwd(), "/plot/tables/", "MODSEL_FULL.png"), cliprect = 'viewport', delay = 5,
                  vwidth = 1200,
                  vheight = 11160, zoom = 2.5)

  
###### b) modsel MTOR all steps (as example workflow) #####

modsel_MTOR_all_steps <- modsel_full_species %>% 
  filter(Species == 'MTOR') %>% 
  select(-X, -Model, -Species) %>% 
  arrange(match(Step, c('p', 'psi', 'g_e')), Formula)

###### c) modsel full species only step g_e (col and ext - final candidate models) ####

# derive table from global modsel table by selecting needed columns and renaming columns
modsel_full_species_export_g_e <- modsel_full_species %>% filter(Step == 'g_e') # safe new object to prevent changes in global table

# add GOF c-hat and p-value to table
minAIC <- modsel_full_species %>% filter(Step == 'g_e') %>% aggregate(AIC ~ Species, FUN = min) # finds minimal AIC for each species 
modsel_GOF <- merge(GOF_df, minAIC, by = c('Species')) # add GOF data
modsel_full_species_export_g_e <- full_join(x = modsel_full_species_export_g_e, y = modsel_GOF, by = c('Species', 'AIC')) # join both data sets

# make the table prettier
modsel_full_species_export_g_e <- modsel_full_species_export_g_e %>%
  select(-X, -`Cumulative Akaike weights`) %>% 
  mutate(`AIC weight` = as.character(format(as.numeric(`AIC weight`), scientific = F))) %>%
  mutate_all(~ifelse(is.na(.), "-", .)) %>% 
  mutate(p_value = if_else(p_value == '0', false = as.character(p_value), true = '>0.001'), `AIC weight` = if_else(`AIC weight` == '0.0000', true = '0', false = `AIC weight`),  Model = str_to_title(Model), Model = str_replace(Model, pattern = '_', replace = '-')) %>% 
  select(-Step) %>%
  rename(`Formula (~Ψ~γ~ε~P)` = Formula, `p-value` = p_value, `c-hat` = c_hat)  
modsel_full_species_export$Step[modsel_full_species_export$Step == 'p'] <- 'P'
modsel_full_species_export$Step[modsel_full_species_export$Step == 'psi'] <- 'Ψ'
modsel_full_species_export$Step[modsel_full_species_export$Step == 'g_e'] <- 'ε,γ'

index <- modsel_full_species_export_g_e %>% 
  mutate(RowNumber = row_number()) %>% 
  group_by(Species) %>% 
  filter(AIC == min(AIC)) %>% 
  select(RowNumber)

# make pretty table for modselFullSpecies
gt_modsel_full_species_g_e = gt(modsel_full_species_export_g_e)
gt_modsel_full_species_g_e  <- gt_modsel_full_species_g_e   %>%
  tab_header(
    title = md("Model Comparison Table for Final Candidate Models for Each Species")) %>%
  tab_spanner(
    label = "GoF",
    columns = c("c-hat", "p-value")) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")),
    locations = list(cells_column_labels(
      columns = everything()), cells_title())) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")),
    locations = list(cells_body(rows = index$RowNumber))) %>%
  cols_align(
    align = "right",
    columns = c("Formula (~Ψ~γ~ε~P)", "AIC", "ΔAIC", "AIC weight", "nPars", "p-value", "c-hat")) %>%
  cols_align(
    align = "left", 
    columns = c("Model", "Species")) %>% 
  cols_width(
    'Formula (~Ψ~γ~ε~P)' ~ px(595)) %>%
  tab_source_note(
    source_note = md("**Models**: Refers to the candidate models in the last modelling step (ε, γ) (more detailed explaination in the methods section).<br>
    **Model Parameters**: Ψ = Initial Occupancy; γ = Colonization Probability;  ε = Extinction Probability; P = Detection Probability.<br>
    **Model Predictors**: day = survey day; time = survey time; rain = 1 - 'rain' and 0 - 'no rain' during survey; wind = wind strength ('calm' - 1, 'light' - 2, 'moderate' - 3, 'high' - 4) during survey; activity = general bird activity; altitude = altitude in m a.s.l., ridge = 1 - 'ridge' or 0 - 'valley'/'midslope'; dbh = diameter-breast height (DBH); treeheight = tree-height around the survey point; year = survey year (see Table 1 for details of measurements).<br>
    **Assessment of Model Fit**: GoF = MacKenzie-Bailey Goodness-of-Fit Test; c-hat = Overdispersion Parameter.<br>
    **Model Comparison Parameter**: nPars = Number of Parameters; AIC = Aikake Information Criterion; ΔAIC = Delta AIC; AIC weight = Akaike weights.<br>
    **Species**: ACHU = Antillean Crested Hummingbird; BANA = Bananaquit; BRQD = Bridled Quail-Dove; 
    CAEL = Caribbean Elaenia; FOTH = Forest Thrush; GTCA = Green-Throated Carib; MTOR = Montserrat Oriole; 
    PETH = Pearly-Eyed Thrasher; PTCA = Purple-Throated Carib; SBTH = Scaly-Breasted Thrasher; SNPI = Scaly-Naped Pigeon; 
    TREM = Brown Trembler.")) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('top'), color = 'black', weight = px(2.25))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('bottom'), color = 'black', weight = px(1.5))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c("top"), color = 'gray', weight = px(1))),
    locations = list(cells_column_labels(), cells_body(rows = index$RowNumber))) %>%
  tab_style(
    style = list(
      cell_borders(sides = c("bottom"), color = 'black', weight = px(2.25))),
    locations = list(cells_body(rows = dim(modsel_full_species_export_g_e)[1]))) %>% 
  tab_style(
      style = list(
        cell_borders(sides = c('bottom'), color = 'black', weight = px(1.5))),
      locations = list(cells_column_labels()))
gt_modsel_full_species_g_e

#save the table (a bit complicated, but ok): create an html file first, that is opened by webshot2 and then takes a photo and stores it as .png
html_file = "MODSEL_FULL_g_e.html"
gtsave(gt_modsel_full_species_g_e, html_file, path = paste0(getwd(), "/plot/tables/")) 
webshot2::webshot(paste0(getwd(), "/plot/tables/", html_file), 
                  file = paste0(getwd(), "/plot/tables/", "MODSEL_FULL_g_e.png"), 
                  vwidth = 1300,
                  vheight = 5200, zoom = 3)

###### e) modSel table for Bachelor Thesis in two parts for better layout

#### part 1 - select first 6 analysed species
gt_modsel_full_species_g_e_1 = gt(modsel_full_species_export_g_e %>% 
                                  filter(Species %in% c(unique(modsel_full_species_export_g_e$Species)[1:6])))

index <- modsel_full_species_export_g_e %>% # update index for proper selection of rows
  filter(Species %in% c(unique(modsel_full_species_export_g_e$Species)[1:6])) %>%
  mutate(RowNumber = row_number()) %>% 
  group_by(Species) %>% 
  filter(AIC == min(AIC))

gt_modsel_full_species_g_e_1  <- gt_modsel_full_species_g_e_1 %>%
  tab_header(title = md("Model Comparison Table for Final Candidate Models for Each Species, Part a)")) %>%
  tab_spanner(
    label = "GoF",
    columns = c("c-hat", "p-value")) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")),
    locations = list(cells_column_labels(
      columns = everything()), cells_title())) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")),
    locations = list(cells_body(rows = index$RowNumber))) %>%
  cols_align(
    align = "right",
    columns = c("Formula (~Ψ~γ~ε~P)", "AIC", "ΔAIC", "AIC weight", "nPars", "p-value", "c-hat")) %>%
  cols_align(
    align = "left", 
    columns = c("Model", "Species")) %>% 
  cols_width(
    'Formula (~Ψ~γ~ε~P)' ~ px(595)) %>%
  tab_footnote(footnote = 'Information on all analysed but here not displayed species PETH, GTCA, SBTH, SNPI, CAEL and BANA are displayed in part b) of the same table', locations = cells_title(groups = 'title')) %>%
  tab_source_note(
    source_note = md("**Models**: Refers to the candidate models in the last modelling step (ε, γ) (more detailed explaination in the methods section).<br>
    **Model Parameters**: Ψ = Initial Occupancy; γ = Colonization Probability;  ε = Extinction Probability; P = Detection Probability.<br>
    **Model Predictors**: day = survey day; time = survey time; rain = 1 - 'rain' and 0 - 'no rain' during survey; wind = wind strength ('calm' - 1, 'light' - 2, 'moderate' - 3, 'high' - 4) during survey; activity = general bird activity; altitude = altitude in m a.s.l., ridge = 1 - 'ridge' or 0 - 'valley'/'midslope'; dbh = diameter-breast height (DBH); treeheight = tree-height around the survey point; year = survey year (see Table 1 for details of measurements).<br>
    **Assessment of Model Fit**: GoF = MacKenzie-Bailey Goodness-of-Fit Test; c-hat = Overdispersion Parameter.<br>
    **Model Comparison Parameter**: nPars = Number of Parameters; AIC = Aikake Information Criterion; ΔAIC = Delta AIC; AIC weight = Akaike weights.<br>
    **Species**: ACHU = Antillean Crested Hummingbird; BRQD = Bridled Quail-Dove; FOTH = Forest Thrush; MTOR = Montserrat Oriole; PTCA = Purple-Throated Carib; 
    TREM = Brown Trembler.")) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('top'), color = 'black', weight = px(2.25))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('bottom'), color = 'black', weight = px(1.5))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c("top"), color = 'gray', weight = px(1))),
    locations = list(cells_column_labels(), cells_body(rows = index$RowNumber))) %>%
  tab_style(
    style = list(
      cell_borders(sides = c("bottom"), color = 'black', weight = px(2.25))),
    locations = list(cells_body(rows = 39))) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('bottom'), color = 'black', weight = px(1.5))),
    locations = list(cells_column_labels()))
gt_modsel_full_species_g_e_1

#save the table (a bit complicated, but ok): create an html file first, that is opened by webshot2 and then takes a photo and stores it as .png
html_file = "MODSEL_FULL_g_e.html"
gtsave(gt_modsel_full_species_g_e_1, html_file, path = paste0(getwd(), "/plot/tables/")) 
webshot2::webshot(paste0(getwd(), "/plot/tables/", html_file), 
                  file = paste0(getwd(), "/plot/tables/", "MODSEL_FULL_g_e_part_a.png"), 
                  vwidth = 1300,
                  vheight = 5200, zoom = 3)


#### part 2 - select last 6 analysed species
gt_modsel_full_species_g_e_2 = gt(modsel_full_species_export_g_e %>% 
                                    filter(Species %in% c(unique(modsel_full_species_export_g_e$Species)[7:12])))

index <- modsel_full_species_export_g_e %>% # update index for proper selection of rows
  filter(Species %in% c(unique(modsel_full_species_export_g_e$Species)[7:12])) %>%
  mutate(RowNumber = row_number()) %>% 
  group_by(Species) %>% 
  filter(AIC == min(AIC))

gt_modsel_full_species_g_e_2  <- gt_modsel_full_species_g_e_2 %>%
  tab_header(title = md("Model Comparison Table for Final Candidate Models for Each Species, Part b)")) %>%
  tab_spanner(
    label = "GoF",
    columns = c("c-hat", "p-value")) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")),
    locations = list(cells_column_labels(
      columns = everything()), cells_title())) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")),
    locations = list(cells_body(rows = index$RowNumber))) %>%
  cols_align(
    align = "right",
    columns = c("Formula (~Ψ~γ~ε~P)", "AIC", "ΔAIC", "AIC weight", "nPars", "p-value", "c-hat")) %>%
  cols_align(
    align = "left", 
    columns = c("Model", "Species")) %>% 
  cols_width(
    'Formula (~Ψ~γ~ε~P)' ~ px(595)) %>%
  tab_footnote(footnote = 'Information on all analysed but here not displayed species ACHU, FOTH, MTOR, PTCA and TREM are displayed in part a) of the same table', locations = cells_title(groups = 'title')) %>%
  tab_source_note(
    source_note = md("**Models**: Refers to the candidate models in the last modelling step (ε, γ) (more detailed explaination in the methods section).<br>
    **Model Parameters**: Ψ = Initial Occupancy; γ = Colonization Probability;  ε = Extinction Probability; P = Detection Probability.<br>
    **Model Predictors**: day = survey day; time = survey time; rain = 1 - 'rain' and 0 - 'no rain' during survey; wind = wind strength ('calm' - 1, 'light' - 2, 'moderate' - 3, 'high' - 4) during survey; activity = general bird activity; altitude = altitude in m a.s.l., ridge = 1 - 'ridge' or 0 - 'valley'/'midslope'; dbh = diameter-breast height (DBH); treeheight = tree-height around the survey point; year = survey year (see Table 1 for details of measurements).<br>
    **Assessment of Model Fit**: GoF = MacKenzie-Bailey Goodness-of-Fit Test; c-hat = Overdispersion Parameter.<br>
    **Model Comparison Parameter**: nPars = Number of Parameters; AIC = Aikake Information Criterion; ΔAIC = Delta AIC; AIC weight = Akaike weights.<br>
    **Species**: BANA = Bananaquit;
    CAEL = Caribbean Elaenia; GTCA = Green-Throated Carib;
    PETH = Pearly-Eyed Thrasher; SBTH = Scaly-Breasted Thrasher; SNPI = Scaly-Naped Pigeon.")) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('top'), color = 'black', weight = px(2.25))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('bottom'), color = 'black', weight = px(1.5))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c("top"), color = 'gray', weight = px(1))),
    locations = list(cells_column_labels(), cells_body(rows = index$RowNumber))) %>%
  tab_style(
    style = list(
      cell_borders(sides = c("bottom"), color = 'black', weight = px(2.25))),
    locations = list(cells_body(rows = 42))) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('bottom'), color = 'black', weight = px(1.5))),
    locations = list(cells_column_labels()))
gt_modsel_full_species_g_e_2

#save the table (a bit complicated, but ok): create an html file first, that is opened by webshot2 and then takes a photo and stores it as .png
html_file = "MODSEL_FULL_g_e.html"
gtsave(gt_modsel_full_species_g_e_2, html_file, path = paste0(getwd(), "/plot/tables/")) 
webshot2::webshot(paste0(getwd(), "/plot/tables/", html_file), 
                  file = paste0(getwd(), "/plot/tables/", "MODSEL_FULL_g_e_part_b.png"), 
                  vwidth = 1300,
                  vheight = 5200, zoom = 3)

###### e) best model (bm) per Species model with occupancy and effect sizes on original scale !!! ######

# load data 
bm_original <- read.csv2(file = 'Montserrat_ResultsTable_original_scale.csv', dec = '.', sep = ';', header = T)
names(bm_original)[2:length(names(bm_original))] <- c('Formula (~Ψ~γ~ε~P)', 'Initial Occupancy Ψ (CI)', 'Ψ (CI)', 'γ (CI)', 'ε (CI)', 'c-hat', 'p-value') # rename columns

model_list <- modsel_full_species_export_g_e %>% group_by(Species) %>% slice(which.min(AIC)) # get the correct of best models order to add them to table

# change order and add best model description for better overview
bm_original <- bm_original %>%
  arrange(Species) %>% 
  mutate(Model = str_to_title(model_list$Model)) %>% 
  select(Species, Model, everything()) %>% 
  mutate(`p-value` = if_else(`p-value` == '0', false = as.character(`p-value`), true = '>0.001'))

# make pretty table
gt_bm_original = gt(bm_original)
gt_bm_original <- 
  gt_bm_original  %>%
  tab_header(
    title = md("The Effect of Elevation on Initial Occupancy Ψ, Colonization γ and Extinction ε Probabilities")) %>%
  tab_spanner(
    label = "GoF",
    columns = c("c-hat", "p-value")) %>%
  tab_spanner(
    label = "Effect of Elevation on",
    columns = c("Ψ (CI)","γ (CI)","ε (CI)")) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")),
    locations = list(cells_column_labels(
      columns = everything()), cells_title())) %>%
  cols_align(
    align = "right",
    columns = c("Formula (~Ψ~γ~ε~P)", "p-value")) %>%
  cols_align(
    align = "left", 
    columns = c("Model", "Species")) %>%
  cols_align(
    align = "center",
    columns = c("Initial Occupancy Ψ (CI)",                          
                "Ψ (CI)",       
                "γ (CI)",
                "ε (CI)",  
                "c-hat")) %>% 
  cols_width("Ψ (CI)" ~ px(200),
             "γ (CI)" ~ px(200), 
             "ε (CI)" ~ px(200), 
             "c-hat" ~ px(100), 
             "p-value" ~ px(100), 
             "Initial Occupancy Ψ (CI)" ~ px(200)) %>% 
  tab_footnote(
    footnote = "For FOTH and SBTH the effect sizes of alt:treeheight are displayed.",
    locations = cells_body(
      columns = c(5),
      rows = c(7,10))) %>% 
  tab_footnote(
    footnote = "CI values were incalculable from the best model and the error occured in sqrt(diag(vcov(obj, fixedOnly = fixedOnly))) : NaNs were produced.",
    locations = cells_body(
      columns = c(7),
      rows = c(5))) %>% 
  tab_source_note(
    source_note = md("**Model Parameters**: Ψ = Initial Occupancy; γ = Colonization Probability;  ε = Extinction Probability; P = Detection Probability; CI = Confidence Interval Limits (lower, upper).<br>
    **Model Predictors**: day = survey day; time = survey time; rain = 1 - 'rain' and 0 - 'no rain' during survey; wind = wind strength 'low' or 'high' during survey; activity = general bird activity; altitude = altitude in m a.s.l.; ridge = 1 - 'ridge' or 0 - 'valley'/'midslope'; dbh = diameter-breast height (DBH); treeheight = tree-height around the survey point; year = survey year (see Table 1 for details of measurements).<br>
    **Assessment of Model Fit**: GoF = MacKenzie-Bailey Goodness-of-Fit Test; c-hat = Overdispersion Parameter.<br>
    **Species**: ACHU = Antillean Crested Hummingbird; BANA = Bananaquit; BRQD = Bridled Quail-Dove; 
    CAEL = Caribbean Elaenia; FOTH = Forest Thrush; GTCA = Green-Throated Carib; MTOR = Montserrat Oriole; 
    PETH = Pearly-Eyed Thrasher; PTCA = Purple-Throated Carib; SBTH = Scaly-Breasted Thrasher; SNPI = Scaly-Naped Pigeon; 
    TREM = Brown Trembler.")) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('top'), color = 'black', weight = px(2.25))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('bottom'), color = 'black', weight = px(1.5))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c("bottom"), color = 'black', weight = px(2.25))),
    locations = list(cells_column_labels(), cells_body(rows = 12)))
gt_bm_original


# save the table (a bit complicated, but ok): create an html file first, that is opened by webshot2 and then takes a photo and stores it as .png
html_file = "COL_EXT_TABLE_original.html"
gtsave(gt_bm_original, html_file, path = paste0(getwd(), "/plot/tables/")) 
webshot2::webshot(paste0(getwd(), "/plot/tables/", html_file), 
                  file = paste0(getwd(), "/plot/tables/", "COL_EXT_TABLE_original.png"), 
                  vwidth = 1500,
                  vheight = 1200, zoom = 3)
webshot2::webshot(paste0(getwd(), "/plot/tables/", html_file), 
                  file = paste0(getwd(), "/plot/tables/", "COL_EXT_TABLE_original_small.png"), 
                  vwidth = 1300,
                  vheight = 1400, zoom = 3)

####### f) table as overview with predictors #####

# load data from csv table 
variables <- read.csv2(file = 'Variables_table.csv', header = T, sep = ';')

# make pretty table using gt 
gt_variables = gt(variables)
gt_variables <- 
  gt_variables  %>%
  tab_header(
    title = md("Definitions of Variables used as Predictors")) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")),
    locations = list(cells_column_labels(
      columns = everything()), cells_title())) %>%
  cols_align(
    align = "left", 
    columns = c("Type", "Covariate", "Metric", "Definition")) %>%
  tab_source_note(
    source_note = md("**Type**: obsCov = Observation Level Covariates, siteCov = Site Level Covariates,  yearlySiteCov = Year Level Covariates.<br>
                     **Transformation**: All numeric variables were scaled to a mean of 0 and standard deviation of 1.")) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('top'), color = 'black', weight = px(2.25))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('bottom'), color = 'black', weight = px(1.5))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c("bottom"), color = 'black', weight = px(2.25))),
    locations = list(cells_column_labels(), cells_body(rows = 10)))
gt_variables

# save the table (a bit complicated, but ok): create an html file first, that is opened by webshot2 and then takes a photo and stores it as .png
html_file = "variables.html"
gtsave(gt_variables, html_file, path = paste0(getwd(), "/plot/tables/")) 
webshot2::webshot(paste0(getwd(), "/plot/tables/", html_file), 
                  file = paste0(getwd(), "/plot/tables/", "variables.png"), 
                  vwidth = 700,
                  vheight = 1200, zoom = 3)
webshot2::webshot(paste0(getwd(), "/plot/tables/", html_file), 
                  file = paste0(getwd(), "/plot/tables/", "variables_wide.png"), 
                  vwidth = 1300,
                  vheight = 1200, zoom = 3)

####### h) yearly occupancy table #####



####### i) produce LATEX code for tables #####

# produce LATEX code for my tables 
stargazer(modsel_full_species, type = 'latex', title = 'Model comparison table of all investigated species', digits = 4, out = 'plot/table.html')


#### 3. produce table just with best models per species ####

# use full species table and filter for AIC best model 
modsel_best_models <- modsel_full_species %>% 
  filter(Step == 'g_e') %>% 
  group_by(Species) %>% 
  slice(which.min(AIC)) %>% 
  select(-X, -Step, -`Delta AIC`)

#### 4. merge graphics of occupancy for all species ####

# load .RDS data containing occupancy plots // pretty solution didn't work out:(
plotnames_occu <- list.files(path = 'plot/data/', pattern = 'occupancy_plot.rds')
occu_plots <- list()
for(i in 1:length(plotnames_occu)){
  occu_plots[[i]] <- readRDS(file = paste0('plot/data/', plotnames_occu[i]))}
names(occu_plots) <- str_remove(plotnames_occu, pattern = '_occupancy_plot.rds')

# design legend

# arrange graphs in one large graph: aligned
occu_plots[[13]] <- ggarrange(occu_plots$ACHU + rremove('xylab') + rremove('x.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)), 
          occu_plots$BANA + rremove('xylab') + rremove('xy.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16)), 
          occu_plots$BRQD + rremove('xylab') + rremove('xy.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16)), 
          occu_plots$CAEL + rremove('xylab') + rremove('x.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)), 
          occu_plots$FOTH + rremove('xylab') + rremove('xy.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16)),
          occu_plots$GTCA + rremove('xylab') + rremove('xy.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16)),
          occu_plots$MTOR + rremove('xylab') + rremove('x.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)),
          occu_plots$PETH + rremove('xylab') + rremove('xy.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16)),
          occu_plots$PTCA + rremove('xylab') + rremove('xy.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16)),
          occu_plots$SBTH + rremove('xylab') + theme(plot.margin = unit(c(0,0.2,0.2,1), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)),
          occu_plots$SNPI + rremove('xylab') + rremove('y.text') + theme(plot.margin = unit(c(0,0.2,0.2,1), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)),
          occu_plots$TREM + rremove('xylab') + rremove('y.text') + theme(plot.margin = unit(c(0,0.2,0.2,1), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)), 
          ncol = 3, nrow = 4, align = 'hv', labels = c('a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)'), label.x = .05) + 
  plot_annotation(title = "Occupancy Trajectories of Forest Birds on Montserrat", theme = theme(plot.title = element_text(size = 20))) +
  theme(plot.margin = margin(r = 15, b = 15))
occu_plots[[13]]
names(occu_plots)[13] <- 'ALL_species_aligned'
ggsave(filename = 'plot/ALL_species_occupancy_plot_aligned.jpeg', 
       plot = occu_plots$ALL_species, scale = 1, dpi = 'retina', units = 'px', 
       width = 4000, height = 4500)

occu_plots[[14]] <- annotate_figure(occu_plots$ALL_species_aligned + theme(plot.margin = margin(r = 15, b = 15, l = 15, t = 15)), left = textGrob('Occupancy', rot = 90, vjust = 1, gp = gpar(cex = 1.3)), bottom = textGrob("Year", vjust = 0, gp = gpar(cex = 1.3)))
occu_plots[[14]]
names(occu_plots)[14] <- 'ALL_species_aligned_axis_title'
ggsave(filename = 'plot/ALL_species_occupancy_plot_aligned_axis_title.jpeg', 
       plot = occu_plots$ALL_species_aligned_axis_title, scale = 1, dpi = 'retina', units = 'px', 
       width = 4000, height = 4500)

# put everything together but without align plots 
(occu_plots[[15]] <- ggarrange(occu_plots$ACHU + rremove('xylab') + rremove('x.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)), 
                               occu_plots$BANA + rremove('xylab') + rremove('xy.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16)), 
                               occu_plots$BRQD + rremove('xylab') + rremove('xy.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16)), 
                               occu_plots$CAEL + rremove('xylab') + rremove('x.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)), 
                               occu_plots$FOTH + rremove('xylab') + rremove('xy.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16)),
                               occu_plots$GTCA + rremove('xylab') + rremove('xy.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16)),
                               occu_plots$MTOR + rremove('xylab') + rremove('x.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)),
                               occu_plots$PETH + rremove('xylab') + rremove('xy.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16)),
                               occu_plots$PTCA + rremove('xylab') + rremove('xy.text') + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'), plot.title = element_text(vjust = -4, size = 16)),
                               occu_plots$SBTH + rremove('xylab') + theme(plot.margin = unit(c(0,0.2,0.2,1), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)),
                               occu_plots$SNPI + rremove('xylab') + rremove('y.text') + theme(plot.margin = unit(c(0,0.2,0.2,1), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)),
                               occu_plots$TREM + rremove('xylab') + rremove('y.text') + theme(plot.margin = unit(c(0,0.2,0.2,1), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)), 
                               ncol = 3, nrow = 4) + 
    plot_annotation(title = "Occupancy Trajectories of Forest Birds on Montserrat", theme = theme(plot.title = element_text(size = 20))) +
    theme(plot.margin = margin(r = 15, b = 15)))
names(occu_plots)[15] <- 'ALL_species_not_aligned'
ggsave(filename = 'plot/ALL_species_occupancy_plot_not_aligned.jpeg', 
       plot = occu_plots$ALL_species, scale = 1, dpi = 'retina', units = 'px', 
       width = 4000, height = 4500)

(occu_plots[[16]] <- annotate_figure(occu_plots$ALL_species_not_aligned + theme(plot.margin = margin(r = 15, b = 15, l = 15, t = 15)), left = textGrob('Occupancy', rot = 90, vjust = 1, gp = gpar(cex = 1.3)), bottom = textGrob("Year", vjust = 0, gp = gpar(cex = 1.3))))
names(occu_plots)[16] <- 'ALL_species_not_aligned_axis_title'
ggsave(filename = 'plot/ALL_species_occupancy_plot_aligned_axis_title.jpeg', 
       plot = occu_plots$ALL_species_aligned_axis_title, scale = 1, dpi = 'retina', units = 'px', 
       width = 4000, height = 4500)

#### 5. merge graphics of col and ext plots for all species with alt effect #####

# load .RDS data containing occupancy plots 
plotnames_col_ext <- list.files(path = 'plot/data/', pattern = 'ext_col_plot.rds')
col_ext_plots <- list()
for(i in 1:length(plotnames_col_ext)){
  col_ext_plots[[i]] <- readRDS(file = paste0('plot/data/', plotnames_col_ext[i]))}
names(col_ext_plots) <- str_remove(plotnames_col_ext, pattern = '_ext_col_plot.rds')

# create data frame with data for legend plot
legend_data <- data.frame(Altitude = seq(from = 132, to = 740, by = 0.2), 
                          Predicted = c(rep(0.1, times = length(seq(from = 132, to = 740, by = 0.2))), rep(0.625, times = length(seq(from = 132, to = 740, by = 0.2)))), 
                          Type = factor(c(rep('Colonization', times = length(seq(from = 132, to = 740, by = 0.2))), rep('Extinction', times = length(seq(from = 132, to = 740, by = 0.2))))), 
                          upper = c(rep(0.37, times = length(seq(from = 132, to = 740, by = 0.2))), rep(0.97, times = length(seq(from = 132, to = 740, by = 0.2)))), 
                          lower = c(rep(0.02, times = length(seq(from = 132, to = 740, by = 0.2))), rep(0.45, times = length(seq(from = 132, to = 740, by = 0.2))))) 
legend_data <- legend_data %>%
  mutate(Altitude = if_else(Type == 'Extinction' & Altitude > 525, NA, Altitude))

# create legend plot
legend_plot <- legend_data %>%
    ggplot(aes(x = Altitude, y = Predicted, colour = Type, fill = Type)) +
    geom_line(linewidth = 3) +
    geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.2) +
    geom_segment(aes(x = 660, y = 0.4, xend = 660, yend = 0.135), arrow = arrow(length = unit(0.2, units = 'cm')), color = 'black', lwd = 0.5) +
    annotate('text', y = 0.5, x = 660, label = 'Predicted \nvalues', size = 4, color = 'black', hjust = 0.5, vjust = 0.5) +
    geom_segment(aes(x = 535, y = 0.97, xend = 545, yend = 0.97), color = 'black', lwd = .5) +
    geom_segment(aes(x = 535, y = 0.45, xend = 545, yend = 0.45), color = 'black', lwd = .5) +
    geom_segment(aes(x = 545, y = 0.97, xend = 545, yend = 0.45), color = 'black', lwd = .5) + 
    annotate('text', y = 0.75, x = 555, label = '95% CI', size = 4, color = 'black', hjust = 0) +
    annotate('text', y = 0.21, x = 370, label = 'Colonization', size = 9, color = '#440154') +
    annotate('text', y = 0.735, x = 325, label = 'Extinction', size = 9, color = '#F1E51D') +
    
    labs(x="Altitude (in m a.s.l.)", y="Predicted Probability", title='Legend') + # change Name
    scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) + # set begin/end to 0.98 for yellow/ext // 0 for purple/col
    scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
    scale_y_continuous(limits = c(0, 1)) + 
    
    theme(panel.background=element_rect(fill="white", colour="black"),
          plot.margin = unit(c(50, 70, 40, 50), "pt"),
          plot.title=element_text(size=20, color='black', margin = margin(b=20)),
          axis.text=element_text(size=15, color="black"), 
          axis.title.y=element_text(size=18, margin = margin(r=15)),
          axis.title.x=element_text(size=18, margin = margin(t=15)),
          legend.text=element_text(size=15),
          legend.title = element_text(size=18),
          legend.position=c(0.2,0.88),
          panel.grid.major = element_line(linewidth=.1, color="grey94"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, colour = "black"))
legend_plot

# col_ext_legend <- get_legend(readRDS(file = 'plot/data/SBTH_ext_col_plot.rds'))

# change properties of plots
(ACHU <- col_ext_plots$ACHU + rremove('xylab') + rremove('x.text') + rremove('legend') + ggtitle(label = 'Antillean Crested Hummingbird') + theme(plot.margin = unit(c(0,0.4,0,0.4), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)))
(BRQD <- col_ext_plots$BRQD + rremove('xylab') + rremove('xy.text') + rremove('legend') + ggtitle(label = 'Bridled Quail-Dove') + theme(plot.margin = unit(c(0,0.4,0,0.4), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)))
(FOTH <- col_ext_plots$FOTH + rremove('xylab') + rremove('xy.text') + rremove('legend') + ggtitle(label = 'Forest Thrush') + theme(plot.margin = unit(c(0,0.4,0,0.4), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)))
(MTOR <- col_ext_plots$MTOR + rremove('xylab') + rremove('x.text') + rremove('legend') + ggtitle(label = 'Montserrat Oriole') + theme(plot.margin = unit(c(0,0.4,0,0.4), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)))
(PETH <- col_ext_plots$PETH + rremove('xylab') + rremove('xy.text') + rremove('legend') + ggtitle(label = 'Pearly-Eyed Thrasher') + theme(plot.margin = unit(c(0, 0.4, 0, 0.4), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)))
(PTCA <- col_ext_plots$PTCA + rremove('xylab') + rremove('xy.text') + rremove('legend') + ggtitle(label = 'Purple-Throated Carib') + theme(plot.margin = unit(c(0, 0.4, 0, 0.4), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)))
(SBTH <- col_ext_plots$SBTH + rremove('xylab') + rremove('legend') + ggtitle(label = 'Scaly-Breasted Thrasher') + theme(plot.margin = unit(c(0, 0.4, 0, 0.4), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)))
(SNPI <- col_ext_plots$SNPI + rremove('xylab') + rremove('y.text') + rremove('legend') + ggtitle(label = 'Scaly-Naped Pigeon') + theme(plot.margin = unit(c(0, 0.4, 0, 0.4), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)))

# put everything together : Version aligned 
col_ext_plots[[9]] <- cowplot::plot_grid(ACHU, BRQD, FOTH, MTOR, PETH, PTCA, SBTH, SNPI,
                                          legend_plot + rremove('xylab') + rremove('y.text') + rremove('legend') + 
                                            theme(plot.margin = unit(c(0,0.4,0,0.4), 'lines'), 
                                                  plot.title = element_text(vjust = -4, size = 16), 
                                                  axis.text = element_text(size = 12)), 
                                          align = 'hv', labels = c('a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)')) +
    theme(plot.margin = margin(r = 15, b = 15))
col_ext_plots[[9]]
names(col_ext_plots)[9] <- 'All_species_aligned'

col_ext_plots[[10]] <- annotate_figure(col_ext_plots[[9]] + theme(plot.margin = margin(r = 15, b = 15, l = 15, t = 15)), left = textGrob('Predicted Probability', rot = 90, vjust = 1, gp = gpar(cex = 1.3)), bottom = textGrob("Elevation (in m a.s.l.)", vjust = 0, gp = gpar(cex = 1.3)))
col_ext_plots[[10]]
names(col_ext_plots)[10] <- 'All_species_aligned_axis_title'
ggsave(filename = 'plot/ALL_species_col_ext_plot_aligned_axis_title.jpeg', 
       plot = col_ext_plots$All_species_aligned_axis_title, height = 2.83*3.5, width = 3.28*3.5, dpi = 'retina')

# delete plots that does not make sense
(ACHU <- col_ext_plots$ACHU + rremove('xylab') + rremove('x.text') + rremove('legend') + ggtitle(label = 'Antillean Crested Hummingbird') + theme(plot.margin = unit(c(0,0.4,0,0.4), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)))
(BRQD <- col_ext_plots$BRQD + rremove('xylab') + rremove('xy.text') + rremove('legend') + ggtitle(label = 'Bridled Quail-Dove') + theme(plot.margin = unit(c(0,0.4,0,0.4), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)))
(MTOR <- col_ext_plots$MTOR + rremove('xylab') + rremove('xy.text') + rremove('legend') + ggtitle(label = 'Montserrat Oriole') + theme(plot.margin = unit(c(0,0.4,0,0.4), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)))
(PTCA <- col_ext_plots$PTCA + rremove('xylab') + rremove('legend') + ggtitle(label = 'Purple-Throated Carib') + theme(plot.margin = unit(c(0, 0.4, 0, 0.4), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)))
(SNPI <- col_ext_plots$SNPI + rremove('xylab') + rremove('y.text') + rremove('legend') + ggtitle(label = 'Scaly-Naped Pigeon') + theme(plot.margin = unit(c(0, 0.4, 0, 0.4), 'lines'), plot.title = element_text(vjust = -4, size = 16), axis.text = element_text(size = 12)))

# put everything together : Version aligned 
col_ext_plots[[11]] <- cowplot::plot_grid(ACHU, BRQD, MTOR, PTCA, SNPI,
                                         legend_plot + rremove('xylab') + rremove('y.text') + rremove('legend') + 
                                           theme(plot.margin = unit(c(0,0.4,0,0.4), 'lines'), 
                                                 plot.title = element_text(vjust = -4, size = 16), 
                                                 axis.text = element_text(size = 12)), 
                                         align = 'hv', labels = c('a)', 'b)', 'c)', 'd)', 'e)')) +
  theme(plot.margin = margin(r = 15, b = 15))
col_ext_plots[[11]]
names(col_ext_plots)[11] <- 'All_species_aligned_short'

col_ext_plots[[12]] <- annotate_figure(col_ext_plots[[11]] + theme(plot.margin = margin(r = 15, b = 15, l = 15, t = 15)), left = textGrob('Predicted Probability', rot = 90, vjust = 1, gp = gpar(cex = 1.3)), bottom = textGrob("Elevation (in m a.s.l.)", vjust = 0, gp = gpar(cex = 1.3)))
col_ext_plots[[12]]
names(col_ext_plots)[12] <- 'All_species_aligned_short_axis_title'
ggsave(filename = 'plot/ALL_species_col_ext_plot_aligned_short_axis_title.jpeg', 
       plot = col_ext_plots$All_species_aligned_short_axis_title, height = 2*3.5, width = 3.28*3.5, dpi = 'retina')

####### Occupancy plot
occupancy_data <- data.frame(Year = 2011:2023, Occupancy = 0.5, lower_cl = 0.4, upper_cl = 0.65)

occupancy_plot <- ggplot(data = occupancy_data, aes(x = Year, y = Occupancy)) +
  geom_point(size = 3, pch = 16, colour="firebrick") +
  geom_errorbar(aes(ymin = lower_cl,ymax = upper_cl), width = 0.2) +
  scale_y_continuous(limits = c(0, 1.1), breaks=seq(0,1,0.2),labels=seq(0,1,0.2))  +
  scale_x_continuous(limits = c(2010, 2024), breaks=seq(2011,2023,2),labels=seq(2011,2023,2)) +
  labs(x="Year", y="Mean occupancy probability", title='Bridled Quail-Dove') +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  
  theme(panel.background=element_rect(fill="white", colour="black"),
        plot.margin = unit(c(50, 70, 40, 50), "pt"),
        plot.title=element_text(size=20, color='black', margin = margin(b=20)),
        axis.text=element_text(size=15, color="black"), 
        axis.title.y=element_text(size=18, margin = margin(r=15)),
        axis.title.x=element_text(size=18, margin = margin(t=15)),
        legend.text=element_text(size=15),
        legend.title = element_text(size=18),
        legend.position=c(0.2,0.88),
        panel.grid.major = element_line(linewidth =.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black")) +
  guides(color = guide_legend(title = "Color", override.aes = list(size = 3, pch = 16)),
         fill = guide_legend(title = "Fill", override.aes = list(alpha = 0.3, shape = NA)))
occupancy_plot

legend <- get_legend(occupancy_plot)

#### 6. access best model objects ####

# read data in a list 
best_models <- vector(mode = 'list', length = length(SPECIES)) # create input vector with mode 'list'
for(s in 1:length(SPECIES)){
  best_models[[s]] <- readRDS(file = sprintf('data/best_models/%s_best_model.rds', SPECIES[s]))
  names(best_models)[s] <- SPECIES[s]}

# extract data and back transform to display in a data.frame 


#### 7. access best model estimates ####

# read data
estimates_full <- vector(mode = 'list', length = length(SPECIES)) # create input vector with mode 'list'
for(s in 1:length(SPECIES)){
  estimates_full[[s]] <- data.frame(read.csv(file = sprintf('data/best_model_output_estimates/%s_best_model_output_estimates_all.csv', SPECIES[s]), header = T, sep = ','))
  names(estimates_full)[s] <- SPECIES[s]}

# assign column with species to each table, exclude MTOR and BANA
for(i in 1:length(SPECIES[1:(length(SPECIES)-1)])){
  estimates_full[[i]] <- estimates_full[[i]] %>% mutate(Species = SPECIES[i]) %>%
    select(Species, everything())}

# make table with all the estimates 
estimates_table <- data.frame(estimates_full[[1]]) # initial data.frame with first species as input data.frame
for(p in 2:length(SPECIES)){
  estimates_table <- rbind(estimates_table, estimates_full[[p]])}

# View(estimates_table %>% filter(Component == 'Extinction' & Parameter == 'alt'))
# View(estimates_table %>% filter(Parameter %in% c('alt', 'alt:treeheight')))
# View(modsel_full_species %>% filter(Step == 'g_e') %>% group_by(Species) %>% slice(which.min(AIC)) %>% ungroup())

##### 9. load col ext data with altitude for each species ####

# load col_ext data in a list
col_ext_names <- list.files(path = 'data/pred_col_ext/', pattern = 'pred_ext_col.csv')
col_ext_data <- vector(mode = 'list', length = length(col_ext_names)) # create input vector with mode 'list'
for(i in 1:length(col_ext_names)){
  col_ext_data[[i]] <- data.frame(read.csv(file = paste0('data/pred_col_ext/', col_ext_names[i]), header = T, sep = ','))
  names(col_ext_data)[i] <- str_remove(col_ext_names[i], pattern = '_pred_ext_col.csv')}

# create data.frame containing all species occupancy data per year
col_ext_table <- data.frame(col_ext_data[[1]], Species = str_remove(col_ext_names[1], pattern = '_pred_ext_col.csv')) # initial data.frame with first species 
for(s in 2:length(col_ext_names)){
  col_ext_table  <- rbind(col_ext_table, col_ext_data[[s]] %>% mutate(Species = str_remove(col_ext_names[s], pattern = '_pred_ext_col.csv')))}

# very simple plot
str(col_ext_table)
col_ext_table <- col_ext_table %>% mutate(Type = factor(Type), Species = factor(Species))
col_ext_table %>% 
  ggplot() +
  geom_line(mapping = aes(y = Predicted, x = Altitude, colour = Type)) + 
  facet_wrap(~Species) + 
  theme_bw()

View(col_ext_table %>% filter(Species == 'PETH') %>% mutate(Predicted = format(Predicted, scientific = F)))

col_ext_table %>% filter(Species == 'SNPI') %>%
  ggplot() +
  geom_line(mapping = aes(y = Predicted, x = Altitude, colour = Type)) +
  theme_bw()



##### 10. load occupancy data per species and year ####

occu_all <- vector(mode = 'list', length = length(SPECIES_sorted)) # create input vector with mode 'list'
for(s in 1:length(SPECIES_sorted)){
  occu_all[[s]] <- data.frame(read.csv(file = sprintf('data/ranef_occupancy_data/%s_occupancy_data.csv', SPECIES_sorted[s]), header = T, sep = ','))
  names(occu_all)[s] <- SPECIES_sorted[s]}

# create data.frame containing all species occupancy data per year
occu_full_species <- data.frame(occu_all[[1]], Species = SPECIES_sorted[1]) # initial data.frame with first species 
for(s in 2:length(SPECIES_sorted)){
  occu_full_species  <- rbind(occu_full_species, occu_all[[s]] %>% mutate(Species = SPECIES_sorted[s]))}

# make some adjustments 
mean <- occu_full_species %>% group_by(Species) %>% summarize(Mean = format(mean(as.numeric(Occupancy)), digits = 4)) %>% select(Mean)
occu_full_species <- occu_full_species %>% 
  select(Species, Year, Occupancy, -upper_cl, -lower_cl) %>% 
  mutate(Species = factor(Species)) %>%
  format(scientific = F)
occu_full_species$Occupancy <- round(as.numeric(occu_full_species$Occupancy), digits = 4)
occu_full_species <- rbind(occu_full_species, data.frame(Species = SPECIES_sorted, Year = 2020, Occupancy = NA))

# wide format and arrange by year
pretty_occu_table <- pivot_wider(occu_full_species, names_from = Species, values_from = Occupancy, id_cols = Year) %>% 
  arrange(Year) %>% 
  mutate_all(~ifelse(is.na(.), "-", .))

# calculate mean 
mean_values <- pretty_occu_table %>%
  summarise(across(-Year, ~format(round(mean(as.numeric(.), na.rm = TRUE), 4), nsmall = 4))) %>%
  mutate(Year = "Overall Mean")

# add mean to table
pretty_occu_table <- bind_rows(pretty_occu_table, mean_values)

# make these funny gt tables
gt_occu = gt(pretty_occu_table)
gt_occu <- gt_occu  %>%
  tab_header(
    title = md("Mean Occupancies for all Years")) %>%
  tab_spanner(
    label = "Species",
    columns = SPECIES_sorted) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")),
    locations = list(cells_column_labels(
      columns = everything()), cells_title())) %>%
  cols_align(
    align = "right",
    columns = SPECIES) %>%
  cols_align(
    align = "center",
    columns = 'Year') %>%
  tab_footnote(
    footnote = "No bird surveys were conducted in 2020.",
    locations = cells_body(
      columns = c(1),
      rows = 10)) %>% 
  tab_source_note(
    source_note = md("**Mean Occupancy**: Derived using the ranef() function from best models (see Table 3) for each species and summarized as mean above all survey points and each year.<br>
    **Overall Mean**: Averaged occupancy for each species above all years, except of 2020 where no bird surveys where conducted.<br>
    **Species**: ACHU = Antillean Crested Hummingbird; BANA = Bananaquit; BRQD = Bridled Quail-Dove; 
    CAEL = Caribbean Elaenia; FOTH = Forest Thrush; GTCA = Green-Throated Carib; MTOR = Montserrat Oriole; 
    PETH = Pearly-Eyed Thrasher; PTCA = Purple-Throated Carib; SBTH = Scaly-Breasted Thrasher; SNPI = Scaly-Naped Pigeon; 
    TREM = Brown Trembler.")) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('top'), color = 'black', weight = px(2.25))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c('bottom'), color = 'black', weight = px(1.5))),
    locations = list(cells_title())) %>% 
  tab_style(
    style = list(
      cell_borders(sides = c("bottom"), color = 'black', weight = px(2.25))),
    locations = list(cells_column_labels(), cells_body(rows = c(13,14)))) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = 14))  # last row bold 
gt_occu

#save the table (a bit complicated, but ok): create an html file first, that is opened by webshot2 and then takes a photo and stores it as .png
html_file = "OCCU_TABLE_original.html"
gtsave(gt_occu, html_file, path = paste0(getwd(), "/plot/tables/")) 
webshot2::webshot(paste0(getwd(), "/plot/tables/", html_file), 
                  file = paste0(getwd(), "/plot/tables/", "OCCU_TABLE.png"), 
                  vwidth = 1500,
                  vheight = 1200, zoom = 3)
webshot2::webshot(paste0(getwd(), "/plot/tables/", html_file), 
                  file = paste0(getwd(), "/plot/tables/", "OCCU_TABLE_small.png"), 
                  vwidth = 1300,
                  vheight = 1200, zoom = 3)
webshot2::webshot(paste0(getwd(), "/plot/tables/", html_file), 
                  file = paste0(getwd(), "/plot/tables/", "OCCU_TABLE_narrow.png"), 
                  vwidth = 900,
                  vheight = 1200, zoom = 3)

# calculate gain in occupancy in % and 0...
occu_full_species$Occupancy[occu_full_species$Species == 'GTCA'][9]*100 / occu_full_species$Occupancy[occu_full_species$Species == 'GTCA'][4] -100
occu_full_species$Occupancy[occu_full_species$Species == 'GTCA'][9] - occu_full_species$Occupancy[occu_full_species$Species == 'GTCA'][4]

# 
round(unmarked::coef(best_models$PETH), digits = 4)
round(unmarked::confint(best_models$PETH, type = 'ext'), digits = 4)

####### FOR JOHANNES

# GOF 
GOF_mb$MTOR # assess GOF_test_objects
GOF_df %>% arrange(Species) # c-hat estimates and p-values for all species

# effect sizes 
best_models$MTOR # assess model summary
View(estimates_table) # all effect sizes on logit and prob scale





