
##### Montserrat Forest Bird Counts - data import and analysis with 'unmarked' #####
##### written in June 2023 by  Filibert Heim, filibert.heim@posteo.de          #####

##### 1: load required packages ####

# install.packages('scales')
library(scales)
# install.packages('data.table')
library(data.table)
library(reshape)
library(lubridate)
library(unmarked)
library(AICcmodavg) # package for model selection and goodness-of-fit tests
library(MuMIn)
library(tidyverse)

##### 2: load the prepared data and make last preparation #####

# set working directory and load prepared data from Steffen
setwd('C:/Users/filib/Documents/Studium/Bachelorarbeit/R_BachelorThesisMontserrat')
# setwd('C:/Users/sop/Documents/Steffen/RSPB/Montserrat') # for Steffen
load(file = 'data/MONTSERRAT_ANNUAL_DATA_INPUT2023.RData')
# load(file = 'MONTSERRAT_ANNUAL_DATA_INPUT2023.RData') # for Steffen

# save 'filter', 'rename' and 'select' as functions from dplyr to avoid problems
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename

###### 2.1: set YEAR and SPECIES that should be analysed ####
YEAR <- 2023 # fill in the most recent year that data is available for
SPECIES <- c('BANA') # fill in species the analysis should be made for
# SPECIES: 'MTOR','FOTH','BRQD','TREM','ACHU','PTCA','PETH','GTCA','SBTH','SNPI','CAEL','BANA'

###### 2.2: check the data and remove unneeded stuff #####

head(activity) # (tidyverse-)data.frame with calculated bird-activity data
head(add) # I don't know what kind of data there is inside
head(birds) # data like its provided by the database
head(birds2023) # bird data from 2023
head(COUNTDATA) # completely prepared bird observation data for analysis!
head(MTOR_count) # RODBC connection to data base
head(nyears) # number of years with bird data from the monitoring 
head(obsCov) # prepared obsCov data, rename columns and remove 'Rain', 'Year' (2nd) - but also without wind data!
head(pointmax) # (tidyverse-)data.frame with extracted maximal activity per point
head(siteCov) # prepared siteCov data for analysis, but there are mistakes in column names!
head(siteCov_db) # loaded siteCovs from database to check colnames of siteCov
head(species) # species list and abundance over all surveys
bird_names <- species
rm(species)
print(SPECIES)
head(startdate) # first date in every year to calculate numeric dates
head(summary2023) # hohle abundance over all visits in 2023
head(SURVEYDATA) # nearly the same data.frame like obsCOV
head(surveys2023) # obsCov from all visits in 2023
head(table1) # I don't really know what kind of data this data.frame contains
head(totals2023) # I don't really know what kind of data this data.frame contains
print(y) # last recent year 
head(rain_y)

###### 2.3: remove unneeded object from global environment

rm(list = c('surveys2023', 'pointmax', 'table1', 'totals2023', 'birds', 
            'birds2023', 'activity', 'summary2023', 'nyears', 'MTOR_count', 
            'SURVEYDATA', 'rain_y', 'add', 'startdate', 'starttime', 'y', 
            'removal'))

###### 2.4: relable variate names/columns #####

obsCov <- select(obsCov, -Rain, -Year)
names(obsCov)[1:10] <- c('year','point','count','obs','skill','rain','wind','day',
                         'time','activity')

names(siteCov)[c(1:9, 17, 18, 21)] <- c('point', 'eastings', 'northings', 'habitat', 
                                        'dbh', 'distance', 'treeheight', 'elevation', 
                                        'slope', 'canopy', 'location', 'alt') # go on with renaming

names(COUNTDATA) <- c('species', 'year', 'point', 'count', 'N')

###### 2.5: digression to calculate missing wind values as obsCov ####

obsCov_db <- obsCov_db %>% 
  select(-Time, -Rain, -MinOfskill, -Date, -CountOfObservername) %>% # delete unneeded columns
  dplyr::rename(point = Point, count = Count) 
obsCov_db$Wind <- obsCov_db$Wind %>% 
  recode('calm' = 1, 'light' = 2, 'moderate' = 3, 'strong' = 4) # convert categories 'calm', light', 'moderate' and 'strong' to numbers
obsCov_db$Wind <- factor(obsCov_db$Wind, levels = c(1, 2, 3, 4), ordered = TRUE) # store as ordered factor with 4 levels!
obsCov <- obsCov %>% # transfer wind data from data base (obsCov_db) to prepared data.frame
  left_join(x = obsCov, y = obsCov_db, by = c('year', 'point', 'count')) %>%
  select(-wind) %>% 
  rename(wind = Wind)
obsCov$wind[c(which(is.na(obsCov$wind)))] <- as.factor(1) # replace all the NA values by 1 for 'calm' - this is not really good, but I cant find another workaround

###### 2.6: digression to postpone surveys if some weren't made ######

# filter for SPECIES and align with obsCov to get information from surveys that weren't made
COUNTDATA <- COUNTDATA %>% filter(species == SPECIES) %>% # choose species
  left_join(obsCov %>% mutate(point = as.numeric(point)), by = c('point', 'count', 'year')) %>%
  mutate(N = ifelse(obs == 'NA', 'NA', N)) %>% # insert 'NA' to COUNTDATA were no surveys were made
  select(-obs, -skill, -rain, -time, -activity, -wind) # remove unneeded columns

# transfer information of surveys that weren't made also to obsCov data
obsCov <- COUNTDATA %>%
  left_join(obsCov %>% mutate(point = as.numeric(point)), by = c("point", 'count', 'year')) %>%
  select(-species, -N, -day.x) %>%
  rename(day = day.y)

# filter COUNTDATA for all the surveys in one year at one point, sort by day and rename counts
COUNTDATA <- COUNTDATA %>% 
  group_by(year, point) %>% 
  arrange(day, .by_group = TRUE) %>% 
  mutate(count = as.numeric(rep(1:3))) %>% 
  select(-day) %>% 
  ungroup()

# filter obsCov for all the surveys in one year at one point, sort by day and rename counts
obsCov <- obsCov %>% 
  group_by(year, point) %>% 
  arrange(day, .by_group = TRUE) %>% 
  mutate(count = as.numeric(rep(1:3))) %>% 
  ungroup()

##### 3: provide and prepare data for unmarkedMultFrame ####

numPrimary <- length(unique(COUNTDATA$year)) # calculates number of years

occdataColext <- COUNTDATA %>%
  mutate(occupancy = ifelse(N > 0, 1, 0)) %>% # convert in detection/non-detection
  mutate(season = paste(year, count, sep = '_')) %>% # connect year, count to string
  select(-species, -N, -year, -count) %>% # remove unneeded columns 
  spread(key = season, value = occupancy) %>% # spread a key-value pair across multiple columns
  arrange(point) # sort rows in order of point number

siteCovsColext <- siteCov %>% mutate(point = as.numeric(point), ridge = as.factor(ridge)) %>% # code ridge as.factor
  arrange(point) # sort rows in order of point number

obsCovsColext <- COUNTDATA %>% filter(species == SPECIES) %>%
  left_join((obsCov %>% mutate(point = as.numeric(point))), by = c('point', 'count', 'year')) %>% 
  select(-N, -species, -skill) %>% # remove columns
  arrange(point, year, count) %>% # sort rows in order of year, point, count
  rename(season = year) %>% # rename year to avoid confusion with yearly-siteCovs
  mutate(rain = as.numeric(rain)) # code rain as.numeric

yearlySiteCovsYear <- COUNTDATA %>% filter(species == SPECIES) %>% 
  group_by(point,year) %>%
  summarise(occu = as.numeric(mean(year) - 2010)) %>%
  spread(key = year, value = occu) %>% # spread a key-value pair across multiple columns
  arrange(point)
yearlySiteCovsColext <- list(year=yearlySiteCovsYear[,2:14])
  
###### 3.2: check dimensions and input data.frames####
head(siteCovsColext)
head(obsCovsColext)
head(yearlySiteCovsColext)
head(occdataColext)

dim(occdataColext)
dim(siteCovsColext)
dim(obsCovsColext)/(3*numPrimary)
dim(yearlySiteCovsYear)

###### 3.3: create unmarkedMultFrame ####

umf <- unmarkedMultFrame(y = occdataColext[,2:40], siteCovs = siteCovsColext,yearlySiteCovs=yearlySiteCovsColext, 
                         obsCovs = obsCovsColext, numPrimary = numPrimary)
summary(umf) # looks good
str(umf) # looks good, but pay attention: not all variables from this data set are coded in an appropriate way (as.factor())

###### 3.4: scale numeric variables to avoid fitting problems ####

siteCovs(umf)[c(2,3,5:17,19,21:37)] <- scale(siteCovs(umf)[c(2,3,5:17,19,21:37)]) 
obsCovs(umf)[c(5:8)] <- scale(obsCovs(umf)[c(5:8)])
yearlySiteCovs(umf)[1] <- scale(yearlySiteCovs(umf)[1])
str(umf)
summary(umf)

##### 4: multi-season, dynamic one species occupancy model ####

###### 4.1: fit models for detection probability p() first for modSel ####

# fit models for detection probability p() manually, go on with fitList(), modSel()
fm1 <- colext(~1, ~1, ~1, ~1, data = umf, se = T)
fm2 <- colext(~1, ~1, ~1, ~day, data = umf, se = T)
fm3 <- colext(~1, ~1, ~1, ~day+time, data = umf, se = T)
fm4 <- colext(~1, ~1, ~1, ~day+time+rain, data = umf, se = T)
fm5 <- colext(~1, ~1, ~1, ~day+time+rain+wind, data = umf, se = T)
fm6 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity, data = umf, se = T)
fm7 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)

# organise models in fitList before printing the modSel-table
p_fitList <- fitList('psi(.)g(.)e(.)p(.)' = fm1, 
                     'psi(.)g(.)e(.)p(~day)' = fm2, 
                     'psi(.)g(.)e(.)p(~day+time)' = fm3, 
                     'psi(.)g(.)e(.)p(~day+time+rain)' = fm4,
                     'psi(.)g(.)e(.)p(~day+time+rain+wind)' = fm5,
                     'psi(.)g(.)e(.)p(~day+time+rain+wind+activity)' = fm6,
                     'psi(.)g(.)e(.)p(~day+time+rain+wind+activity+ridge)' = fm7)
(p_modSel <- modSel(p_fitList))
# best submodel for p(): ~day+time+rain+wind+activity+ridge, AIC difference to second best 3.12 - go on with this best one

# create data.frame from modSel table
p_modSel_df <- as(p_modSel, Class = 'data.frame') %>% 
  mutate(step = 'p')

###### 4.2: fit models for initial occupancy psi() first for modSel ####

fm8 <- colext(~alt, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# add treeheight
fm9 <- colext(~treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm10 <- colext(~alt+treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# add dbh
fm11 <- colext(~dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm12 <- colext(~alt+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm13 <- colext(~treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm14 <- colext(~alt+treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
# add interaction between alt and treeheight
fm15 <- colext(~alt:treeheight, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm16 <- colext(~alt:treeheight+dbh, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)

# put the fitted models in a fitList() and rank them by AIC in modSel()
psi_fitList <- fitList('psi(.)g(.)e(.)p(~day+time+rain+wind+activity+ridge)' = fm7, 
                       'psi(~alt)g(.)e(.)p(~day+time+rain+wind+activity+ridge)' = fm8, 
                       'psi(~treeheight)g(.)e(.)p(~day+time+rain+wind+activity+ridge)' = fm9, 
                       'psi(~alt+treeheight)g(.)e(.)p(~day+time+rain+wind+activity+ridge)' = fm10,
                       'psi(~dbh)g(.)e(.)p(~day+time+rain+wind+activity+ridge)' = fm11, 
                       'psi(alt+dbh)g(.)e(.)p(~day+time+rain+wind+activity+ridge)' = fm12, 
                       'psi(treeheight+dbh)g(.)e(.)p(~day+time+rain+wind+activity+ridge)' = fm13, 
                       'psi(~alt+treeheight+dbh)g(.)e(.)p(~day+time+rain+wind+activity+ridge)' = fm14, 
                       'psi(alt:treeheight)g(.)e(.)p(~day+time+rain+wind+activity+ridge)' = fm15, 
                       'psi(~alt:treeheight+dbh)g(.)e(.)p(~day+time+rain+wind+activity+ridge)' = fm16)
print(psi_modSel <- modSel(psi_fitList)) 
# best sub-model for psi(): ~1, AIC difference to second best is 1.87 (~alt) and third best is 1.96 (alt:treeheight) - go on with the best one

# create data.frame from modSel table
psi_modSel_df <- as(psi_modSel, Class = 'data.frame') %>% 
  mutate(step = 'psi')

###### 4.3: fit models for extinction and colonization probability for modSel ####

fm17 <- colext(~1, ~1, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm18 <- colext(~1, ~alt, ~1, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm19 <- colext(~1, ~1, ~alt, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm20 <- colext(~1, ~alt, ~alt, ~day+time+rain+wind+activity+ridge, data = umf, se = T)
fm21 <- colext(~1, ~year, ~year, ~day+time+rain+wind+activity+ridge, data = umf, se = T)  ## this model will exclude the possibility that observed changes are just annual changes

# put the fitted models in a fitList() and rank them by AIC in modSel()
g_e_fitList <- fitList(constant = fm17, expansion = fm18, contraction = fm19,
                            shift = fm20, year = fm21, nullmodel = fm1)
print(g_e_modSel <- modSel(g_e_fitList))
# best model: constant model (fm17) has lowest AIC, with difference in AIC of 2.00 units to expansion/contraction models

###### 4.4: export all required model selection tables ####

# create data.frame from modSel table
g_e_modSel_df <- as(g_e_modSel, Class = 'data.frame') %>% 
  mutate(step = 'g_e')

# merge data.frames from BANA modSel tables to one data frame and export it as .csv
modSel_export <- rbind(g_e_modSel_df %>% select(model, formula, step, nPars, AIC, delta, AICwt, cumltvWt), 
      psi_modSel_df %>% select(model, formula, step, nPars, AIC, delta, AICwt, cumltvWt), 
      p_modSel_df %>% select(model, formula, step, nPars, AIC, delta, AICwt, cumltvWt))
modSel_export <- modSel_export %>% mutate(species = SPECIES) %>% # add species name
  select(species, everything()) %>% # change order of columns
  arrange(match(step, c('p','psi','g_e')), AIC)
write.csv(modSel_export, file = sprintf('output/data/modSel/%s_modSel_full.csv', SPECIES)) # export as .csv file

# save best model as best_model
best_model <- fm17 # fill in best model
saveRDS(best_model, file = sprintf('output/data/best_models/%s_best_model.rds', SPECIES))

##### 5: Goodness-of-Fit Test - MacKenzie-Bailey Goodness-of-Fit Test by function mb.gof.test() ####

# calculate GoF-Test
(best_model_gof <- mb.gof.test(best_model, nsim = 1000)) ## p-value of 0, c-hat estimate of 2.74

# save output as rds data
saveRDS(best_model_gof, file = sprintf('output/data/GOF/%s_gof_mb.rds', SPECIES))

##### 6: view results and extract data from best single model ####

# view models
summary(best_model) # best model based on AIC

# get names of values from model
names(best_model)

# extract everything (psi, col, ext, det) - on the logit-scale!
toExportpsi <- summary(best_model@estimates@estimates$psi) %>%
  mutate(Parameter=names(best_model@estimates@estimates$psi@estimates)) %>%
  mutate(Component="Initial Occupancy") 

toExportdet <- summary(best_model@estimates@estimates$det) %>%
  mutate(Parameter=names(best_model@estimates@estimates$det@estimates)) %>%
  mutate(Component="Detection probability")

toExportcol <- summary(best_model@estimates@estimates$col) %>%
  mutate(Parameter=names(best_model@estimates@estimates$col@estimates)) %>%
  mutate(Component="Colonization")

toExportext <- summary(best_model@estimates@estimates$ext) %>%
  mutate(Parameter=names(best_model@estimates@estimates$ext@estimates)) %>%
  mutate(Component="Extinction")

toExport <- bind_rows(toExportpsi, toExportcol, toExportext, toExportdet)

# export everything backTransformed

toExport_transformed <- cbind(Estimates = plogis(coef(best_model)), 
                              SE_estimates_backTransformed = plogis(SE(best_model)),
                              rbind(plogis(confint(best_model, type = 'psi')), 
                                      plogis(confint(best_model, type = 'col')), 
                                      plogis(confint(best_model, type = 'ext')), 
                                      plogis(confint(best_model, type = 'det'))))

# merge two data.frames for export
toExport_final <- cbind(toExport, toExport_transformed) %>% 
  rename(Estimate_backTransformed = Estimates, 
         LCI_backTransformed = `0.025`, 
         UCI_backTransformed = `0.975`) %>% 
  mutate(Species = SPECIES) %>% 
  select(Species, Component, Parameter, everything())
toExport_final[,4:ncol(toExport_final)] <- round(toExport_final[,4:ncol(toExport_final)], digits = 4)
fwrite(toExport_final, sprintf('output/data/best_model_output_estimates/%s_best_model_output_estimates_all.csv', SPECIES))

##### 7: extraction of values and plotting ####

###### 7.1 plot relationship between colonization and extinction with altitude #####
# code mainly taken from: https://cran.r-project.org/web/packages/unmarked/vignettes/colext.html and steffen

# create input data.frame with values for the prediction, 
nd <- data.frame(day = 0, time = 0, rain = 0, wind = 1, activity = max(umf@obsCovs$activity, na.rm = T), # use maximum bird activity, lowest wind speed, mean of time and day = 0 used
                 ridge = 0, treeheight = 0, # ridge with valley or midslope used, mean scaled treeheight used, should be 0 (steffen used 15, why?)
                 dbh = 0, # mean scaled dbh used, should be 0 
                 alt = seq(from = min(umf@siteCovs$alt), 
                           to = max(umf@siteCovs$alt), by = 0.02),
                 altitude = rescale(seq(from = min(umf@siteCovs$alt), to = max(umf@siteCovs$alt), by = 0.02), 
                                    to = c(min(siteCovsColext$alt), max(siteCovsColext$alt)), # to = output range as c()
                                    from = c(min(umf@siteCovs$alt), max(umf@siteCovs$alt)))) # from = input range as c()
# alt: scaled altitude for prediction, altitude: rescaled altitude for plotting

#### no prediction needed because of intercept only model for ext and col submodels!

## predict values by input data from nd, add column with re-scaled altitude
#pred_ext <- predict(best_model, type = 'ext', newdata = nd) %>%
#  mutate(Type = 'Extinction', Altitude = nd$altitude)
#pred_col <- predict(best_model, type = 'col', newdata = nd) %>%
#  mutate(Type = 'Colonization', Altitude = nd$altitude)

## connect tables and save as one .csv for plotting
#pred_ext_col <- bind_rows(pred_ext, pred_col)
#fwrite(pred_ext_col, file = sprintf('output/data/pred_col_ext/%s_pred_ext_col.csv', SPECIES))

## plot colonization~altitude and extinction~altitude together
#pred_ext_col_plot <- pred_ext_col %>%
#  ggplot(aes(x = Altitude, y = Predicted, colour = Type, fill = Type)) +
#  geom_line(linewidth = 3) +
#  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.2) +
#  labs(x="Elevation (m above sea level)", y="Predicted probability", title='Bananaquit') + # change Name
#  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) + # set begin/end to 0.98 for yellow/ext // 0 for purple/col
#  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
#  scale_y_continuous(limits = c(0, 1)) + 
  
#  theme(panel.background=element_rect(fill="white", colour="black"),
#        plot.margin = unit(c(50, 70, 40, 50), "pt"),
#         plot.title=element_text(size=20, color='black', margin = margin(b=20)),
#        axis.text=element_text(size=15, color="black"), 
#        axis.title.y=element_text(size=18, margin = margin(r=15)),
#        axis.title.x=element_text(size=18, margin = margin(t=15)),
#        legend.text=element_text(size=15),
#        legend.title = element_text(size=18),
#        legend.position=c(0.2,0.88),
#        panel.grid.major = element_line(size=.1, color="grey94"),
#        panel.grid.minor = element_blank(),
#        panel.border = element_rect(fill=NA, colour = "black"))
#pred_ext_col_plot

## save as image 
#ggsave(filename = sprintf('output/plot/%s_ext_col_plot.jpg', SPECIES), 
#       plot = last_plot(), scale = 1, dpi = 'retina', units = 'px', 
#       width = 3000, height = 2200)
#saveRDS(pred_ext_col_plot, file = sprintf('output/plot/data/%s_ext_col_plot.rds', SPECIES))

###### 7.2 plot estimated occupancy for time series ####

# using empirical bayes estimates of occupancy averaged across all points
occupancy <- ranef(best_model)

# manipulate data structure, summarize and store data in data.frame
occupancy_data <- as.data.frame(bup(occupancy, stat = 'mean')) %>% # posterior mean by stat = ''
  gather(key = 'Year',value = 'occu') %>% # transfer from wide to long format 
  mutate(Year = as.numeric(str_replace(Year,'V',''))) %>% # delete V before number of season
  mutate(Year = Year + 2010) %>% # calculate year 
  group_by(Year) %>%
  summarise(Occupancy = mean(occu)) # calculate the mean occupancy about all points in one year 

# calculate confidence intervals and store them in occupancy_data data.frame for plotting
occupancy_confint <- confint(occupancy, level = 0.95) # 95% CI
occupancy_data$lower_cl <- apply(occupancy_confint[,1,],2, mean)
occupancy_data$upper_cl <- apply(occupancy_confint[,2,],2, mean)

# remove data from 2020 because in this year surveys didn't take place 
occupancy_data <- occupancy_data %>% filter(!Year == '2020')

# save occupancy_data for later plotting
fwrite(occupancy_data, file = sprintf('output/data/ranef_occupancy_data/%s_occupancy_data.csv', SPECIES))

# plot occupancy over time from occupancy_data data.frame 
occupancy_plot <- ggplot(data = occupancy_data, aes(x = Year, y = Occupancy)) +
  geom_point(size = 3, pch = 16, colour="firebrick") +
  geom_errorbar(aes(ymin = lower_cl,ymax = upper_cl), width = 0.2) +
  scale_y_continuous(limits = c(0, 1.1), breaks=seq(0,1,0.2),labels=seq(0,1,0.2))  +
  scale_x_continuous(limits = c(2010, 2024), breaks=seq(2011,2023,2),labels=seq(2011,2023,2)) +
  labs(x="Year", y="Mean occupancy probability", title='Bananaquit') +
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
        panel.grid.major = element_line(size=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))
occupancy_plot

# save as image 
ggsave(filename = sprintf('output/plot/%s_occupancy_plot.jpg', SPECIES), 
       plot = occupancy_plot, scale = 1, dpi = 'retina', units = 'px', 
       width = 3000, height = 2200)
saveRDS(occupancy_plot, file = sprintf('output/plot/data/%s_occupancy_plot.rds', SPECIES))

##### 8: export some prepared data #### 

saveRDS(umf, file = sprintf('output/data/prepared_data/%s_unmarked_mult_frame.rds', SPECIES))

##### 9: notes ####

