### data exploration ### 
### written by filibert Heim, filibert.heim@posteo.de in April 2024 #### 

# load packages
library(RODBC)
library(reshape)
library(tidyverse)
library(lubridate)
library(unmarked)
library(AICcmodavg) # package for model selection and goodness-of-fit tests
library(MuMIn)
library(terra)
library(scales)
library(ggpubr)
filter<-dplyr::filter
select<-dplyr::select
rename<-dplyr::rename

# load data from Steffen
setwd('C:/Users/filib/Documents/Studium/Bachelorarbeit/R_BachelorThesisMontserrat/GitHub_Steffen_Montserrat/')
load(file = 'MONTSERRAT_ANNUAL_DATA_INPUT2023.RData')
ls()

# load data from database
setwd('C:/Users/filib/Documents/Studium/Bachelorarbeit/Datengrundlage/')
mt_db <- odbcConnectAccess2007('Montserrat_Birds_2023.accdb')
rain_y <- sqlQuery(mt_db, 'SELECT * FROM yearlyRainCov')
siteCov_db <- sqlQuery(mt_db, 'SELECT * FROM SiteCov')
obsCov_db <- sqlQuery(mt_db, 'SELECT * FROM obsCovariates_simple')
obsCov_db$Point <- as.character(obsCov_db$Point)
obsCov_db$Wind <- as.factor(obsCov_db$Wind)
odbcClose(mt_db)

# load spatial data
dsm_mt <- rast(x = 'C:/Users/filib/Documents/Studium/Bachelorarbeit/Datengrundlage/spatial data/digital terrain model - ALOS/ALPSMLC30_N016W063_DSM.tif')
ch_pa <- vect(x = 'C:/Users/filib/Documents/Studium/Bachelorarbeit/Datengrundlage/spatial data/Centre Hills/WDPA_WDOECM_Apr2024_Public_555622117_shp-polygons.shp')


#### check obsCovs 
# prepare data set for checking 
obsCov_test <- obsCov_db %>% filter(Point %in% c(7,10,13,15,16,18,19,25,26,27,28,33,35,36,37,44,45,47,48,50,51,56,57,58,59,60,61,62,68,69,70,71,72,73,74,75,77,81,83,84,85,86,90,93,94,95,97,98,101,104,105,106,110,112,113,114,115,116,119,120,121,123,124,126,129,131,134,137,140,151,152,153,154,155,156,157,158,159,160,161,170,171,172,173,174,175))
# check time - there are mistakes in times, first reasonable time is 5:27am - thus, 5:30
obsCov_test %>% filter(Time < as.POSIXct("1899-12-30 06:00:00")) %>%
  arrange(Time) %>%
  ggplot() + 
  geom_histogram(mapping = aes(x = Time))

# check end time - latest time is 13:30 - thus, 13:30 
obsCov_test %>% filter(Time > as.POSIXct("1899-12-30 12:00:00")) %>%
  arrange(Time) %>%
  ggplot() + 
  geom_histogram(mapping = aes(x = Time))

# check time period - end March until 
obsCov_test %>% 
  mutate(Month = obsCov_test %>% pull(Date) %>% format(format = "%m")) %>% 
  mutate(Day = obsCov_test %>% pull(Date) %>% format(format = "%d")) %>%
  ggplot() + 
  geom_bar(mapping = aes(x = Day)) + 
  facet_wrap(~year + Month, nrow = 3)
ggsave('C:/Users/filib/Documents/Studium/Bachelorarbeit/R_BachelorThesisMontserrat/output/plot/data_exploration/SurveyDatePlot.jpeg', 
       last_plot(), width = 15, height = 5)

obsCov_test %>% # this odd vales are probably caused by different time zones 
  group_by(year) %>% 
  mutate(Period = max(Date) - min(Date)) %>% 
  select(year, Period) %>% 
  unique() %>% 
  ggplot() + 
  geom_point(mapping = aes(x = year, y = Period)) + 
  expand_limits(x = 2011, y = 0) + 
  ylab(label = '# of Days between first and last survey')

(period <- obsCov_test %>%  
  filter(Date < "2016-06-01 EEST" | Date > "2016-08-20 EEST") %>% # remove wrong submitted dates
  group_by(year) %>% 
  mutate(Period = max(Date) - min(Date)) %>% 
  select(year, Period) %>% 
  unique())
mean(round(as.numeric(period$Period))) # mean number of days between the first and the last survey

obsCov_test %>% # plot the different diff times 
  filter(Date < "2016-06-01 EEST" | Date > "2016-08-20 EEST") %>%
  group_by(year) %>% 
  mutate(Period = max(Date) - min(Date)) %>% 
  select(year, Period) %>% 
  unique() %>% 
  ggplot() + 
  geom_point(mapping = aes(x = year, y = Period)) + 
  expand_limits(x = 2011, y = 20)+ 
  ylab(label = '# of Days between first and last survey')
ggsave('C:/Users/filib/Documents/Studium/Bachelorarbeit/R_BachelorThesisMontserrat/output/plot/data_exploration/SurveyPeriodPlott.jpeg', 
       last_plot(), width = 7, height = 5)


# investigate how many counts there have been made for each point and year
missing_surveys <- print(obsCov_test %>% group_by(year, Point) %>% 
  count() %>% 
  filter(n != 3), n = 23) 

length(unique(obsCov_test$Point))

# calculate number of conducted point counts and overall observations
surveys_raw <- dim(obsCov_test %>% filter(year != 2020))[1]
missing_surveys <- missing_surveys %>% 
  mutate(n_missing = 3 - n)
surveys_missing <- sum(missing_surveys$n_missing) # calculate the number of missing surveys
(print(surveys_raw - surveys_missing)*12)

# check for representation of altitudinal range of the points 
names(siteCov)[c(1:9, 17, 18, 21)] <- c('point', 'eastings', 'northings', 'habitat', 
                                        'dbh', 'distance', 'treeheight', 'elevation', 
                                        'slope', 'canopy', 'location', 'alt')
siteCov %>% mutate(point = as.numeric(point)) %>% 
  arrange(point) %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = alt), bins = 15) +
  theme_bw()

# check for representation of altitudinal range of the Centre Hills protected area
plot(dsm_mt)
plot(ch_pa, add = T)

ch_alt <- data.frame(terra::extract(x = dsm_mt, y = ch_pa, ID = F))
(plot_alt_ch <- ggplot() + 
  geom_histogram(data = ch_alt, mapping = aes(x = ALPSMLC30_N016W063_DSM), fill = 'darkred', alpha = 0.4) +
  labs(x = 'Altitude asl', y = 'Frequencies of grid cells', title = 'Altitudinal ranges within the Centre Hills, Montserrat') +
  theme_bw())

(plot_alt_points <- ggplot() + 
    geom_histogram(data = siteCov, mapping = aes(x = alt), fill = 'darkblue', alpha = 0.4) +
    labs(x = 'Altitude asl', y = 'Frequencies of monitoring points', title = 'Altitudinal ranges of monitoring points within the Centre Hills, Montserrat') +
    theme_bw())

ggarrange(plot_alt_ch, plot_alt_points)
ggsave(filename = 'output/plot/data_exploration/AltitudinalRepresentationPointsPlot.jpeg', last_plot())

# check for elevational range of points
alt <- data.frame(max = max(siteCov$alt), min = min(siteCov$alt), diff = max(siteCov$alt - min(siteCov$alt)), type = 'alt') %>% 
  rbind(c(max(siteCov$elevation), min(siteCov$elevation), max(siteCov$elevation - min(siteCov$elevation)), 'elevation'))

  
                   