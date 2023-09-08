library(excel.link)
library(zoo)
library(plyr)
library(dplyr)
library(ggpubr)
library(binom)
library(gridExtra)
library(ggplot2)
library(viridis)

tanz_data_all_2017 <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Tanzania/TZ_ANC_data_region_2014_2017.csv')
data<-tanz_data_all_2017
tanz_data_process <- function(data, level = c('District','Region'), remove_before=NULL){
  if(is.null(remove_before)){
    remove_before <- min(data$yearmon)
  }
  level_call <- ifelse(level=='District','Council','Region')
  
  all <- data.frame(month = as.yearmon(data$yearmon),
                    tested = data$tested,
                    positive = data$positive,
                    site = data[,level_call])
  lt20 <- data.frame(month = as.yearmon(data$yearmon),
                       tested = data$test_LT20,
                       positive = data$pos_LT20,
                     site = data[,level_call])
  ge20 <- data.frame(month = as.yearmon(data$yearmon),
                     tested = data$test_GE20,
                     positive = data$pos_GE20,
                     site = data[,level_call])
  
  split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])
  
  all_list <- split_tibble(all[all$month>=remove_before,],'site')
  lt20_list <- split_tibble(lt20[lt20$month>=remove_before,],'site')
  ge20_list <- split_tibble(ge20[ge20$month>=remove_before,],'site')
  
  return(list(all = all_list, lt20 = lt20_list, ge20 = ge20_list))
}

tanz_data_list <- tanz_data_process(tanz_data_all_2017,'Region')
tanz_data_list_15to17 <- tanz_data_process(tanz_data_all_2017,remove_before = as.yearmon('Jan 2015'),level='Region')
rm(tanz_data_list_14to17)
names(tanz_data_list$lt20[1])
tanz_rainfall <- read_csv("./tanz/processed_inputs/ANC_data_Rainfall_ad1.csv")%>%
  mutate(yearmon=as.yearmon(yearmon))
##Lake Malawi Districts
lake_malawi_raw <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Tanzania/ANC_rainfall_Lake_Malawi_dists.csv')
tanz_lakemal_list_15to17 <- tanz_data_process(lake_malawi_raw,level = 'District', remove_before = as.yearmon('Jan 2015'))


###Updated data March 2023###
#### read in ANC data ##
Full_data_2016<-readxl::read_excel("./tanz/Patrick/Data/ANC_Data_Compiled.xlsx",sheet = '2016_2019')
Full_data_2020<-readxl::read_excel("./tanz/Patrick/Data/ANC_Data_Compiled.xlsx",sheet = '2020_2021')
Full_data_2022<-readxl::read_excel("./tanz/Patrick/Data/ANC_Data_Compiled.xlsx",sheet = '2022')
var.names <- c('region','council','hg','periodid','year','month','first_anc_lt12w','first_anc_gt12w','reattend_anc','reattend_total','tested','positive','llin','ipt2','ipt3','ipt4','hb_measured','hb_lt8.5')
names(Full_data_2016) <- var.names[-4]
names(Full_data_2020) <- var.names[-4]
names(Full_data_2022) <- var.names
Full_data_new <- plyr::rbind.fill(Full_data_2016,Full_data_2020,Full_data_2022)
#Full_data$Region_mb=Full_data$Region
Full_data_new$region[Full_data_new$region=="Songwe Region"]="Mbeya Region"

##Combine old and new data sets
tanz_data_all_2017_region <- read_csv("./tanz/Patrick/processed_inputs/TZ_ANC_data_region_2014_2017.csv")%>%
  select(1:6)%>%
  filter(yearmon <= as.yearmon('Dec 2015')) 
tanz_data_all_16to22_region <- read_csv("./tanz/Patrick/processed_inputs/TZ_ANC_data_region_2016_2022.csv")

tanz_data_all_14to22_region <- bind_rows(tanz_data_all_2017_region,tanz_data_all_16to22_region)%>%
  dplyr::rename(region = Region)

##Combine old and new data sets - district councils
tanz_data_all_2017_district <- read_csv("./tanz/Patrick/processed_inputs/TZ_ANC_data_district_2014_2017.csv")%>%
  select(1:7)%>%
  filter(yearmon <= as.yearmon('Dec 2015')) 
tanz_data_all_16to22_district <- read_csv("./tanz/Patrick/processed_inputs/TZ_ANC_data_district_2016_2022.csv")

tanz_data_all_14to22_district <- bind_rows(tanz_data_all_2017_district,tanz_data_all_16to22_district)

data <- tanz_data_all_14to22_region
level <- 'Region'
remove_before <- as.yearmon('Jan 2015')
tanz_data_process_allonly <- function(data, level = c('District','Region'), remove_before=NULL){
  if(is.null(remove_before)){
    remove_before <- min(data$yearmon)
  }
  level_call <- ifelse(level=='District','council','region')
  
  all <- data.frame(month = as.yearmon(data$yearmon),
                    tested = data$tested,
                    positive = data$positive,
                    site = data[,level_call])
  all <- all %>%
    filter(positive<=tested)
 
  split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])
  
  all_list <- split_tibble(all[all$month>=remove_before,],level_call)
  

  return(all_list)
}

tanz_data_list_15to22 <- tanz_data_process_allonly(tanz_data_all_14to22_region,remove_before = as.yearmon('Jan 2015'),level='Region')
saveRDS(tanz_data_list_15to22,'./tanz/tanz_data_list_15to22.rds')
tanga_only <- tanz_data_all_14to22_district%>%
  filter(Region=='Tanga Region')%>%
  rename(region = Region,
         council = Council)
tanga_data_list_15to22 <- tanz_data_process_allonly(tanga_only,remove_before = as.yearmon('Jan 2015'),level='District')
 