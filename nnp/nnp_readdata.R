library(bayesplot)
library(tidyverse)
library(excel.link)
library(zoo)
library(ggpubr)
library(binom)
library(gridExtra)
library(ggplot2)
library(viridis)

source('shared/addCIs.R')
#######################
########NIGERIA########
#######################
#Read in Nigeria data and rename most used variables
NG_ANC_mother <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/ANC-based surveillance/data_anc_mother_nigeria.xlsx')%>%
  dplyr::rename(primigrav = q1_preg,
         prev_pregs = q2_preg,
         mal_symp = q1_mal,
         rdt = q2_mal) %>%
  dplyr::mutate(month = as.yearmon(month),
         grav = prev_pregs+1,
         ward = toupper(ward))
#Read in Nigeria data and rename most used variables
NG_CS_child_2020 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_child_nigeria_2020.xlsx')%>%
  dplyr::rename(rdt = q85c_result,
         mal_symp = q86a_symptoms) %>%
  mutate(age = ifelse(age==0,age_months/12,age))

NG_CS_hh_2020 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_hh_nigeria_2020.xlsx')%>%
  dplyr::rename(individual_id = hh_id,
         date = start) %>%
  mutate(date = as.Date(date))

#Merge child to hh to get date of survey
NG_CS_2020 <- merge(NG_CS_child_2020,NG_CS_hh_2020,by='submission_id',all.x=TRUE,suffixes=c('.child','.hh'))

##Read in 2021 data
NG_CS_child_2021 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_child_nigeria_2021.xlsx')%>%
  dplyr::rename(rdt = q90c_result_6to59months,
         mal_symp = q91_a_symptoms) %>%
  mutate(age = ifelse(age==0,age_months/12,age))

NG_CS_hh_2021 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_hh_nigeria_2021.xlsx')%>%
  dplyr::rename(individual_id = hh_id,
         date = start) %>%
  mutate(date = as.Date(date))

#Merge child to hh to get date of survey
NG_CS_2021 <- merge(NG_CS_child_2021,NG_CS_hh_2021,by='submission_id',all.x=TRUE,suffixes=c('.child','.hh'))

##Read in 2022 data
NG_CS_child_2022 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_child_nigeria_2022.xlsx')%>%
  dplyr::rename(rdt = q90c_result_6to59months,
         mal_symp = q91a_symptoms) %>%
  mutate(age = ifelse(age==0,age_months/12,age))

NG_CS_hh_2022 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Cross-sectional survey/data_nnp_survey_hh_nigeria_2022.xlsx')%>%
  dplyr::rename(individual_id = hh_id,
         date = start) %>%
  mutate(date = as.Date(date))

#Merge child to hh to get date of survey
NG_CS_2022 <- merge(NG_CS_child_2022,NG_CS_hh_2022,by='submission_id',all.x=TRUE,suffixes=c('.child','.hh'))

#Combine 3 years of data
NG_CS_all <- plyr::rbind.fill(NG_CS_2020,NG_CS_2021,NG_CS_2022)
NG_CS_all$district_code <- as.integer(substr(NG_CS_all$cluster_n,1,1))
NG_CS_all$district_n <- unlist(lapply(NG_CS_all$district_code, FUN = function(x) switch(x,'Ejigbo','Ife North','Asa','Moro')))

#By site and gravidity#
NG_CS_all_grouped_site <- NG_CS_all %>%
  dplyr::rename(site=district_n)%>%
  mutate(month=as.yearmon(date),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA)))
  )%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(site))%>%
  group_by(site,month,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
NG_CS_all_grouped_site <- addCIs(NG_CS_all_grouped_site,NG_CS_all_grouped_site$positive,NG_CS_all_grouped_site$total)
saveRDS(NG_CS_all_grouped_site,'nnp/data/NG_CS_all_grouped_site.rds')
NG_CS_all_grouped_site<-readRDS('nnp/data/NG_CS_all_grouped_site.rds')
#Group by grav
NG_ANC_mother_grouped_sitegrav <- NG_ANC_mother %>%
  dplyr::rename(site = lga) %>%
  dplyr::mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")))%>%
  filter(!is.na(site))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
NG_ANC_mother_grouped_sitegrav <- addCIs(NG_ANC_mother_grouped_sitegrav,NG_ANC_mother_grouped_sitegrav$positive,NG_ANC_mother_grouped_sitegrav$total)
saveRDS(NG_ANC_mother_grouped_sitegrav,'nnp/data/NG_ANC_mother_grouped_sitegrav.rds')
NG_ANC_mother_grouped_sitegrav<-readRDS('nnp/data/NG_ANC_mother_grouped_sitegrav.rds')

##secundigrav
NG_ANC_mother_grouped_sitegrav_sg <- NG_ANC_mother %>%
  dplyr::rename(site = lga) %>%
  dplyr::mutate(grav_cat=cut(grav,breaks=c(0,1,2,Inf),labels=c("Gravidities 1","Gravidities 2","Gravidities 3+")))%>%
  filter(!is.na(site))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
NG_ANC_mother_grouped_sitegrav_sg <- addCIs(NG_ANC_mother_grouped_sitegrav_sg,NG_ANC_mother_grouped_sitegrav_sg$positive,NG_ANC_mother_grouped_sitegrav_sg$total)
saveRDS(NG_ANC_mother_grouped_sitegrav_sg,'nnp/data/NG_ANC_mother_grouped_sitegrav_sg.rds')


#######################
######MOZAMBIQUE#######
#######################
#Read in Mozambique data and rename most used variables
MZ_ANC_mother <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Mozambique/ANC-based surveillance/ANC_Mozambique_2022.09_FINAL.xlsx',
                                   xl.sheet = 'ANC')%>%
  dplyr::rename(primigrav = q4_primagravidae,
                prev_pregs = q5_num_pregnancy,
                mal_symp = q6_mal_symptoms,
                rdt = q7_rdt_result,
                age = q2_age_n) %>%
  dplyr::mutate(month = as.yearmon(as.Date(date_interview_n)),
                grav = ifelse(primigrav=='Yes',1,ifelse(prev_pregs==0,NA,prev_pregs+1)))

#Read in Mozambique data and rename most used variables
MZ_CS_rdt_base <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl180123/Mozambique/Cross-sectional survey/Moz Baseline CSS Datasets 2020.xlsx',
                               xl.sheet = 'RDT')%>%
  mutate(date = as.Date(date))
MZ_CS_rdt_mid <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl180123/Mozambique/Cross-sectional survey/Moz Midline CSS Datasets 2021.xlsx',
                              xl.sheet = 'RDT')%>%
  mutate(date = as.Date(date,"%m/%d/%Y"))
MZ_CS_rdt_end <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl180123/Mozambique/Cross-sectional survey/Moz Endline CSS Datasets 2022.xlsx',
                              xl.sheet = 'RDT')%>%
  mutate(date_interview = as.Date(date_interview))%>%
  dplyr::rename(date = date_interview)
MZ_CS_rdt <- rbind(MZ_CS_rdt_base,MZ_CS_rdt_mid,MZ_CS_rdt_end)

#Group CS data by site and calculate prevalence and CI by month
MZ_CS_all_grouped_site <- MZ_CS_rdt %>%
  dplyr::rename(site = district) %>%
  mutate(month=as.yearmon(date),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA)))
  )%>%
  filter(!is.na(rdt)) %>%
  filter(site %in% c('Changara','Chemba','Guro'))%>%
  mutate(month=as.yearmon(ifelse(month=='Sep 2020',as.yearmon('Oct 2020'),as.yearmon(month))))%>% ##Group Sep and Oct 2020 together for the pre-intervention CX survey
  group_by(site,month,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
MZ_CS_all_grouped_site <- addCIs(MZ_CS_all_grouped_site,MZ_CS_all_grouped_site$positive,MZ_CS_all_grouped_site$total)
saveRDS(MZ_CS_all_grouped_site,'nnp/data/MZ_CS_all_grouped_site_180123.rds')
MZ_CS_all_grouped_site <- readRDS('nnp/data/MZ_CS_all_grouped_site_180123.rds')

#Group ANC data by grav and site
MZ_ANC_mother_grouped_sitegrav <- MZ_ANC_mother %>%
  rename(site = district_n) %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA))))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age>=12&age<50) %>%
  filter(!is.na(grav))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
MZ_ANC_mother_grouped_sitegrav <- addCIs(MZ_ANC_mother_grouped_sitegrav,MZ_ANC_mother_grouped_sitegrav$positive,MZ_ANC_mother_grouped_sitegrav$total)
saveRDS(MZ_ANC_mother_grouped_sitegrav,'nnp/data/MZ_ANC_mother_grouped_sitegrav_0822.rds')
MZ_ANC_mother_grouped_sitegrav <- readRDS('nnp/data/MZ_ANC_mother_grouped_sitegrav_0822.rds')

### secundigrav
MZ_ANC_mother_grouped_sitegrav_sg <- MZ_ANC_mother %>%
  dplyr::rename(site = district_n) %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,2,Inf),labels=c("Gravidities 1","Gravidities 2","Gravidities 3+")),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA))))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age>=12&age<50) %>%
  filter(!is.na(grav))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
MZ_ANC_mother_grouped_sitegrav_sg <- addCIs(MZ_ANC_mother_grouped_sitegrav_sg,MZ_ANC_mother_grouped_sitegrav_sg$positive,MZ_ANC_mother_grouped_sitegrav_sg$total)
saveRDS(MZ_ANC_mother_grouped_sitegrav_sg,'nnp/data/MZ_ANC_mother_grouped_sitegrav_0822_sg.rds')

#########################
######Burkina Faso#######
#########################
#Read in Burkina Faso data and rename most used variables
BF_ANC_mother <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Burkina Faso/ANC-based surveillance/ANC_BF_2022.12.28_FINAL.xlsx',
                              xl.sheet = 'ANC_woman')%>%
  dplyr::rename(primigrav = Q04_first_preg,
                prev_pregs = Q05_nbr_prior_preg,
                mal_symp = Q06_mal_symp,
                rdt = RDT,
                age = Age) %>%
  dplyr::mutate(month = as.yearmon(as.Date(Date)),
                grav = ifelse(primigrav=='Oui',1,prev_pregs+1))
#Read in Burkina Faso data and rename most used variables
BF_CS_child_2019 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl180123/Burkina Faso/Cross-sectional survey/NNET_CSS_2019.xls',
                                 xl.sheet = 'Child_data')%>%
  dplyr::rename(rdt = Result_RDT) %>%
  dplyr::mutate(month = as.yearmon(as.Date('2019-07-01')))

BF_CS_child_2020 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl180123/Burkina Faso/Cross-sectional survey/NNET_CSS_2020.xls',
                                 xl.sheet = 'Child_data')%>%
  dplyr::rename(rdt = Result_RDT) %>%
  dplyr::mutate(month = as.yearmon(as.Date('2020-06-01')))

BF_CS_child_2021 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl180123/Burkina Faso/Cross-sectional survey/NNET_CSS_2021.xlsx',
                                 xl.sheet = 'Child_data')%>%
  dplyr::rename(rdt = Result_RDT) %>%
  dplyr::mutate(month = as.yearmon(as.Date('2021-06-01')))

BF_CS_child_2022 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl180123/Burkina Faso/Cross-sectional survey/NNET_CSS_2022.xlsx',
                                 xl.sheet = 'HH_data')%>%
  dplyr::rename(rdt = Q87_tdr,
                District = district) %>%
  dplyr::mutate(month = as.yearmon(as.Date('2022-06-01')))

BF_CS_child_all <- rbind.fill(BF_CS_child_2019,BF_CS_child_2020,BF_CS_child_2021,BF_CS_child_2022)

#Group CS data and calculate prevalence and CI by month
BF_CS_all_grouped_site <- BF_CS_child_all %>%
  rename(site = District) %>%
  mutate(
    rdt=as.numeric(ifelse(rdt=='Négatif',0,1))
  )%>%
  filter(!is.na(rdt)) %>%
  group_by(site,month,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
BF_CS_all_grouped_site <- addCIs(BF_CS_all_grouped_site,BF_CS_all_grouped_site$positive,BF_CS_all_grouped_site$total)
saveRDS(BF_CS_all_grouped_site,'nnp/data/BF_CS_all_grouped_site_0822.rds')
BF_CS_all_grouped_site <-readRDS('nnp/data/BF_CS_all_grouped_site_0822.rds')

#Group ANC data by site
#Group by grav
BF_ANC_mother_grouped_sitegrav <- BF_ANC_mother %>%
  dplyr::rename(site = District) %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")),
         rdt=as.numeric(ifelse(rdt=='Positif',1,0)))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age!=88) %>%
  filter(!is.na(grav))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
BF_ANC_mother_grouped_sitegrav <- addCIs(BF_ANC_mother_grouped_sitegrav,BF_ANC_mother_grouped_sitegrav$positive,BF_ANC_mother_grouped_sitegrav$total)
saveRDS(BF_ANC_mother_grouped_sitegrav,'nnp/data/BF_ANC_mother_grouped_sitegrav.rds')
BF_ANC_mother_grouped_sitegrav <- readRDS('nnp/data/BF_ANC_mother_grouped_sitegrav.rds')

##secundigrav
BF_ANC_mother_grouped_sitegrav_sg <- BF_ANC_mother %>%
  dplyr::rename(site = District) %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,2,Inf),labels=c("Gravidities 1","Gravidities 2","Gravidities 3+")),
         rdt=as.numeric(ifelse(rdt=='Positif',1,0)))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age!=88) %>%
  filter(!is.na(grav))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
BF_ANC_mother_grouped_sitegrav_sg <- addCIs(BF_ANC_mother_grouped_sitegrav_sg,BF_ANC_mother_grouped_sitegrav_sg$positive,BF_ANC_mother_grouped_sitegrav_sg$total)
saveRDS(BF_ANC_mother_grouped_sitegrav_sg,'nnp/data/BF_ANC_mother_grouped_sitegrav_sg.rds')

#########################
#########Zambia##########
#########################
#Read in Zambia data and rename most used variables
ZM_ANC_mother <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/NNP_zambia_dl210823/ANC Surveillance/ANC Surveillance_clean_attendees.csv') %>%
  # exclude those refusing malaria portion
  filter(consent_granted %in% c("Malaria portion only", "Full study")) %>%
  # merge aggregate treatment seeking data per mother
  mutate(
    month = as.yearmon(as.Date(survey_date)),
    
    # RDT testing and result among women of reproductive age (same as rdt for ANC attendees); 
    # added for dataset consistency 
    rdt = as.numeric(na_if(standard_rdt_result, "Not Valid")=="Pf Positive"),
    primigrav = first_pregnancy,
    prev_pregs = num_of_prior_pregnancies,
    rdt = RDT,
    age = age_in_years,
    grav = ifelse(primigrav=='Yes',1,prev_pregs+1),
    # categorize attendee residency
    # within Zambia
    res_zam = as.numeric(live_in_this_hf_catchment=="Yes" | 
                           residence %in% c("Another Part Of Chadiza District",
                                            "In Another District In Eastern Province",
                                            "In Zambia, But Outside Of Eastern Province")),
    res_zam = replace(res_zam, live_in_this_hf_catchment=="Refused", NA),
    # whether near border (in Chadiza)
    res_zam_border = as.numeric(live_in_this_hf_catchment=="Yes" & 
                                  health_facility %in% c("Chanida Health Post", "Kapirimpika Health Post",
                                                         "Mkumbuzi Rural Health Centre", "Mtaya Health Post",
                                                         "Miti Rural Health Centre")),
    res_zam_border = replace(res_zam_border, res_zam=="No" | is.na(res_zam), NA),
    # together in one variable
    res_cat = case_when(res_zam==T & res_zam_border==F            ~ "Zambia, non-border",
                        res_zam==T & res_zam_border==T            ~ "Zambia, border",
                        residence=="Mozambique"                   ~ "Mozambique",
                        residence %in% c("Malawi", "Other Area")  ~ "Other",
                        TRUE                                      ~ NA_character_),
    
  ) %>%
  # drop attendees with uncategorizable residence
  filter(!is.na(res_cat)) %>%
  # drop obs before 2020-01-01; error dates before start of data collection
  filter(month >= "2020-01-01")

#Read in Zambia data and rename most used variables

#Group CS data and calculate prevalence and CI by month
# BF_CS_all_grouped_site <- BF_CS_child_all %>%
#   rename(site = District) %>%
#   mutate(
#     rdt=as.numeric(ifelse(rdt=='Négatif',0,1))
#   )%>%
#   filter(!is.na(rdt)) %>%
#   group_by(site,month,.drop=FALSE)%>%
#   dplyr::summarise(positive=sum(rdt),total=n())
# BF_CS_all_grouped_site <- addCIs(BF_CS_all_grouped_site,BF_CS_all_grouped_site$positive,BF_CS_all_grouped_site$total)
# saveRDS(BF_CS_all_grouped_site,'nnp/data/BF_CS_all_grouped_site_0822.rds')
# BF_CS_all_grouped_site <-readRDS('nnp/data/BF_CS_all_grouped_site_0822.rds')

#Group ANC data by site
#Group by grav
table(ZM_ANC_mother$outreach_zone)
table(ZM_ANC_mother$age)
table(ZM_ANC_mother$age_in_years)
table(ZM_ANC_mother$age_in_months)

ZM_ANC_mother_grouped_sitegrav <- ZM_ANC_mother %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")),
         month = ifelse(month==as.yearmon('Feb 2020'),as.yearmon('Mar 2020'),month))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(grav))%>%
  group_by(month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
ZM_ANC_mother_grouped_sitegrav <- addCIs(ZM_ANC_mother_grouped_sitegrav,ZM_ANC_mother_grouped_sitegrav$positive,ZM_ANC_mother_grouped_sitegrav$total)
saveRDS(ZM_ANC_mother_grouped_sitegrav,'nnp/data/ZM_ANC_mother_grouped_sitegrav.rds')
ZM_ANC_mother_grouped_sitegrav <- readRDS('nnp/data/ZM_ANC_mother_grouped_sitegrav.rds')

##secundigrav
BF_ANC_mother_grouped_sitegrav_sg <- BF_ANC_mother %>%
  dplyr::rename(site = District) %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,2,Inf),labels=c("Gravidities 1","Gravidities 2","Gravidities 3+")),
         rdt=as.numeric(ifelse(rdt=='Positif',1,0)))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age!=88) %>%
  filter(!is.na(grav))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
BF_ANC_mother_grouped_sitegrav_sg <- addCIs(BF_ANC_mother_grouped_sitegrav_sg,BF_ANC_mother_grouped_sitegrav_sg$positive,BF_ANC_mother_grouped_sitegrav_sg$total)
saveRDS(BF_ANC_mother_grouped_sitegrav_sg,'nnp/data/BF_ANC_mother_grouped_sitegrav_sg.rds')

########################################
####Create Plots for Each Country#######
########################################
windows(10,8)
#NIGERIA#
#Create data tables to annotate figure
NG_iv <- data.frame(site = c('Ejigbo','Asa','Moro','Ife North'),
                    month = as.yearmon(rep(as.Date('2020-11-01'),4)))
NG_labs <- c('Ejigbo - Standard','Asa - IG2','Moro - RG','Ife North - PBO')
names(NG_labs) <- c('Ejigbo','Asa','Moro','Ife North')
#Grav colored, facetted by site
ggplot(data=NG_ANC_mother_grouped_sitegrav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  annotate("rect", xmin = as.yearmon('Oct 2020'), xmax = as.yearmon('Nov 2020'), ymin = 0, ymax = 1,alpha = .1,fill = "#0072B2")+
  annotate("rect", xmin = as.yearmon('Jul 2021'), xmax = as.yearmon('Nov 2021'), ymin = 0, ymax = 1,alpha = .1,fill = "#0072B2")+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  facet_wrap(vars(site),nrow=2,ncol=2,labeller = labeller(site = NG_labs))+
  geom_point(data=NG_CS_all_grouped_site,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=NG_CS_all_grouped_site,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")+
  geom_vline(data=NG_iv,aes(xintercept=month),color='darkgrey',size=2,alpha = 0.5)+
  theme(legend.position="bottom",
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1.2, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 12, face = "bold"))


#MOZAMBIQUE#
#Create data tables to annotate figure
MZ_iv <- data.frame(site = c('Chemba','Guro','Changara'),
                    month = as.yearmon(rep(as.Date('2020-11-01'),3)))
MZ_labs <- c('Chemba - Standard','Guro - IG2','Changara - PBO')
names(MZ_labs) <- c('Chemba','Guro','Changara')

#Grav colored
ggplot(data=MZ_ANC_mother_grouped_sitegrav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  annotate("rect", xmin = as.yearmon('Jan 2021'), xmax = as.yearmon('Jun 2021'), ymin = 0, ymax = 1,alpha = .1,fill = "#0072B2")+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  geom_point(data=MZ_CS_all_grouped_site,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=MZ_CS_all_grouped_site,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")+
  facet_wrap(vars(site),nrow=2,ncol=2,labeller = labeller(site = MZ_labs))+
  geom_vline(data=MZ_iv,aes(xintercept=month),color='darkgrey',size=2,alpha = 0.5)+
  theme(legend.position="bottom",
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1.2, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 12, face = "bold"))

#BURKINA FASO#
#Create data tables to annotate figure
BF_iv <- data.frame(site = c('Gaoua','Banfora','Orodara'),
                    month = as.yearmon(c(as.Date('2019-08-01'),as.Date('2019-10-01'),as.Date('2019-06-01'))))
BF_labs <- c('Gaoua - Standard','Banfora - IG2','Orodara - PBO')
names(BF_labs) <- c('Gaoua','Banfora','Orodara')

ggplot(data=BF_ANC_mother_grouped_sitegrav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  annotate("rect", xmin = as.yearmon('Jun 2019'), xmax = as.yearmon('Oct 2019'), ymin = 0, ymax = 1,alpha = .1,fill = "#0072B2")+
  annotate("rect", xmin = as.yearmon('Jun 2020'), xmax = as.yearmon('Oct 2020'), ymin = 0, ymax = 1,alpha = .1,fill = "#0072B2")+
  annotate("rect", xmin = as.yearmon('Jun 2021'), xmax = as.yearmon('Oct 2021'), ymin = 0, ymax = 1,alpha = .1,fill = "#0072B2")+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  geom_point(data=BF_CS_all_grouped_site,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=BF_CS_all_grouped_site,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")+
  facet_wrap(vars(site),nrow=2,ncol=2,labeller = labeller(site = BF_labs))+
  geom_vline(data=BF_iv,aes(xintercept=month),color='darkgrey',size=2,alpha = 0.5)+
  theme(legend.position="bottom",
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1.2, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 12, face = "bold"))

ggplot(data=ZM_ANC_mother_grouped_sitegrav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom",
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1.2, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 12, face = "bold"))

#######################################
######Process data for pMCMC###########
#######################################
##Primigrav data sets
##Nigeria
NG_anc <- NG_ANC_mother_grouped_sitegrav
#readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/NG_ANC_mother_grouped_sitegrav.rds')

NG_pg_asa <- NG_anc[NG_anc$site=='Asa'&NG_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_pg_asa,'nnp/data/data_raw_NG_pg_asa.RDS')

NG_pg_ejigbo <- NG_anc[NG_anc$site=='Ejigbo'&NG_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_pg_ejigbo,'nnp/data/data_raw_NG_pg_ejigbo.RDS')

NG_pg_ifenorth <- NG_anc[NG_anc$site=='Ife North'&NG_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_pg_ifenorth,'nnp/data/data_raw_NG_pg_ifenorth.RDS')

NG_pg_moro <- NG_anc[NG_anc$site=='Moro'&NG_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_pg_moro,'nnp/data/data_raw_NG_pg_moro.RDS')

##Burkina Faso
BF_anc <- BF_ANC_mother_grouped_sitegrav
#readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/BF_ANC_mother_grouped_sitegrav.rds')

BF_pg_banfora <- BF_anc[BF_anc$site=='Banfora'&BF_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_pg_banfora,'nnp/data/data_raw_BF_pg_banfora.RDS')


BF_pg_gaoua <- BF_anc[BF_anc$site=='Gaoua'&BF_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_pg_gaoua,'nnp/data/data_raw_BF_pg_gaoua.RDS')

BF_pg_orodara <- BF_anc[BF_anc$site=='Orodara'&BF_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_pg_orodara,'nnp/data/data_raw_BF_pg_orodara.RDS')

##Mozambique
MZ_anc <- MZ_ANC_mother_grouped_sitegrav
#readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/MZ_ANC_mother_grouped_sitegrav.rds')

MZ_pg_changara <- MZ_anc[MZ_anc$site=='Changara'&MZ_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_pg_changara,'nnp/data/data_raw_MZ_pg_changara.RDS')

MZ_pg_chemba <- MZ_anc[MZ_anc$site=='Chemba'&MZ_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_pg_chemba,'nnp/data/data_raw_MZ_pg_chemba.RDS')

MZ_pg_guro <- MZ_anc[MZ_anc$site=='Guro'&MZ_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_pg_guro,'nnp/data/data_raw_MZ_pg_guro.RDS')

##Grav 2-3 data sets
##Nigeria
NG_anc <- NG_ANC_mother_grouped_sitegrav
#readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/NG_ANC_mother_grouped_sitegrav.rds')

NG_g23_asa <- NG_anc[NG_anc$site=='Asa'&NG_anc$grav_cat=='Gravidities 2-3',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_g23_asa,'nnp/data/data_raw_NG_g23_asa.RDS')

NG_g23_ejigbo <- NG_anc[NG_anc$site=='Ejigbo'&NG_anc$grav_cat=='Gravidities 2-3',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_g23_ejigbo,'nnp/data/data_raw_NG_g23_ejigbo.RDS')

NG_g23_ifenorth <- NG_anc[NG_anc$site=='Ife North'&NG_anc$grav_cat=='Gravidities 2-3',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_g23_ifenorth,'nnp/data/data_raw_NG_g23_ifenorth.RDS')

NG_g23_moro <- NG_anc[NG_anc$site=='Moro'&NG_anc$grav_cat=='Gravidities 2-3',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_g23_moro,'nnp/data/data_raw_NG_g23_moro.RDS')

##Burkina Faso
BF_anc <- BF_ANC_mother_grouped_sitegrav
#readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/BF_ANC_mother_grouped_sitegrav.rds')

BF_g23_banfora <- BF_anc[BF_anc$site=='Banfora'&BF_anc$grav_cat=='Gravidities 2-3',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_g23_banfora,'nnp/data/data_raw_BF_g23_banfora.RDS')


BF_g23_gaoua <- BF_anc[BF_anc$site=='Gaoua'&BF_anc$grav_cat=='Gravidities 2-3',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_g23_gaoua,'nnp/data/data_raw_BF_g23_gaoua.RDS')

BF_g23_orodara <- BF_anc[BF_anc$site=='Orodara'&BF_anc$grav_cat=='Gravidities 2-3',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_g23_orodara,'nnp/data/data_raw_BF_g23_orodara.RDS')

##Mozambique
MZ_anc <- MZ_ANC_mother_grouped_sitegrav
#readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/MZ_ANC_mother_grouped_sitegrav.rds')

MZ_g23_changara <- MZ_anc[MZ_anc$site=='Changara'&MZ_anc$grav_cat=='Gravidities 2-3',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_g23_changara,'nnp/data/data_raw_MZ_g23_changara.RDS')

MZ_g23_chemba <- MZ_anc[MZ_anc$site=='Chemba'&MZ_anc$grav_cat=='Gravidities 2-3',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_g23_chemba,'nnp/data/data_raw_MZ_g23_chemba.RDS')

MZ_g23_guro <- MZ_anc[MZ_anc$site=='Guro'&MZ_anc$grav_cat=='Gravidities 2-3',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_g23_guro,'nnp/data/data_raw_MZ_g23_guro.RDS')

##Grav 4+ data sets
##Nigeria
NG_anc <- NG_ANC_mother_grouped_sitegrav
#readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/NG_ANC_mother_grouped_sitegrav.rds')

NG_g4_asa <- NG_anc[NG_anc$site=='Asa'&NG_anc$grav_cat=='Gravidities 4+',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_g4_asa,'nnp/data/data_raw_NG_g4_asa.RDS')

NG_g4_ejigbo <- NG_anc[NG_anc$site=='Ejigbo'&NG_anc$grav_cat=='Gravidities 4+',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_g4_ejigbo,'nnp/data/data_raw_NG_g4_ejigbo.RDS')

NG_g4_ifenorth <- NG_anc[NG_anc$site=='Ife North'&NG_anc$grav_cat=='Gravidities 4+',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_g4_ifenorth,'nnp/data/data_raw_NG_g4_ifenorth.RDS')

NG_g4_moro <- NG_anc[NG_anc$site=='Moro'&NG_anc$grav_cat=='Gravidities 4+',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_g4_moro,'nnp/data/data_raw_NG_g4_moro.RDS')

##Burkina Faso
BF_anc <- BF_ANC_mother_grouped_sitegrav
#readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/BF_ANC_mother_grouped_sitegrav.rds')

BF_g4_banfora <- BF_anc[BF_anc$site=='Banfora'&BF_anc$grav_cat=='Gravidities 4+',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_g4_banfora,'nnp/data/data_raw_BF_g4_banfora.RDS')


BF_g4_gaoua <- BF_anc[BF_anc$site=='Gaoua'&BF_anc$grav_cat=='Gravidities 4+',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_g4_gaoua,'nnp/data/data_raw_BF_g4_gaoua.RDS')

BF_g4_orodara <- BF_anc[BF_anc$site=='Orodara'&BF_anc$grav_cat=='Gravidities 4+',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_g4_orodara,'nnp/data/data_raw_BF_g4_orodara.RDS')

##Mozambique
MZ_anc <- MZ_ANC_mother_grouped_sitegrav
#readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/MZ_ANC_mother_grouped_sitegrav.rds')

MZ_g4_changara <- MZ_anc[MZ_anc$site=='Changara'&MZ_anc$grav_cat=='Gravidities 4+',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_g4_changara,'nnp/data/data_raw_MZ_g4_changara.RDS')

MZ_g4_chemba <- MZ_anc[MZ_anc$site=='Chemba'&MZ_anc$grav_cat=='Gravidities 4+',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_g4_chemba,'nnp/data/data_raw_MZ_g4_chemba.RDS')

MZ_g4_guro <- MZ_anc[MZ_anc$site=='Guro'&MZ_anc$grav_cat=='Gravidities 4+',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_g4_guro,'nnp/data/data_raw_MZ_g4_guro.RDS')

##Multigrav
NG_mg_asa <- NG_anc[NG_anc$site=='Asa'&NG_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_mg_asa,'nnp/data/data_raw_NG_mg_asa.RDS')

NG_mg_ejigbo <- NG_anc[NG_anc$site=='Ejigbo'&NG_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_mg_ejigbo,'nnp/data/data_raw_NG_mg_ejigbo.RDS')

NG_mg_ifenorth <- NG_anc[NG_anc$site=='Ife North'&NG_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_mg_ifenorth,'nnp/data/data_raw_NG_mg_ifenorth.RDS')

NG_mg_moro <- NG_anc[NG_anc$site=='Moro'&NG_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_mg_moro,'nnp/data/data_raw_NG_mg_moro.RDS')

##Burkina Faso
BF_mg_banfora <- BF_anc[BF_anc$site=='Banfora'&BF_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_mg_banfora,'nnp/data/data_raw_BF_mg_banfora.RDS')

BF_mg_gaoua <- BF_anc[BF_anc$site=='Gaoua'&BF_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_mg_gaoua,'nnp/data/data_raw_BF_mg_gaoua.RDS')

BF_mg_orodara <- BF_anc[BF_anc$site=='Orodara'&BF_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_mg_orodara,'nnp/data/data_raw_BF_mg_orodara.RDS')

##Mozambique
MZ_mg_changara <- MZ_anc[MZ_anc$site=='Changara'&MZ_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_mg_changara,'nnp/data/data_raw_MZ_mg_changara.RDS')

MZ_mg_chemba <- MZ_anc[MZ_anc$site=='Chemba'&MZ_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_mg_chemba,'nnp/data/data_raw_MZ_mg_chemba.RDS')

MZ_mg_guro <- MZ_anc[MZ_anc$site=='Guro'&MZ_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_mg_guro,'nnp/data/data_raw_MZ_mg_guro.RDS')

##MAll gravidities
NG_all_asa <- NG_anc[NG_anc$site=='Asa',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_all_asa,'nnp/data/data_raw_NG_all_asa.RDS')

NG_all_ejigbo <- NG_anc[NG_anc$site=='Ejigbo',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_all_ejigbo,'nnp/data/data_raw_NG_all_ejigbo.RDS')

NG_all_ifenorth <- NG_anc[NG_anc$site=='Ife North',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_all_ifenorth,'nnp/data/data_raw_NG_all_ifenorth.RDS')

NG_all_moro <- NG_anc[NG_anc$site=='Moro',] %>%
  group_by(month)%>%
  dplyr::summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_all_moro,'nnp/data/data_raw_NG_all_moro.RDS')

##Burkina Faso
BF_all_banfora <- BF_anc[BF_anc$site=='Banfora',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_all_banfora,'nnp/data/data_raw_BF_all_banfora.RDS')

BF_all_gaoua <- BF_anc[BF_anc$site=='Gaoua',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_all_gaoua,'nnp/data/data_raw_BF_all_gaoua.RDS')

BF_all_orodara <- BF_anc[BF_anc$site=='Orodara',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_all_orodara,'nnp/data/data_raw_BF_all_orodara.RDS')

##Mozambique
MZ_all_changara <- MZ_anc[MZ_anc$site=='Changara',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_all_changara,'nnp/data/data_raw_MZ_all_changara.RDS')

MZ_all_chemba <- MZ_anc[MZ_anc$site=='Chemba',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_all_chemba,'nnp/data/data_raw_MZ_all_chemba.RDS')

MZ_all_guro <- MZ_anc[MZ_anc$site=='Guro',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_all_guro,'nnp/data/data_raw_MZ_all_guro.RDS')

table(MZ_ANC_mother$month,MZ_ANC_mother$mal_symp,MZ_ANC_mother$district_n)
table(NG_ANC_mother$month,NG_ANC_mother$mal_symp,NG_ANC_mother$lga)
prop.table(table(NG_ANC_mother$mal_symp))
table(BF_ANC_mother$month,BF_ANC_mother$mal_symp,BF_ANC_mother$District)
prop.table(table(BF_ANC_mother$mal_symp))


##Summary table
all_ANC_mother_grouped_sitegrav <- bind_rows(BF_ANC_mother_grouped_sitegrav,MZ_ANC_mother_grouped_sitegrav,NG_ANC_mother_grouped_sitegrav)
all_ANC_mother_grouped_sitegrav%>%
  group_by(site)%>%
  summarise(start=min(month),
            end=max(month),
            total=sum(total),
            tot.pos = sum(positive),
            prev = round(sum(positive)/sum(total) * 100,1))

total_prev <- all_ANC_mother_grouped_sitegrav%>%
  group_by()%>%
  dplyr::summarise(positive = sum(positive),
                   total = sum(total))
total_prev_cis <- addCIs(total_prev,total_prev$positive,total_prev$total)
