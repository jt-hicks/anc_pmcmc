library(foreign)
library(readstata13)
library(sas7bdat)
library(ggplot2)
library(dplyr)
library(zoo)
library(gridExtra)
library(ggpubr)
library(viridis)

##Read in data and select relevant columns
W_africa_data<-as_tibble(read.csv("trial/Data/1_main_study_file - version for PW and Sequoia.csv"))
WA_select_var<-W_africa_data%>%select(studyid,country,site,age,gravidity,fundal_height,
                                      intervention,rdt_treat0,screen0_date,v0_result,pcr_result)


##Summarise number of rdt tests and positives by month, gravidity, and location
WA_grouped<-WA_select_var%>%
  dplyr::filter(!is.na(gravidity))%>%
  dplyr::mutate(month=as.yearmon(screen0_date,"%d/%m/%Y"))%>%
  group_by(site,month,gravidity)%>%
  dplyr::summarise(positive=sum(v0_result=="positive"),tested=sum(v0_result%in%c("positive","negative")))%>%
  ungroup()%>%
  dplyr::mutate(site = ifelse(site=='Burkina','Burkina Faso',site),
         test_type = 'Slide')%>%
  filter(!(site=='Gambia'&month==as.yearmon('Jun 2010')))
WA_grouped %>%
  group_by(site)%>%
  summarise(start=min(month),
            end=max(month),
            total=sum(tested),
            tot.pos = sum(positive),
            prev = round(sum(positive)/sum(tested) * 100,1))
WA_grouped_rdt<-WA_select_var%>%
  dplyr::filter(!is.na(gravidity))%>%
  mutate(month=as.yearmon(screen0_date,"%d/%m/%Y"))%>%
  group_by(site,month,gravidity)%>%
  summarise(positive=sum(rdt_treat0=="positive"),tested=sum(rdt_treat0%in%c("positive","negative")))%>%
  ungroup()%>%
  mutate(site = ifelse(site=='Burkina','Burkina Faso',site),
         test_type = 'RDT')
WA_grouped_pcr<-WA_select_var%>%
  dplyr::filter(!is.na(gravidity))%>%
  mutate(month=as.yearmon(screen0_date,"%d/%m/%Y"))%>%
  group_by(site,month,gravidity)%>%
  summarise(positive=sum(pcr_result=="positive"),tested=sum(pcr_result%in%c("positive","negative")))%>%
  ungroup()%>%
  mutate(site = ifelse(site=='Burkina','Burkina Faso',site),
         test_type = 'PCR')
WA_all_tests <- rbind(WA_grouped,WA_grouped_rdt,WA_grouped_pcr) %>%
  mutate(test_type = factor(test_type,levels=c('Slide','RDT','PCR')),
         grav_cat = factor(ifelse(gravidity==1,'pg','sg'),levels=c('pg','sg'),labels=c('Primigrav','Secundigrav')),
         prev = positive/tested,
         site = factor(site,levels = c('Ghana', 'Burkina Faso', 'Mali', 'Gambia')))
WA_all_tests <- addCIs(WA_all_tests,WA_all_tests$positive,WA_all_tests$tested)

##Split dataframes into two lists by site, one for each grav level
WA_grouped_pg <- WA_grouped%>%
  filter(gravidity==1)%>%
  group_by(site)
WA_pg_data_list <- split(WA_grouped_pg, f = WA_grouped_pg$site)
WA_grouped_sg <- WA_grouped%>%
  filter(gravidity==2)%>%
  group_by(site)
WA_sg_data_list <- split(WA_grouped_sg, f = WA_grouped_sg$site)

saveRDS(WA_pg_data_list,'trial/Data/WA_pg_data_list.rds')
saveRDS(WA_sg_data_list,'trial/Data/WA_sg_data_list.rds')

##Summarise number of rdt tests and positives by month,  and location
WA_grouped_all<-WA_select_var%>%
  dplyr::filter(!is.na(gravidity))%>%
  mutate(month=as.yearmon(screen0_date,"%d/%m/%Y"))%>%
  group_by(site,month)%>%
  summarise(positive=sum(v0_result=="positive"),tested=sum(v0_result%in%c("positive","negative")))%>%
  ungroup()%>%
  mutate(site = ifelse(site=='Burkina','Burkina Faso',site))

WA_all_data_list <- split(WA_grouped_all, f = WA_grouped_all$site)
saveRDS(WA_all_data_list,'trial/Data/WA_all_data_list.rds')

##########Malawi###############
library(haven)
ISTpMalawi_PWalker2 <- read_dta("trial/Data/ISTpMalawi-PWalker2.dta")

malawi_select_var<-ISTpMalawi_PWalker2%>%select(IDNumber,site,dtbook,grav,pcr_bk,smr_bk,rdt_bk,
                                                anymal_bk)
table(malawi_select_var$site)
table(malawi_select_var$smr_bk)
table(is.na(malawi_select_var$grav))
str(malawi_select_var)
##Summarise number of rdt tests and positives by month, gravidity, and location
malawi_grouped<-malawi_select_var%>%
  dplyr::filter(!is.na(grav)&!is.na(smr_bk))%>%
  mutate(month=as.yearmon(dtbook),
         grav_cat = ifelse(grav==1,'pg','mg'))%>%
  group_by(month,grav_cat)%>%
  summarise(positive=sum(smr_bk==1),tested=sum(smr_bk%in%c(1,0)))%>%
  mutate(site='Malawi')%>%
  ungroup()

WA_grouped_mod <- WA_grouped %>%
  mutate(grav_cat = ifelse(gravidity==1,'pg','sg'))%>%
  select(site,month,grav_cat,positive,tested)

##W Kenya
anc_wkenya <- haven::read_sas('./trial/Data/book_screen_ses_ms.sas7bdat')
table(anc_wkenya$RDT)
anc_wkenya_select_var<-anc_wkenya%>%select(IDNumber,datevisit,Primi,Parity,gravidity,RDT,patent)
table(anc_wkenya_select_var$gravidity,anc_wkenya_select_var$Parity,useNA='always')
table(anc_wkenya_select_var$Primi,useNA='always')
names(anc_wkenya)
wkenya_group <- anc_wkenya%>%
  dplyr::filter(!is.na(Primi)&!is.na(patent))%>%
  mutate(month=as.yearmon(datevisit),
         grav_cat = ifelse(Primi=='Y','pg','mg'))%>%
  group_by(month,grav_cat)%>%
  summarise(positive=sum(patent==1),tested=sum(patent%in%c(0,1)))%>%
  mutate(site='Kenya')%>%
  ungroup()%>%
  mutate(grav_cat=factor(grav_cat,levels=c('pg','sg','mg')),
         grav_dich=factor(ifelse(grav_cat=='pg','pg','mg'),levels=c('pg','mg')),
         prev = positive/tested)

##Split dataframes into two lists by site, one for each grav level
EA_grouped <- plyr::rbind.fill(malawi_grouped,wkenya_group) %>%
  select(month,grav_cat,positive,tested,site)
EA_grouped %>%
  group_by(site)%>%
  summarise(start=min(month),
            end=max(month),
            total=sum(tested),
            tot.pos = sum(positive),
            prev = round(sum(positive)/sum(tested) * 100,1))
EA_grouped_pg <- EA_grouped%>%
  filter(grav_cat=='pg')%>%
  group_by(site)
EA_pg_data_list <- split(EA_grouped_pg, f = EA_grouped_pg$site)
EA_grouped_mg <- EA_grouped%>%
  filter(grav_cat=='mg')%>%
  group_by(site)
EA_mg_data_list <- split(EA_grouped_mg, f = EA_grouped_mg$site)

saveRDS(EA_pg_data_list,'trial/Data/EA_pg_data_list.rds')
saveRDS(EA_mg_data_list,'trial/Data/EA_mg_data_list.rds')


all_prev <- plyr::rbind.fill(WA_grouped_mod,malawi_grouped,wkenya_group)%>%
  mutate(grav_cat=factor(grav_cat,levels=c('pg','sg','mg')),
         grav_dich=factor(ifelse(grav_cat=='pg','pg','mg'),levels=c('pg','mg'),labels=c('Primigrav','Multigrav')),
         prev = positive/tested,
         site = factor(site,levels = c('Ghana', 'Burkina Faso', 'Mali', 'Gambia','Malawi','Kenya')))
all_prev_cis <- addCIs(all_prev,all_prev$positive,all_prev$tested)
all_prev_nok <- rbind(WA_grouped_mod,malawi_grouped)%>%
  mutate(grav_cat=factor(grav_cat,levels=c('pg','sg','mg')),
         grav_dich=factor(ifelse(grav_cat=='pg','pg','mg'),levels=c('pg','mg')),
         prev = positive/tested)

istp_rainfall <- read_csv('./trial/Data/All_ISTp_rainfall.csv')%>%
  filter(between(as.yearmon(Month),min(all_prev$month),max(all_prev$month)))%>%
  rename(site=Country)%>%
  group_by(site)%>%
  mutate(rel_rainfall = Rainfall/max(Rainfall),
         site = factor(site,levels = c('Ghana', 'Burkina Faso', 'Mali', 'Gambia','Malawi','Kenya')))

istp_rainfall_nok <- read_csv('./trial/Data/All_ISTp_rainfall.csv')%>%
  filter(between(as.yearmon(Month),min(all_prev_nok$month),max(all_prev_nok$month)))%>%
  filter(Country!='Kenya')%>%
  rename(site=Country)%>%
  group_by(site)%>%
  mutate(rel_rainfall = Rainfall/max(Rainfall))

istp_rainfall_k <- read_csv('./trial/Data/All_ISTp_rainfall.csv')%>%
  filter(between(as.yearmon(Month),min(wkenya_group$month),max(wkenya_group$month)))%>%
  filter(Country=='Kenya')%>%
  rename(site=Country)%>%
  group_by(site)%>%
  mutate(rel_rainfall = Rainfall/max(Rainfall))

istp_rainfall_wa <- read_csv('./trial/Data/All_ISTp_rainfall.csv')%>%
  filter(between(as.yearmon(Month),min(WA_all_tests$month),max(WA_all_tests$month)))%>%
  filter(!(Country %in% c('Kenya','Malawi')))%>%
  rename(site=Country)%>%
  group_by(site)%>%
  mutate(rel_rainfall = Rainfall/max(Rainfall),
         site = factor(site,levels = c('Ghana', 'Burkina Faso', 'Mali', 'Gambia')))

istp_rainfall_ea <- read_csv('./trial/Data/All_ISTp_rainfall.csv')%>%
  filter(between(as.yearmon(Month),min(EA_grouped$month),max(EA_grouped$month)))%>%
  filter(Country %in% c('Kenya','Malawi'))%>%
  rename(site=Country)%>%
  group_by(site)%>%
  mutate(rel_rainfall = Rainfall/max(Rainfall),
         site = factor(site,levels = c('Kenya','Malawi')),
         Country = site)
write.csv(istp_rainfall_ea,'./trial/Data/EA_ISTp_rainfall.csv' )

table(all_prev$site,all_prev$grav_cat,useNA = 'always')
windows(10,8)
str(all_prev_cis$site)
levels(all_prev_cis$site)
ggplot(all_prev_cis)+
  geom_line(data=istp_rainfall,aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),col='#999999',size=1)+
  geom_point(aes(x=as.Date(month),y=prev,col=site))+
  geom_errorbar(aes(x=as.Date(month),ymin = lower,ymax = upper,color=site),width=0)+
  facet_grid(site~grav_dich)+
  labs(x='Date',y='ANC prevalence',linetype='Gravidity',color='Location')+
  theme(legend.position = "none")

ggplot(all_prev_nok)+
  geom_line(data=istp_rainfall_nok,aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),col='black',size=1)+
  geom_line(aes(x=as.Date(month),y=prev,col=site,linetype=grav_cat),size=1)+
  facet_grid(site~grav_dich)+
  labs(x='Date',y='ANC prevalence',linetype='Gravidity',color='Location')

ggplot(wkenya_group)+
  geom_line(data=istp_rainfall_k,aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),col='black',size=1)+
  geom_line(aes(x=as.Date(month),y=prev,linetype=grav_cat),size=1,col='#CE5126')+
  facet_grid(grav_dich~.)+
  labs(x='Date',y='ANC prevalence',linetype='Gravidity')

####WA compare test types
WA_slide_pcr <- WA_all_tests%>%
  filter(test_type!='RDT'&!is.na(prev))
ggplot(WA_slide_pcr)+
  geom_line(data=istp_rainfall_wa,aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),col='#999999',size=1)+
  geom_point(aes(x=as.Date(month),y=prev,col=site,shape=test_type,group=test_type,alpha=test_type),position = position_dodge(width = 20),size=3)+
  geom_errorbar(aes(x=as.Date(month),ymin = lower,ymax = upper,color=site,group=test_type,alpha=test_type),width=0,position = position_dodge(width = 20),size=0.6)+
  facet_grid(site~grav_cat)+
  scale_color_discrete(guide = "none")+
  scale_shape_manual(values=c(1,17))+
  scale_alpha_manual(values=c(0.8,1),guide = "none")+
  labs(x='Date',y='ANC prevalence',color='Location',shape='Test')
