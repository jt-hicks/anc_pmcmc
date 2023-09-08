library(R2OpenBUGS)
library(coda)
library(car)
library(MASS)
library(plotrix)
library(binom)
library(scales)
library(zoo)
library(dplyr)
library(plyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(readxl)

## fn to return prevalence from log_odds
get_prev_from_log_odds<-function(log_odds){
  return(exp(log_odds)/(1+exp(log_odds)))
}

## fn to return odds from prevalence
get_odds_from_prev<-function(prev){
  return(prev/(1-prev))
}

addCIs_anc<-function(df,Ys.cs,Ns.cs,Ys.anc,Ns.anc){
  df$mean.cs<-NA
  df$upper.cs<-NA
  df$lower.cs<-NA
  CIs.cs<-binom.confint(Ys.cs,Ns.cs,method="exact")
  df$mean.cs[Ns.cs>0]<-CIs.cs$mean[Ns.cs>0]
  df$upper.cs[Ns.cs>0]<-CIs.cs$upper[Ns.cs>0]
  df$lower.cs[Ns.cs>0]<-CIs.cs$lower[Ns.cs>0]
  df$mean.anc<-NA
  df$upper.anc<-NA
  df$lower.anc<-NA
  CIs.anc<-binom.confint(Ys.anc,Ns.anc,method="exact")
  df$mean.anc[Ns.anc>0]<-CIs.anc$mean[Ns.anc>0]
  df$upper.anc[Ns.anc>0]<-CIs.anc$upper[Ns.anc>0]
  df$lower.anc[Ns.anc>0]<-CIs.anc$lower[Ns.anc>0]
  return(df)
}
##Load in data data grouped by site and month
NG_cs <- readRDS('nnp/data/NG_CS_all_grouped_site.rds')
NG_anc <- readRDS('nnp/data/NG_ANC_mother_grouped_sitegrav.rds')
MZ_cs <- readRDS('nnp/data/MZ_CS_all_grouped_site_180123.rds')
MZ_anc <- readRDS('nnp/data/MZ_ANC_mother_grouped_sitegrav_0822.rds')
BF_cs <- readRDS('nnp/data/BF_CS_all_grouped_site_0822.rds')
BF_anc <- readRDS('nnp/data/BF_ANC_mother_grouped_sitegrav.rds')

BF_anc_sg <- readRDS('nnp/data/BF_ANC_mother_grouped_sitegrav_sg.rds')
MZ_anc_sg <- readRDS('nnp/data/MZ_ANC_mother_grouped_sitegrav_0822_sg.rds')
NG_anc_sg <- readRDS('nnp/data/NG_ANC_mother_grouped_sitegrav_sg.rds')

NG_cs$country <- 'Nigeria'
NG_anc$country <- 'Nigeria'
MZ_cs$country <- 'Mozambique'
MZ_anc$country <- 'Mozambique'
BF_cs$country <- 'Burkina Faso'
BF_anc$country <- 'Burkina Faso'

NG_anc_sg$country <- 'Nigeria'
MZ_anc_sg$country <- 'Mozambique'
BF_anc_sg$country <- 'Burkina Faso'

all_cs <- rbind(NG_cs,MZ_cs,BF_cs)
table(all_cs$month,all_cs$country)
all_anc <- rbind(NG_anc,MZ_anc,BF_anc)
table(all_anc$month,all_anc$country)
all_cs$month_adj <- all_cs$month
all_anc_sg <- rbind(NG_anc_sg,MZ_anc_sg,BF_anc_sg)
##Create an adjusted month variable so that the first month of CX data matches the first month of ANC data
all_cs[all_cs$country=='Burkina Faso'&all_cs$month=='Jun 2020',]$month_adj <- as.yearmon('Sep 2020')
all_cs[all_cs$country=='Mozambique'&(all_cs$month=='Oct 2020'),]$month_adj <- as.yearmon('Dec 2020')
all_cs[all_cs$country=='Nigeria'&(all_cs$month=='Oct 2020'),]$month_adj <- as.yearmon('Nov 2020')
all_cs <- all_cs[,-2] %>%
  dplyr::rename(month=month_adj)

all_both <- merge(all_cs,all_anc, by = c('site','month','country'), suffixes = c('.cs','.anc'))
all_both_total <- all_both %>%
  group_by(country,site,month,.drop=FALSE)%>%
  dplyr::summarise(positive.anc=sum(positive.anc),total.anc=sum(total.anc),positive.cs=mean(positive.cs),total.cs=mean(total.cs),
                   mean.cs=mean(mean.cs),upper.cs=mean(upper.cs),lower.cs=mean(lower.cs),
                   mean.anc=mean(mean.anc),upper.anc=mean(upper.anc),lower.anc=mean(lower.anc))%>%
  mutate(grav_cat = 'All pregnancies')
all_both <- rbind(all_both,all_both_total)
excel.link::xl.save.file(all_both,'nnp/data/nnp_prev4corr_main.xlsx')

all_both_sg <- merge(all_cs,all_anc_sg, by = c('site','month','country'), suffixes = c('.cs','.anc'))
all_both_total_sg <- all_both_sg %>%
  group_by(country,site,month,.drop=FALSE)%>%
  dplyr::summarise(positive.anc=sum(positive.anc),total.anc=sum(total.anc),positive.cs=mean(positive.cs),total.cs=mean(total.cs),
                   mean.cs=mean(mean.cs),upper.cs=mean(upper.cs),lower.cs=mean(lower.cs),
                   mean.anc=mean(mean.anc),upper.anc=mean(upper.anc),lower.anc=mean(lower.anc))%>%
  mutate(grav_cat = 'All pregnancies')
all_both_sg <- rbind(all_both_sg,all_both_total_sg)
excel.link::xl.save.file(all_both_sg,'nnp/data/nnp_prev4corr_sg.xlsx')


vaneijk <- read_excel('nnp/data/paper_data_strat_G.xlsx')
primi_ve <- vaneijk[vaneijk$subG == 'P',3:6]%>%
  dplyr::rename(positive.cs=child_Y,
         total.cs=child_N,
         positive.anc=preg_Y,
         total.anc=preg_N)
multi_ve <- vaneijk[vaneijk$subG == 'M',3:6]%>%
  dplyr::rename(positive.cs=child_Y,
         total.cs=child_N,
         positive.anc=preg_Y,
         total.anc=preg_N)
all_ve <- vaneijk %>%
  subset(subG %in% c('P','M')) %>%
  group_by(Site)%>%
  dplyr::summarise(inf_prev_n=mean(child_Y),
            inf_N=mean(child_N),
            ANC_prev_n=sum(preg_Y),
            ANC_N=sum(preg_N))%>%
  dplyr::select(-Site) %>%
  dplyr::rename(positive.cs=inf_prev_n,
         total.cs=inf_N,
         positive.anc=ANC_prev_n,
         total.anc=ANC_N)
primi_ve <- addCIs_anc(primi_ve,primi_ve$positive.cs,primi_ve$total.cs,primi_ve$positive.anc,primi_ve$total.anc)%>%
  mutate(country='Van Eijk')
multi_ve <- addCIs_anc(multi_ve,multi_ve$positive.cs,multi_ve$total.cs,multi_ve$positive.anc,multi_ve$total.anc)%>%
  mutate(country='Van Eijk')
all_ve <- addCIs_anc(all_ve,all_ve$positive.cs,all_ve$total.cs,all_ve$positive.anc,all_ve$total.anc)%>%
  mutate(country='Van Eijk')
excel.link::xl.save.file(all_ve,'nnp/data/vaneijk_prev4corr_main.xlsx')

all_nnp_total <- all_both_total %>%
  dplyr::mutate(country=dplyr::recode(country,`Nigeria`='NNP - Nigeria',
         `Burkina Faso`='NNP - Burkina Faso',
         `Mozambique`='NNP - Mozambique'))
all_nnp_pg <- all_both[all_both$grav_cat=='Gravidities 1',]%>%
  dplyr::mutate(country=dplyr::recode(country,`Nigeria`='NNP - Nigeria',
                        `Burkina Faso`='NNP - Burkina Faso',
                        `Mozambique`='NNP - Mozambique'))
all_nnp_mg <- all_both[all_both$grav_cat %in% c('Gravidities 2-3','Gravidities 4+'),]%>%
  group_by(country,site,month)%>%
  dplyr::summarise(positive.cs=mean(positive.cs),
            total.cs=mean(total.cs),
            positive.anc=sum(positive.anc),
            total.anc=sum(total.anc))%>%
  dplyr::mutate(country=dplyr::recode(country,`Nigeria`='NNP - Nigeria',
                        `Burkina Faso`='NNP - Burkina Faso',
                        `Mozambique`='NNP - Mozambique'))
all_nnp_mg <- addCIs(all_nnp_mg,all_nnp_mg$positive.cs,all_nnp_mg$total.cs)%>%
  dplyr::rename(mean.cs=mean,
                upper.cs=upper,
                lower.cs=lower)
all_nnp_mg <- addCIs(all_nnp_mg,all_nnp_mg$positive.anc,all_nnp_mg$total.anc)%>%
  dplyr::rename(mean.anc=mean,
                upper.anc=upper,
                lower.anc=lower)

all_nnp_sg <- all_both_sg[all_both_sg$grav_cat=='Gravidities 2',]%>%
  dplyr::mutate(country=dplyr::recode(country,`Nigeria`='NNP - Nigeria',
                                      `Burkina Faso`='NNP - Burkina Faso',
                                      `Mozambique`='NNP - Mozambique'))

names(all_nnp)
names(all_mipmon)
mipmon <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Mozambique/Analysis/mipmon/mipmon4corr.RDS')
all_mipmon_total <- mipmon[mipmon$grav=='All',]%>%
  group_by(site,month)%>%
  dplyr::summarise(positive.cs=mean(positive.cs),
            total.cs=mean(total.cs),
            positive.anc=sum(positive.anc),
            total.anc=sum(total.anc))
all_mipmon_pg <- mipmon[mipmon$grav=='primi',]
all_mipmon_mg <- mipmon[mipmon$grav=='multi',]
mipmon_sg <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Mozambique/Analysis/mipmon/mipmon4corr_sg.RDS')
all_mipmon_sg <- mipmon_sg[mipmon_sg$grav_sg=='secundi',]

all_mipmon_total <- addCIs_anc(all_mipmon_total,all_mipmon_total$positive.cs,all_mipmon_total$total.cs,all_mipmon_total$positive.anc,all_mipmon_total$total.anc)%>%
  mutate(country='MiPMon')
all_mipmon_pg <- addCIs_anc(all_mipmon_pg,all_mipmon_pg$positive.cs,all_mipmon_pg$total.cs,all_mipmon_pg$positive.anc,all_mipmon_pg$total.anc)%>%
  mutate(country='MiPMon')
all_mipmon_mg <- addCIs_anc(all_mipmon_mg,all_mipmon_mg$positive.cs,all_mipmon_mg$total.cs,all_mipmon_mg$positive.anc,all_mipmon_mg$total.anc)%>%
  mutate(country='MiPMon')
all_mipmon_sg <- addCIs_anc(all_mipmon_sg,all_mipmon_sg$positive.cs,all_mipmon_sg$total.cs,all_mipmon_sg$positive.anc,all_mipmon_sg$total.anc)%>%
  mutate(country='MiPMon')
excel.link::xl.save.file(all_mipmon_total,'nnp/data/mipmon_prev4corr_main.xlsx')

all_data_total <- plyr::rbind.fill(all_nnp_total,all_ve,all_mipmon_total)
all_data_pg <- plyr::rbind.fill(all_nnp_pg,primi_ve,all_mipmon_pg)
all_data_mg <- plyr::rbind.fill(all_nnp_mg,multi_ve,all_mipmon_mg)
all_data_sg <- plyr::rbind.fill(all_nnp_sg,all_mipmon_sg)

all_data_total4model <- all_data_total %>%
  filter(country != 'NNP - Burkina Faso') %>%
  filter(!(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month == as.yearmon('Nov 2021')))
all_data_pg4model <- all_data_pg %>%
  filter(country != 'NNP - Burkina Faso') %>%
  filter(!(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month == as.yearmon('Nov 2021')))
all_data_mg4model <- all_data_mg %>%
  filter(country != 'NNP - Burkina Faso') %>%
  filter(!(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month == as.yearmon('Nov 2021')))
all_data_pgmg4model <- merge(all_data_pg4model,all_data_mg4model, by = c('site','month','country','positive.cs','total.cs','mean.cs','upper.cs','lower.cs'), suffixes = c('.pg','.mg'))

all_data_sg4model <- all_data_sg %>%
  filter(country != 'NNP - Burkina Faso') %>%
  filter(!(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month == as.yearmon('Nov 2021')))


country_levels <- c("NNP - Mozambique","NNP - Nigeria","MiPMon","Van Eijk",'NNP - SMC implemented')
country_levels_nosmc <- c("NNP - Mozambique","NNP - Nigeria","MiPMon","Van Eijk")
all_data_total4fig <- all_data_total %>%
  mutate(country = factor(ifelse((country=='NNP - Burkina Faso')|(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month == as.yearmon('Nov 2021')),'NNP - SMC implemented',country),
         levels = country_levels))

all_data_pg4fig <- all_data_pg %>%
  mutate(country = factor(ifelse((country=='NNP - Burkina Faso')|(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month == as.yearmon('Nov 2021')),'NNP - SMC implemented',country),
                          levels = country_levels))
all_data_mg4fig <- all_data_mg %>%
  mutate(country = factor(ifelse((country=='NNP - Burkina Faso')|(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month == as.yearmon('Nov 2021')),'NNP - SMC implemented',country),
                          levels = country_levels))

N_sites=nrow(all_data_total4model)

####All pregnancies####
windows(height=5,width=10)
windows(height=5,width=3.3)
par(mar=c(5,5,1,1))

### now define the model in bugs terms to fit to the data
preg_child_model<-function(){
  for (i in 1:N) {
    pos_child[i] ~ dbin(p_child[i],total_child[i])
    pos_preg[i] ~ dbin(p_preg[i],total_preg[i])
    logit(p_child[i]) <- log_odds_child[i]
    logit(p_preg[i]) <- log_odds_child[i]+log_OR_p_v_c[i]
    log_odds_child[i]~dnorm(av_lo_child,sigma_c_inv)
    log_OR_p_v_c[i]<-RE_intercept[i]+gradient*(log_odds_child[i]-av_lo_child)
    RE_intercept[i]~dnorm(intercept,sigma_int_inv)
  }
  intercept~dnorm(0,0.001)
  gradient ~ dnorm(0,0.001)
  av_lo_child~dnorm(0,0.001)
  sigma_c_inv~dgamma(0.001,0.001)
  sigma_int_inv~dgamma(0.001,0.001)
  sigma_child<-1/sqrt(sigma_c_inv)
  sigma_intercept<-1/sqrt(sigma_int_inv)
}


## write to directory
model.file <- file.path(tempdir(),"model.txt")
write.model(preg_child_model, model.file)

## put in bugs data format
data_list<-list(pos_child=all_data_total4model$positive.cs,
                total_child=all_data_total4model$total.cs,
                pos_preg=all_data_total4model$positive.anc,
                total_preg=all_data_total4model$total.anc,
                N=N_sites)

## fn for random initial conditions
inits<-function(){
  list(av_lo_child=rnorm(1),intercept=rnorm(1),gradient=rnorm(1),sigma_c_inv=runif(1),sigma_int_inv=runif(1))
}
## params of interest
params<-c("av_lo_child","intercept","gradient","sigma_child","sigma_intercept")

## run the model
run_model <- bugs(data_list, inits, params,
                model.file, n.iter=50000,n.burnin=10000)

##attach the output
attach.bugs(run_model)
### should recapture params - something slightly funny with gradient- need to check SD formats and priors
quantile(gradient,c(0.025,0.5,0.975))
quantile(intercept,c(0.025,0.5,0.975))
quantile(av_lo_child,c(0.025,0.5,0.975))
get_prev_from_log_odds(quantile(av_lo_child,c(0.025,0.5,0.975)))
## lets plot the relationship (apols for horrible base R)
prev_child=seq(0.01,0.99,by=0.01)
logodds_child=log(get_odds_from_prev(prev_child))
prev_preg_lower=array(dim=length(logodds_child))
prev_preg_upper=array(dim=length(logodds_child))
for(i in 1:length(logodds_child)){
  prev_preg_lower[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.025))
  prev_preg_upper[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.975))
}
prev_preg_median=get_prev_from_log_odds(logodds_child+median(gradient)*(logodds_child-median(av_lo_child))+median(intercept))
pred_all <- data.frame(prev_child,prev_preg_lower,prev_preg_upper,prev_preg_median)

##Try to make similar graph in ggplot
colors <- c(viridis(4,begin=0.2,end=0.9),"#D95F02")
"#414487FF" "#27808EFF" "#35B779FF" "#BBDF27FF"
colors_nosmc <- c(viridis(4,begin=0.2,end=0.9))

grav_all <- ggplot(all_data_total4fig)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_all,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_all,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='All pregnancies',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
windows(7,7)
grav_all
##plot without smc
grav_all_nosmc <- ggplot(all_data_total4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_all,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_all,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='All pregnancies',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
grav_all_nosmc

########################################  
##Primigrav and Multigrav simultaneously
########################################
### now define the model in bugs terms to fit to the data
preg_mg_child_model<-function(){
  ##For each site (N)
  for (i in 1:N) {
    #Probability distribution of prevalence measures
    pos_child[i] ~ dbin(p_child[i],total_child[i])
    pos_preg_pg[i] ~ dbin(p_preg_pg[i],total_preg_pg[i])
    pos_preg_mg[i] ~ dbin(p_preg_mg[i],total_preg_mg[i])
    
    #Log odds conversions from probability
    logit(p_child[i]) <- log_odds_child[i]
    logit(p_preg_pg[i]) <- log_odds_child[i]+log_OR_pp_v_c[i]
    logit(p_preg_mg[i]) <- log_odds_child[i]+log_OR_pm_v_c[i]
    
    
    log_odds_child[i]~dnorm(av_lo_child,sigma_c_inv)
    
    #Primigrav-specific gradient
    log_OR_pp_v_c[i]<-RE_intercept_pg[i]+gradient_pg*(log_odds_child[i]-av_lo_child)
    #Multigrav-specific gradient
    log_OR_pm_v_c[i]<-RE_intercept_mg[i]+gradient_mg*(log_odds_child[i]-av_lo_child)
    
    #Random-effects intercept for each gravidity type
    RE_intercept_pg[i]~dnorm(intercept_pg,sigma_int_inv)
    RE_intercept_mg[i]~dnorm(intercept_mg,sigma_int_inv)
  }
  intercept_pg~dnorm(0,0.001)
  intercept_mg~dnorm(0,0.001)
  gradient_pg ~ dnorm(0,0.001)
  gradient_mg ~ dnorm(0,0.001)
  av_lo_child~dnorm(0,0.001)
  sigma_c_inv~dgamma(0.001,0.001)
  sigma_int_inv~dgamma(0.001,0.001)
  sigma_child<-1/sqrt(sigma_c_inv)
  sigma_intercept<-1/sqrt(sigma_int_inv)
}

## write to directory
model.file <- file.path(tempdir(),"model.txt")
write.model(preg_mg_child_model, model.file)

## put in bugs data format
data_list<-list(pos_child=all_data_pgmg4model$positive.cs,
                total_child=all_data_pgmg4model$total.cs,
                pos_preg_pg=all_data_pgmg4model$positive.anc.pg,
                total_preg_pg=all_data_pgmg4model$total.anc.pg,
                pos_preg_mg=all_data_pgmg4model$positive.anc.mg,
                total_preg_mg=all_data_pgmg4model$total.anc.mg,
                N=N_sites)

## fn for random initial conditions
inits<-function(){
  list(av_lo_child=rnorm(1),intercept_pg=rnorm(1),intercept_mg=rnorm(1),gradient_pg=rnorm(1),gradient_mg=rnorm(1),sigma_c_inv=runif(1),sigma_int_inv=runif(1))
}
## params of interest
params<-c("av_lo_child","intercept_pg","intercept_mg","gradient_pg","gradient_mg","sigma_child","sigma_intercept")

## run the model
run_model_pgmg <- bugs(data_list, inits, params,
                  model.file, n.iter=50000,n.burnin=10000)

##attach the output
attach.bugs(run_model_pgmg)
### should recapture params - something slightly funny with gradient- need to check SD formats and priors
quantile(gradient_pg,c(0.025,0.5,0.975))
quantile(gradient_mg,c(0.025,0.5,0.975))
quantile(intercept_pg,c(0.025,0.5,0.975))
quantile(intercept_mg,c(0.025,0.5,0.975))
quantile(av_lo_child,c(0.025,0.5,0.975))
get_prev_from_log_odds(quantile(av_lo_child,c(0.025,0.5,0.975)))
## lets plot the relationship (apols for horrible base R)
prev_child=seq(0.01,0.99,by=0.01)
logodds_child=log(get_odds_from_prev(prev_child))
prev_preg_lower=array(dim=length(logodds_child))
prev_preg_upper=array(dim=length(logodds_child))
for(i in 1:length(logodds_child)){
  prev_preg_lower[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept_pg+gradient_pg*(logodds_child[i]-av_lo_child),0.025))
  prev_preg_upper[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept_pg+gradient_pg*(logodds_child[i]-av_lo_child),0.975))
}
prev_preg_median=get_prev_from_log_odds(logodds_child+median(gradient_pg)*(logodds_child-median(av_lo_child))+median(intercept_pg))
pred_pgmg_pg <- data.frame(prev_child,prev_preg_lower,prev_preg_upper,prev_preg_median)
prev_preg_lower=array(dim=length(logodds_child))
prev_preg_upper=array(dim=length(logodds_child))
for(i in 1:length(logodds_child)){
  prev_preg_lower[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept_mg+gradient_mg*(logodds_child[i]-av_lo_child),0.025))
  prev_preg_upper[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept_mg+gradient_mg*(logodds_child[i]-av_lo_child),0.975))
}
prev_preg_median=get_prev_from_log_odds(logodds_child+median(gradient_mg)*(logodds_child-median(av_lo_child))+median(intercept_mg))
pred_pgmg_mg <- data.frame(prev_child,prev_preg_lower,prev_preg_upper,prev_preg_median)

gravpgmg_pg <- ggplot(all_data_pgmg4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc.pg*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc.pg*100,ymax=upper.anc.pg*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc.pg*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_pgmg_pg,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_pgmg_pg,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='Primigravida',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
windows(7,7)
gravpgmg_pg

gravpgmg_mg <- ggplot(all_data_pgmg4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc.mg*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc.mg*100,ymax=upper.anc.mg*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc.mg*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_pgmg_mg,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_pgmg_mg,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='Multigravida',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
windows(7,7)
gravpgmg_mg

###########  
##Primigrav
###########
##Model (preg_child_model) has been defined above
## write to directory
model.file <- file.path(tempdir(),"model.txt")
write.model(preg_child_model, model.file)
## put in bugs data format
data_list_pg<-list(pos_child=all_data_pg4model$positive.cs,
                total_child=all_data_pg4model$total.cs,
                pos_preg=all_data_pg4model$positive.anc,
                total_preg=all_data_pg4model$total.anc,
                N=N_sites)

## fn for random initial conditions
inits<-function(){
  list(av_lo_child=rnorm(1),intercept=rnorm(1),gradient=rnorm(1),sigma_c_inv=runif(1),sigma_int_inv=runif(1))
}
## params of interest
params<-c("av_lo_child","intercept","gradient","sigma_child","sigma_intercept")

## run the model
run_model_pg <- bugs(data_list_pg, inits, params,
                  model.file, n.iter=50000,n.burnin=10000)

##attach the output
attach.bugs(run_model_pg)
### should recapture params - something slightly funny with gradient- need to check SD formats and priors
quantile(gradient,c(0.025,0.5,0.975))
quantile(intercept,c(0.025,0.5,0.975))
quantile(av_lo_child,c(0.025,0.5,0.975))
get_prev_from_log_odds(quantile(av_lo_child,c(0.025,0.5,0.975)))
## lets plot the relationship (apols for horrible base R)
prev_child=seq(0.01,0.99,by=0.01)
logodds_child=log(get_odds_from_prev(prev_child))
prev_preg_lower=array(dim=length(logodds_child))
prev_preg_upper=array(dim=length(logodds_child))
for(i in 1:length(logodds_child)){
  prev_preg_lower[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.025))
  prev_preg_upper[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.975))
}
prev_preg_median=get_prev_from_log_odds(logodds_child+median(gradient)*(logodds_child-median(av_lo_child))+median(intercept))
pred_pg <- data.frame(prev_child,prev_preg_lower,prev_preg_upper,prev_preg_median)

##save a sample of 1000 simulations for use in pmcmc model fitting:
pg_corr_sample <- sample_n(as.data.frame(run_model_pg$sims.matrix),1000)
saveRDS(pg_corr_sample,'nnp/Corr/pg_corr_sample.RDS')
pg_corr_sample <- readRDS('nnp/Corr/pg_corr_sample.RDS')
##Try to make similar graph in ggplot
# grav_pg <- ggplot(all_data_pg4fig)+
#   geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
#   geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
#   geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
#   scale_color_manual(name="Study",values=colors)+
#   scale_y_continuous(limits = c(0,100))+
#   scale_x_continuous(limits = c(0,100))+
#   geom_abline(size=0.8,linetype='dashed')+
#   geom_ribbon(data=pred_pg,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
#   geom_line(data=pred_pg,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
#   theme_bw()+
#   theme(legend.position = 'bottom')+
#   labs(title='Primigravida',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
# windows(7,7)
# grav_pg

##Plot points without smc
grav_pg <- ggplot(all_data_pg4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_pg,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_pg,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='Primigravida',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
windows(7,7)
grav_pg
###########  
##Multi Grav
###########
##Model has been defined above
## put in bugs data format
data_list_mg<-list(pos_child=all_data_mg4model$positive.cs,
                   total_child=all_data_mg4model$total.cs,
                   pos_preg=all_data_mg4model$positive.anc,
                   total_preg=all_data_mg4model$total.anc,
                   N=N_sites-1)

## fn for random initial conditions
inits<-function(){
  list(av_lo_child=rnorm(1),intercept=rnorm(1),gradient=rnorm(1),sigma_c_inv=runif(1),sigma_int_inv=runif(1))
}
## params of interest
params<-c("av_lo_child","intercept","gradient","sigma_child","sigma_intercept")

## run the model
run_model_mg <- bugs(data_list_mg, inits, params,
                     model.file, n.iter=50000,n.burnin=10000)

##attach the output
attach.bugs(run_model_mg)
### should recapture params - something slightly funny with gradient- need to check SD formats and priors
quantile(gradient,c(0.025,0.5,0.975))
quantile(intercept,c(0.025,0.5,0.975))
quantile(av_lo_child,c(0.025,0.5,0.975))
get_prev_from_log_odds(quantile(av_lo_child,c(0.025,0.5,0.975)))
## lets plot the relationship (apols for horrible base R)
prev_child=seq(0.01,0.99,by=0.01)
logodds_child=log(get_odds_from_prev(prev_child))
prev_preg_lower=array(dim=length(logodds_child))
prev_preg_upper=array(dim=length(logodds_child))
for(i in 1:length(logodds_child)){
  prev_preg_lower[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.025))
  prev_preg_upper[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.975))
}
prev_preg_median=get_prev_from_log_odds(logodds_child+median(gradient)*(logodds_child-median(av_lo_child))+median(intercept))
pred_mg <- data.frame(prev_child,prev_preg_lower,prev_preg_upper,prev_preg_median)

##Try to make similar graph in ggplot
# grav_mg <- ggplot(all_data_mg4fig)+
#   geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
#   geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
#   geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
#   scale_color_manual(name="Study",values=colors)+
#   scale_y_continuous(limits = c(0,100))+
#   scale_x_continuous(limits = c(0,100))+
#   geom_abline(size=0.8,linetype='dashed')+
#   geom_ribbon(data=pred_mg,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
#   geom_line(data=pred_mg,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
#   theme_bw()+
#   theme(legend.position = 'bottom')+
#   labs(title='Multigravida',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
# grav_mg

###########  
##Secundigrav
###########
##Model (preg_child_model) has been defined above
## write to directory
N_sites=nrow(all_data_sg4model)

model.file <- file.path(tempdir(),"model.txt")
write.model(preg_child_model, model.file)
## put in bugs data format
data_list_sg<-list(pos_child=all_data_sg4model$positive.cs,
                   total_child=all_data_sg4model$total.cs,
                   pos_preg=all_data_sg4model$positive.anc,
                   total_preg=all_data_sg4model$total.anc,
                   N=N_sites)

## fn for random initial conditions
inits<-function(){
  list(av_lo_child=rnorm(1),intercept=rnorm(1),gradient=rnorm(1),sigma_c_inv=runif(1),sigma_int_inv=runif(1))
}
## params of interest
params<-c("av_lo_child","intercept","gradient","sigma_child","sigma_intercept")

## run the model
run_model_sg <- bugs(data_list_sg, inits, params,
                     model.file, n.iter=50000,n.burnin=10000)

##attach the output
attach.bugs(run_model_sg)
### should recapture params - something slightly funny with gradient- need to check SD formats and priors
quantile(gradient,c(0.025,0.5,0.975))
quantile(intercept,c(0.025,0.5,0.975))
quantile(av_lo_child,c(0.025,0.5,0.975))
get_prev_from_log_odds(quantile(av_lo_child,c(0.025,0.5,0.975)))
## lets plot the relationship (apols for horrible base R)
prev_child=seq(0.01,0.99,by=0.01)
logodds_child=log(get_odds_from_prev(prev_child))
prev_preg_lower=array(dim=length(logodds_child))
prev_preg_upper=array(dim=length(logodds_child))
for(i in 1:length(logodds_child)){
  prev_preg_lower[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.025))
  prev_preg_upper[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.975))
}
prev_preg_median=get_prev_from_log_odds(logodds_child+median(gradient)*(logodds_child-median(av_lo_child))+median(intercept))
pred_sg <- data.frame(prev_child,prev_preg_lower,prev_preg_upper,prev_preg_median)

##save a sample of 1000 simulations for use in pmcmc model fitting:
sg_corr_sample <- sample_n(as.data.frame(run_model_sg$sims.matrix),1000)
saveRDS(sg_corr_sample,'nnp/Corr/sg_corr_sample.RDS')

##No smc
grav_mg <- ggplot(all_data_mg4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_mg,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_mg,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='Multigravida',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
windows(7,7)
grav_mg
mg_corr_sample <- sample_n(as.data.frame(run_model_mg$sims.matrix),1000)
saveRDS(mg_corr_sample,'nnp/Corr/mg_corr_sample.RDS')

pgmg_corr_sample <- sample_n(as.data.frame(run_model_pgmg$sims.matrix),1000)
saveRDS(pgmg_corr_sample,'nnp/Corr/pgmg_corr_sample.RDS')

windows(height=5.5,width=13)

grav_all_nosmc + grav_pg + grav_mg + 
  plot_annotation(title = 'PG and MG Fitted Independently') + plot_layout(ncol=3,guides = 'collect')
ggsave('nnp/Corr/correlation_ind_150223.pdf',height=5.5,width=13)
grav_all_nosmc + gravpgmg_pg + gravpgmg_mg + 
  plot_annotation(title = 'PG and MG Fitted Together') + plot_layout(ncol=3,guides = 'collect')
ggsave('nnp/Corr/correlation_comb_150223.pdf',height=5.5,width=13)

##Format predicted relationship data
pred_all$model <- 'All ANC'
pred_pg$model <- 'Independent'
pred_mg$model <- 'Independent'
pred_pgmg_pg$model <- 'Together'
pred_pgmg_mg$model <- 'Together'

pred_all$grav <- 'All ANC'
pred_pg$grav <- 'Primigrav'
pred_mg$grav <- 'Multigrav'
pred_pgmg_pg$grav <- 'Primigrav'
pred_pgmg_mg$grav <- 'Multigrav'

pred_combined <- rbind(pred_all,pred_pg,pred_mg,pred_pgmg_pg,pred_pgmg_mg)
colors_model <- c(viridis(2,begin=0.25,end=0.75))
c(viridis(2,begin=0.5,end=0.75))
"#414487FF" "#BBDF27FF"
"#482576FF" "#7AD151FF"
"#3B528BFF" "#5DC863FF"
"#21908CFF" "#5DC863FF"
grav_all_4comp <- ggplot(all_data_total4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100),col='grey',size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100),col='grey',width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100),col='grey',height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_all,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_all,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='All pregnancies',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
windows(7,7)
grav_all_4comp


model_comp_pg <- ggplot(all_data_pg4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100),col='grey',size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100),col='grey',width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100),col='grey',height=0)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_combined[pred_combined$grav=='Primigrav',],aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100,fill=model),alpha=0.2)+
  geom_line(data=pred_combined[pred_combined$grav=='Primigrav',],aes(x=prev_child*100,y=prev_preg_median*100,col=model),size=1)+
  scale_color_manual(name="Model",values=colors_model)+
  scale_fill_manual(name="Model",values=colors_model)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='Primigravida',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
windows(7,7)
model_comp_pg

model_comp_mg <- ggplot(all_data_mg4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100),col='grey',size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100),col='grey',width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100),col='grey',height=0)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_combined[pred_combined$grav=='Multigrav',],aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100,fill=model),alpha=0.2)+
  geom_line(data=pred_combined[pred_combined$grav=='Multigrav',],aes(x=prev_child*100,y=prev_preg_median*100,col=model),size=1)+
  scale_color_manual(name="Model",values=colors_model)+
  scale_fill_manual(name="Model",values=colors_model)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='Multigravida',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
windows(7,7)
model_comp_mg

model_comp_sg <- ggplot(all_data_sg4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100),col='grey',size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100),col='grey',width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100),col='grey',height=0)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_sg,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_sg,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  # scale_color_manual(name="Model",values=colors_model)+
  # scale_fill_manual(name="Model",values=colors_model)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='Secundigravida',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
windows(7,7)
model_comp_sg

windows(height=5.5,width=13)
grav_all_4comp + model_comp_pg + model_comp_mg + 
  plot_annotation(title = 'Model Comparison') + plot_layout(ncol=3,guides = 'collect')
ggsave('nnp/Corr/correlation_compari_150223.pdf',height=5.5,width=13)

run_model$DIC
saveRDS(run_model,'nnp/Corr/corr_run_all_150223.RDS')
run_model_pg$DIC
saveRDS(run_model_pg,'nnp/Corr/corr_run_pg_150223.RDS')
run_model_mg$DIC
saveRDS(run_model_mg,'nnp/Corr/corr_run_mg_150223.RDS')
run_model_pgmg$DIC
saveRDS(run_model_pgmg,'nnp/Corr/corr_run_pgmg_150223.RDS')
run_model$last.values
plot(run_model)

##Extract posterior distribution of coefficients and combine multiple models for comparison
run_model_df <- as.data.frame(run_model$sims.matrix)%>%
  mutate(grav = 'All ANC',
         model = 'All ANC',
         code = '1')
run_model_pg_df <- as.data.frame(run_model_pg$sims.matrix)%>%
  mutate(grav = 'Primigrav',
         model = 'Independent',
         code = '2')
run_model_mg_df <- as.data.frame(run_model_mg$sims.matrix)%>%
  mutate(grav = 'Multigrav',
         model = 'Independent',
         code = '3')
run_model_pgmg_df <- as.data.frame(run_model_pgmg$sims.matrix)%>%
  mutate(model = 'Together',
         code = '4')
run_model_pgmg_pg_df <- run_model_pgmg_df%>%
  dplyr::select(!c('intercept_mg','gradient_mg'))%>%
  mutate(grav = 'Primigrav')%>%
  dplyr::rename(intercept = intercept_pg,
         gradient = gradient_pg)
run_model_pgmg_mg_df <- run_model_pgmg_df%>%
  dplyr::select(!c(intercept_pg,gradient_pg))%>%
  mutate(grav = 'Multigrav')%>%
  dplyr::rename(intercept = intercept_mg,
         gradient = gradient_mg)

run_model_combined <- rbind(run_model_df,run_model_pg_df,run_model_mg_df,
                            run_model_pgmg_pg_df,run_model_pgmg_mg_df)%>%
  mutate(av_prev_child = get_prev_from_log_odds(av_lo_child))

library(ggridges)
windows(10,7)
intercept_comp <- ggplot(run_model_combined)+
  geom_density_ridges2(aes(x=intercept,y=model,fill=grav),scale=1)+
  facet_grid(grav~.,scales = 'free_y')+
  scale_fill_viridis(discrete = TRUE,begin=0.2,end=0.9)+
  labs(title='Intercept Posterior Distribution',x='Intercept Value',y='Model')
gradient_comp <- ggplot(run_model_combined)+
  geom_density_ridges2(aes(x=gradient,y=model,fill=grav),scale=1)+
  facet_grid(grav~.,scales = 'free_y')+
  scale_fill_viridis(discrete = TRUE,begin=0.2,end=0.9)+
  labs(title='Gradient Posterior Distribution',x='Gradient Value',y='Model')
windows(8,5.5)
intercept_comp + gradient_comp + 
  plot_annotation(title = 'Model Comparison') + plot_layout(ncol=2,guides = 'collect')
ggsave('nnp/Corr/parameter_compari_150223.pdf',height=5.5,width=8)

child_prev_plot <- ggplot(run_model_combined[!(run_model_combined$model=='Together'&run_model_combined$grav=='Multigrav'),])+
  geom_density_ridges2(aes(x=av_prev_child,y=code),scale=1)+
  scale_y_discrete(labels=c('All ANC', 'Primigrav\nIndepedent', 'Multigrav\nIndepedent', 'Together'))+
  labs(x='Child Prevalence',y='Model')
sigma_child_plot <- ggplot(run_model_combined[!(run_model_combined$model=='Together'&run_model_combined$grav=='Multigrav'),])+
  geom_density_ridges2(aes(x=sigma_child,y=code),scale=1)+
  scale_y_discrete(labels=c('All ANC', 'Primigrav\nIndepedent', 'Multigrav\nIndepedent', 'Together'))+
  labs(x='Sigma - childhood prevalance',y='Model')
sigma_intercept_plot <- ggplot(run_model_combined[!(run_model_combined$model=='Together'&run_model_combined$grav=='Multigrav'),])+
  geom_density_ridges2(aes(x=sigma_intercept,y=code),scale=1)+
  scale_y_discrete(labels=c('All ANC', 'Primigrav\nIndepedent', 'Multigrav\nIndepedent', 'Together'))+
  labs(x='Sigma - Intercept',y='Model')
windows(13,5.5)
child_prev_plot + sigma_child_plot + sigma_intercept_plot + 
  plot_annotation(title = 'Model Comparison') + plot_layout(ncol=3)
ggsave('nnp/Corr/parameter_compari_150223-2.pdf',height=5.5,width=13)
