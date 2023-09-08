##summarize results##
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(zoo)
library(patchwork)
library(RColorBrewer)
library(bayesplot)
library(binom)
source('shared/addCIs.R')
library(lubridate)
library(zoo)
library(geofacet)
library(sf)
library(terra)
library(spData)
library(spDataLarge)
library(tmap)
library(excel.link)
library(readxl)

theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))

lindi_4 <- readRDS('./tanz/smc_sim/lindi_4.RDS')
lindi_5 <- readRDS('./tanz/smc_sim/lindi_5.RDS')
mtwara_4 <- readRDS('./tanz/smc_sim/mtwara_4.RDS')
mtwara_5 <- readRDS('./tanz/smc_sim/mtwara_5.RDS')
ruvuma_4 <- readRDS('./tanz/smc_sim/ruvuma_4.RDS')
ruvuma_5 <- readRDS('./tanz/smc_sim/ruvuma_5.RDS')

lindi_4$rounds <- 'Four'
lindi_4$site <- 'Lindi'
lindi_5$rounds <- 'Five'
lindi_5$site <- 'Lindi'
mtwara_4$rounds <- 'Four'
mtwara_4$site <- 'Mtwara'
mtwara_5$rounds <- 'Five'
mtwara_5$site <- 'Mtwara'
ruvuma_4$rounds <- 'Four'
ruvuma_4$site <- 'Ruvuma'
ruvuma_5$rounds <- 'Five'
ruvuma_5$site <- 'Ruvuma'

lindi_none <- lindi_5%>%
  select(c(1:4,8:9))
lindi_none$rounds <- 'None'
mtwara_none <- mtwara_5%>%
  select(c(1:4,8:9))
mtwara_none$rounds <- 'None'
ruvuma_none <- ruvuma_5%>%
  select(c(1:4,8:9))
ruvuma_none$rounds <- 'None'


all_three_smc_sites <- rbind(lindi_4,lindi_5,mtwara_4,mtwara_5,ruvuma_4,ruvuma_5)
all_three_smc_sites <- all_three_smc_sites %>%
  select(c(1,5:9))%>%
  rename(incu5_med = incu5_smc_med,
         incu5_lower = incu5_smc_lower,
         incu5_upper = incu5_smc_upper)
nrow(smc_formatted)
range(smc_formatted$t)
smc_formatted <- rbind(all_three_smc_sites,lindi_none,mtwara_none,ruvuma_none)%>%
  mutate(rounds=factor(rounds,levels=c('None','Four','Five')))
smc_formatted$date <- rep(seq.Date(from=as.Date('2015-01-01'),length.out = 2892,by='day'),9)

ANC_rainfall_region<-read_csv("./tanz/processed_inputs/ANC_data_Rainfall_ad1.csv")%>%
  mutate(yearmon=as.yearmon(yearmon))

rainfall_lmr <- ANC_rainfall_region%>%
  filter(Region %in% c('Lindi','Mtwara','Ruvuma'))%>%
  mutate(scale=ifelse(Region=='Lindi',9,ifelse(Region=='Mtwara',7.5,4.5)))%>%
  group_by(Region)%>%
  mutate(std_rainfall=(Rainfall - mean(Rainfall))/sd(Rainfall),
         rel_rainfall=Rainfall*scale/max(Rainfall))%>%
  rename(site=Region)

windows(10,6)
scales::show_col(viridis(3))
viridis(3)
colors <- c(None='white',Four='white',Five='white')
colors <- c(None='white',Four='white',Five='white')
colors <- c(None='#FDE725FF',Four='#21908CFF',Five='#440154FF')

rain <- ggplot(smc_formatted)+
  geom_line(aes(x=date,y=incu5_med*1000,color=rounds),size=0.8)+
  geom_line(data=rainfall_lmr,aes(x=as.Date(as.yearmon(yearmon)),y=rel_rainfall),color="#1f78b4",size=0.8)+
  # geom_ribbon(aes(x=date,ymin=incu5_lower*1000,ymax=incu5_upper*1000,fill=rounds),alpha=0.2)+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  facet_grid(site~.,scales = 'free_y')+
  coord_cartesian(xlim=as.Date(c('2015-01-01','2022-31-12')))+
  labs(x='Year',
       y='Relative Rainfall', color='Number SMC Rounds',fill='Number SMC Rounds')
colors <- c(None='#666666',Four='#CE5126',Five='#6a3d9a')

ggsave('./tanz/smc_sim/rainfall.pdf',plot=rain,width=10,height = 6)
no_smc <-  ggplot(smc_formatted[smc_formatted$rounds=='None',])+
  geom_line(data=rainfall_lmr,aes(x=as.Date(as.yearmon(yearmon)),y=rel_rainfall),color="#1f78b4",size=0.8)+
  geom_line(aes(x=date,y=incu5_med*1000,color=rounds),size=0.8)+
  # geom_ribbon(aes(x=date,ymin=incu5_lower*1000,ymax=incu5_upper*1000,fill=rounds),alpha=0.2)+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  facet_grid(site~.,scales = 'free_y')+
  coord_cartesian(xlim=as.Date(c('2015-01-01','2022-31-12')))+
  labs(x='Year',
       y='Clinical cases per day per 1000 children under 5 years old', color='Number SMC Rounds',fill='Number SMC Rounds')
ggsave('./tanz/smc_sim/no_smc.pdf',plot=no_smc,width=10,height = 6)

four_smc <-  ggplot(smc_formatted[smc_formatted$rounds%in%c('None','Four'),])+
  geom_line(data=rainfall_lmr,aes(x=as.Date(as.yearmon(yearmon)),y=rel_rainfall),color="#1f78b4",size=0.8)+
  geom_line(aes(x=date,y=incu5_med*1000,color=rounds),size=0.8)+
  # geom_ribbon(aes(x=date,ymin=incu5_lower*1000,ymax=incu5_upper*1000,fill=rounds),alpha=0.2)+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  facet_grid(site~.,scales = 'free_y')+
  coord_cartesian(xlim=as.Date(c('2015-01-01','2022-31-12')))+
  labs(x='Year',
       y='Clinical cases per day per 1000 children under 5 years old', color='Number SMC Rounds',fill='Number SMC Rounds')
ggsave('./tanz/smc_sim/four_smc.pdf',plot=four_smc,width=10,height = 6)

five_smc <-  ggplot(smc_formatted)+
  geom_line(data=rainfall_lmr,aes(x=as.Date(as.yearmon(yearmon)),y=rel_rainfall),color="#1f78b4",size=0.8)+
  geom_line(aes(x=date,y=incu5_med*1000,color=rounds),size=0.8)+
  # geom_ribbon(aes(x=date,ymin=incu5_lower*1000,ymax=incu5_upper*1000,fill=rounds),alpha=0.2)+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  facet_grid(site~.,scales = 'free_y')+
  coord_cartesian(xlim=as.Date(c('2015-01-01','2022-31-12')))+
  labs(x='Year',
       y='Clinical cases per day per 1000 children under 5 years old', color='Number SMC Rounds',fill='Number SMC Rounds')
ggsave('./tanz/smc_sim/five_smc.pdf',plot=five_smc,width=10,height = 6)
