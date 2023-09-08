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
results <- tanz_lakemal_lt20_2015to2017_results
results <- tanz_lt20_2015to2017_results
create_diag_figs_bulk <- function(results){
  ar.df <- data.frame(sites = names(results),
                      acceptance_rate = sapply(1:length(results), function(x){
                        1 - coda::rejectionRate(as.mcmc(results[[x]]$mcmc))[[1]]})
  )

  ess.df <- bind_rows(lapply(1:length(results), 
                                   function(x){
                                     data.frame(sites = names(results[x]),
                                                t(coda::effectiveSize(as.mcmc(results[[x]]$mcmc))))
                                     }))%>%
    melt(id.vars = 'sites')
  
  mcmc.df <- bind_rows(lapply(1:length(results), 
                      function(x){
                        df <- results[[x]]$mcmc
                        df$sites <- names(results[x])
                        df$step <- 1:nrow(df)
                        return(df)
                      }))

  ar.plot <- ggplot(ar.df,aes(x=sites,y=acceptance_rate))+
    geom_col()+
    labs(title='Acceptance Rate by Site',x='Site',y='Acceptance Rate')+
    geom_text(aes(label = round(acceptance_rate,digits=2)), vjust = -0.5)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))

  ess.plot <- ggplot(ess.df,aes(x=variable,y=value))+
    geom_col()+
    facet_wrap(vars(sites))+
    labs(title='ESS Values by Parameter and Site',x='Parameter',y='Effective Size')+
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  prior.trace <- ggplot(mcmc.df)+
    geom_line(aes(x=step,y=log_prior))+
    facet_wrap(vars(sites),scales='free_y')+
    labs(title='Log Prior Trace by Site',y='Log Prior')+
    theme(axis.title.x = element_blank())
  
  posterior.trace <- ggplot(mcmc.df)+
    geom_line(aes(x=step,y=log_posterior))+
    facet_wrap(vars(sites),scales='free_y')+
    labs(title='Log Posterior Trace by Site',y='Log Posterior')+
    theme(axis.title.x = element_blank())

  likelihood.trace <- ggplot(mcmc.df)+
    geom_line(aes(x=step,y=log_likelihood))+
    facet_wrap(vars(sites),scales='free_y')+
    labs(title='Log Likelihood Trace by Site',y='Log Likelihood')+
    theme(axis.title.x = element_blank())
  
  EIR_SD.trace <- ggplot(mcmc.df)+
    geom_line(aes(x=step,y=EIR_SD))+
    facet_wrap(vars(sites),scales='free_y')+
    labs(title='EIR_SD Trace by Site',y='EIR_SD')+
    theme(axis.title.x = element_blank())
  
  init_EIR.trace <- ggplot(mcmc.df)+
    geom_line(aes(x=step,y=log_init_EIR))+
    facet_wrap(vars(sites),scales='free_y')+
    labs(title='Log init_EIR Trace by Site',y='Log init_EIR')+
    theme(axis.title.x = element_blank())
  
  EIR_SD.density <- ggplot(mcmc.df)+
    geom_violin(aes(x=sites,y=EIR_SD),fill='darkgray')+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(title='EIR_SD density by Site',x='Site',y='EIR_SD')
  
  log_init_EIR.density <- ggplot(mcmc.df)+
    geom_violin(aes(x=sites,y=log_init_EIR),fill='darkgray')+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(title='Log init_EIR density by Site',x='Site',y='Log init_EIR')
  
  init_EIR.density <- ggplot(mcmc.df)+
    geom_violin(aes(x=sites,y=exp(log_init_EIR)),fill='darkgray')+
    scale_y_log10()+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(title='init_EIR density by Site',x='Site',y='init_EIR')
    
  diags <- list(ar.plot = ar.plot, ess.plot = ess.plot, 
                prior.trace = prior.trace, posterior.trace = posterior.trace, 
                likelihood.trace = likelihood.trace, EIR_SD.trace =EIR_SD.trace,
                init_EIR.trace = init_EIR.trace, EIR_SD.density = EIR_SD.density,
                log_init_EIR.density = log_init_EIR.density, init_EIR.density = init_EIR.density)
  
  
  return(diags)
}

diagnostic_plots <- create_diag_figs_bulk(tanz_lt20_2015to2017_results)
diagnostic_plots <- create_diag_figs_bulk(tanz_all_2015to2022_results)
diagnostic_plots$ar.plot
diagnostic_plots$posterior.trace
diagnostic_plots$EIR_SD.density
diagnostic_plots$init_EIR.density
data_list <- tanz_data_list_15to17$lt20
level
rainfall <- tanz_rainfall
results <- tanz_lt20_2015to2017_results
level <- 'Region'
create_summary_plots <- function(results,data_list,rainfall,
                                 level=c('Region','Council'),
                                 date_limits = as.Date(c(NA,NA))){
  mcmc.df <- bind_rows(lapply(1:length(results), 
                              function(x){
                                df <- results[[x]]$mcmc[101:1000,]
                                df$sites <- names(results[x])
                                df$step <- 1:nrow(df)
                                return(df)
                              }))
  prev.df <- bind_rows(lapply(1:length(results), 
                              function(x){
                                history.df <- as.data.frame(t(results[[x]]$history['prev_05', 101:1000, -1]))
                                prev_history <- history.df%>%
                                  dplyr::mutate(t=c(1:nrow(history.df)))%>%
                                  melt(id='t')%>%
                                  dplyr::rename(time=t)%>%
                                  group_by(time)%>%
                                  dplyr::summarise(median=median(value),
                                            mean=mean(value),
                                            upper=quantile(value,probs=0.975),
                                            lower=quantile(value,probs=0.025))%>%
                                 dplyr::mutate(sites = names(results[x]),
                                        month = as.Date(data_list[[x]]$month))
                                return(prev_history)
                              }))

  inc.df <- bind_rows(lapply(1:length(results), 
                              function(x){
                                history.df <- as.data.frame(t(results[[x]]$history['clininc_all', 101:1000, -1]))
                                inc_history <- history.df%>%
                                  mutate(t=c(1:nrow(history.df)))%>%
                                  melt(id='t')%>%
                                  dplyr::rename(time=t)%>%
                                  group_by(time)%>%
                                  dplyr::summarise(median=median(value),
                                                   mean=mean(value),
                                                   upper=quantile(value,probs=0.975),
                                                   lower=quantile(value,probs=0.025))%>%
                                  mutate(sites = names(results[x]),
                                         month = as.Date(data_list[[x]]$month))
                                return(inc_history)
                              }))

  rainfall$sites <- sapply(1:nrow(rainfall), function(x) gsub(' Region','',rainfall[x,level]))
  rainfall$month <- as.Date(as.yearmon(rainfall$yearmon))
  
  inc.rainfall.df <- bind_rows(lapply(1:length(results), 
                             function(x){
                               history.df <- as.data.frame(t(results[[x]]$history['clininc_all', 101:1000, -1]))
                               inc_history <- history.df%>%
                                 mutate(t=c(1:nrow(history.df)))%>%
                                 melt(id='t')%>%
                                 dplyr::rename(time=t)%>%
                                 group_by(time)%>%
                                 dplyr::summarise(inc.median=median(value),
                                                  inc.mean=mean(value),
                                                  inc.upper=quantile(value,probs=0.975),
                                                  inc.lower=quantile(value,probs=0.025))%>%
                                 mutate(sites = names(results[x]),
                                        month = as.Date(data_list[[x]]$month))
                               inc_plus_rainfall <- left_join(inc_history,rainfall,by=c('month','sites'))%>%
                                 mutate(rainfall_norm = Rainfall/max(Rainfall),
                                        rainfall_maxinc = Rainfall * (max(inc.median)/max(Rainfall)))
                               history.df.prev <- as.data.frame(t(results[[x]]$history['prev_05', 101:1000, -1]))
                               prev_history <- history.df.prev%>%
                                 mutate(t=c(1:nrow(history.df.prev)))%>%
                                 melt(id='t')%>%
                                 dplyr::rename(time=t)%>%
                                 group_by(time)%>%
                                 dplyr::summarise(prev.median=median(value),
                                                  prev.mean=mean(value),
                                                  prev.upper=quantile(value,probs=0.975),
                                                  prev.lower=quantile(value,probs=0.025))%>%
                                 mutate(sites = names(results[x]),
                                        month = as.Date(data_list[[x]]$month))
                               all <- left_join(inc_plus_rainfall,prev_history,by=c('month','sites'))%>%
                                 mutate(prev_maxinc = prev.median * (max(inc.median)),
                                        upper_maxinc = prev.upper *max(inc.median),
                                        lower_maxinc = prev.lower *max(inc.median),
                                        rainfall_maxprev = Rainfall * (max(prev.median)/max(Rainfall)),
                                        inc_maxinc = inc.median/ max(inc.median),
                                        upper_maxinc = inc.upper /max(inc.median),
                                        lower_maxinc = inc.lower /max(inc.median),
                                        rainfall_maxprev = Rainfall * (max(prev.median)/max(Rainfall)))
                               return(all)
                             }))
  inc.rainfall.df <- addCIs(inc.rainfall.df,inc.rainfall.df$positive,inc.rainfall.df$tested)
  
  eir.df <- bind_rows(lapply(1:length(results), 
                             function(x){
                               history.df <- as.data.frame(t(results[[x]]$history['EIR', 101:1000, -1]))
                               eir_history <- history.df%>%
                                 mutate(t=c(1:nrow(history.df)))%>%
                                 melt(id='t')%>%
                                 dplyr::rename(time=t)%>%
                                 group_by(time)%>%
                                 dplyr::summarise(median=median(value),
                                                  mean=mean(value),
                                                  upper=quantile(value,probs=0.975),
                                                  lower=quantile(value,probs=0.025))%>%
                                 mutate(sites = names(results[x]),
                                        month = data_list[[x]]$month)
                               return(eir_history)
                             }))
  EIR_SD.density <- ggplot(mcmc.df)+
    geom_violin(aes(x=sites,y=EIR_SD),fill='#6D6A67')+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(title='EIR_SD density by Site',x='Site',y='EIR_SD')
  
  log_init_EIR.density <- ggplot(mcmc.df)+
    geom_violin(aes(x=sites,y=log_init_EIR),fill='#6D6A67')+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(title='Log init_EIR density by Site',x='Site',y='Log init_EIR')
  
  init_EIR.density <- ggplot(mcmc.df)+
    geom_violin(aes(x=sites,y=exp(log_init_EIR)),fill='#6D6A67')+
    scale_y_log10()+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(title='init_EIR density by Site',x='Site',y='init_EIR')
  # "#2A788EFF" "#5DC863FF"
  if(level=='Region'){
    province_grid<-read_csv("./tanz/processed_inputs/Province_grid.csv")
    province_grid$name=province_grid$code=province_grid$NAME_1
    scaler <- 0.0075
    
  }else{
    province_grid<-as.data.frame(read_csv("./tanz/processed_inputs/Lake_Malawi_grid.csv"))
    province_grid <- province_grid[province_grid$name != 'Songea District Council',]
    scaler <- max(inc.rainfall.df$median)
  }

  obs_prev_plot <- ggplot(inc.rainfall.df[!(inc.rainfall.df$sites=='Mbeya'&inc.rainfall.df$yearmon=='Oct 2020'),])+
    geom_line(aes(x=month,y=rainfall_norm*0.3),col="#1582AD",size=0.8)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper),alpha=0.2)+
    # geom_line(aes(x=month,y=median),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(aes(x=month,y=mean),pch = 19,color='#6D6A67')+
    geom_errorbar(aes(x=month,ymin=lower,ymax=upper),width = 0,color='#6D6A67')+
    # scale_color_manual(values=colors)+
    # scale_fill_manual(values=colors)+
    facet_geo(~ sites, grid = province_grid%>%
                select(row,col,code,name))+
    scale_x_date(date_labels = "'%y",limits = date_limits)+
    labs(x='Date',y='ANC Prevalence')
  est_prev_plot <- ggplot(inc.rainfall.df[!(inc.rainfall.df$sites=='Mbeya'&inc.rainfall.df$yearmon=='Oct 2020'),])+
    geom_line(aes(x=month,y=rainfall_norm*0.3),col="#1582AD",size=0.8)+
    geom_point(aes(x=month,y=mean),pch = 19,color='#6D6A67')+
    geom_errorbar(aes(x=month,ymin=lower,ymax=upper),width = 0,color='#6D6A67')+
    geom_ribbon(aes(x=month,ymin=prev.lower,ymax=prev.upper),alpha=0.2,fill="#EFBB12")+
    geom_line(aes(x=month,y=prev.median),size=1,color="#EFBB12")+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    # scale_color_manual(values=colors)+
    # scale_fill_manual(values=colors)+
    facet_geo(~ sites, grid = province_grid%>%
                select(row,col,code,name))+
    scale_x_date(date_labels = "'%y",limits = date_limits)+
    labs(x='Date',y='ANC Prevalence')
  inc_plot <- ggplot(inc.df)+
    geom_ribbon(aes(x=month,ymin=lower,ymax=upper),alpha=0.2)+
    geom_line(aes(x=month,y=median),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    # geom_point(data=inc.rainfall.df,aes(x=month,y=mean),pch = 19)+
    # geom_errorbar(data=inc.rainfall.df,aes(x=month,ymin=lower,ymax=upper),width = 0)+
    # scale_color_manual(values=colors)+
    # scale_fill_manual(values=colors)+
    facet_geo(~ sites, grid = province_grid%>%
                select(row,col,code,name))+
    scale_x_date(date_labels = "'%y",limits = date_limits)+
    labs(x='Date',y='Incidence')
    
  corr_plot <- ggplot(inc.rainfall.df)+
    geom_point(aes(x=Rainfall,y=inc.median))+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    # geom_point(data=inc.rainfall.df,aes(x=month,y=mean),pch = 19)+
    # geom_errorbar(data=inc.rainfall.df,aes(x=month,ymin=lower,ymax=upper),width = 0)+
    # scale_color_manual(values=colors)+
    # scale_fill_manual(values=colors)+
    facet_geo(~ sites, grid = province_grid%>%
                select(row,col,code,name), scale = 'free')+
    labs(x='rainfall',y='Incidence')+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  inc_rf_corr <- inc.rainfall.df%>%
    group_by(sites)%>%
    dplyr::summarise(corr=cor(inc.median,Rainfall))
  
  # corr_map <- tm_shape(tz) +
  #   tm_borders()
  inc.rainfall <- ggplot(inc.rainfall.df)+
    geom_line(aes(x=month,y=rainfall_norm*scaler),col="#1582AD",size=0.8)+
    geom_ribbon(aes(x=month,ymin=inc.lower,ymax=inc.upper),alpha=0.4,fill="#CE5126")+
    geom_line(aes(x=month,y=inc.median),size=0.8,col="#CE5126")+
    scale_x_date(date_labels = "'%y",limits = date_limits)+
    facet_geo(~ sites, grid = province_grid%>%
                               select(row,col,code,name))+
    ylab("Estimated incidence/normalised rainfall")+xlab("Year")+
    coord_cartesian(ylim=c(0,scaler))

  inc.rainfall.3 <- ggplot(inc.rainfall.df)+
    geom_line(aes(x=month,y=rainfall_norm),col="#1582AD",size=0.8)+
    geom_ribbon(aes(x=month,ymin=lower_maxinc,ymax=upper_maxinc),alpha=0.4,fill="#CE5126")+
    geom_line(aes(x=month,y=inc_maxinc),size=0.8,col="#CE5126")+
    facet_geo(~ sites, grid = province_grid%>%
                               select(row,col,code,name))+
    # geom_ribbon(aes(x=month,ymin=lower_maxinc,ymax=upper_maxinc),alpha=0.2,fill="#414487FF")+
    # geom_line(aes(x=month,y=prev_maxinc),size=1,color="#414487FF")+
    scale_x_date(date_labels = "'%y",limits = date_limits)+
    ylab("Normalised incidence/normalised rainfall")+xlab("Year")+
    coord_cartesian(ylim=c(0,1))+
    theme(axis.text.y = element_blank())
  
  inc.rainfall.2 <- ggplot(inc.rainfall.df)+
    geom_line(aes(x=month,y=rainfall_norm*max(inc.df$upper)),col="#1582AD",size=0.8)+
    geom_ribbon(aes(x=month,ymin=lower,ymax=upper),alpha=0.2,fill="#CE5126")+
    geom_line(aes(x=month,y=median),size=0.8,col="#CE5126")+
    facet_geo(~ sites, grid = province_grid%>%
                select(row,col,code,name),scale = 'free_y')+
    scale_x_date(date_labels = "%Y",limits = date_limits)+
    ylab("Estimated incidence/normalised rainfall")+xlab("Year")
  return(list(EIR_SD.density = EIR_SD.density,
              log_init_EIR.density =log_init_EIR.density,
              init_EIR.density = init_EIR.density,
              obs_prev_plot = obs_prev_plot,
              est_prev_plot = est_prev_plot,
              inc_plot = inc_plot,
              inc.rainfall = inc.rainfall,
              inc.rainfall.3 = inc.rainfall.3,
              inc.rainfall.2 = inc.rainfall.2,
              inc_rf_corr = inc_rf_corr
              ))
}
ANC_rainfall_region<-read_csv("./tanz/processed_inputs/ANC_data_Rainfall_ad1.csv")%>%
  mutate(yearmon=as.yearmon(yearmon))
province_grid<-read_csv("./tanz/processed_inputs/Province_grid.csv")
province_grid$name=province_grid$code=province_grid$NAME_1

##Tanzania
tanzania_summary_plots <- create_summary_plots(results = tanz_lt20_2015to2017_results,
                                                 data_list = tanz_data_list_15to17$lt20,
                                                 rainfall = tanz_rainfall,
                                                 level = 'Region')
windows(5,5)
obs <- tanzania_summary_plots$obs_prev_plot
ggsave('tanz/figures/obs_prev_tanz_230222.pdf',plot = obs,width = 8,height=8)
prev <- tanzania_summary_plots$est_prev_plot
ggsave('tanz/figures/est_prev_tanz_230222.pdf',plot = prev,width = 8,height=8)

inc <- tanzania_summary_plots$inc.rainfall
ggsave('tanz/figures/est_inc_tanz_280222.pdf',plot = inc,width = 8,height=8)
tanzania_summary_plots$inc.rainfall.3

## Lake Malawi
lakemalawi_diagnostic_plots <- create_diag_figs_bulk(tanz_lakemal_lt20_2015to2017_results[-11])
windows(10,7)
lakemalawi_diagnostic_plots$ar.plot
lakemalawi_diagnostic_plots$ess.plot
lakemalawi_diagnostic_plots$prior.trace
lakemalawi_diagnostic_plots$posterior.trace
lakemalawi_diagnostic_plots$EIR_SD.density
lakemalawi_diagnostic_plots$init_EIR.density


lakemalawi_summary_plots <- create_summary_plots(results = tanz_lakemal_lt20_2015to2017_results[-11],
                                                 data_list = tanz_lakemal_list_15to17$lt20[-11],
                                                 rainfall = lake_malawi_raw,
                                                 level = 'Council')
lake_obs <- lakemalawi_summary_plots$obs_prev_plot
lake_corr <- lakemalawi_summary_plots$inc_rf_corr
windows(10,10)
tz_regions <- read_sf("./tanz/tza_admbnda_adm1_20181019/tza_admbnda_adm1_20181019.shp")%>%
  mutate(sites=ifelse(ADM1_EN=='Dar-es-salaam','Dar Es Salaam',ADM1_EN))
tz_corr <- merge(tz_regions,inc_rf_corr,by='sites')
tm_shape(tz_corr) + tm_borders() + tm_fill(col='corr') + tmap_options(check.and.fix = TRUE)
library(sf)
ggsave('tanz/figures/obs_prev_lake_230222.pdf',plot = obs,width = 8,height=8)

lakemalawi_summary_plots$inc.rainfall
windows(10,10)
lakemalawi_summary_plots$inc.rainfall.2
windows(10,7)
lakemalawi_summary_plots$inc.rainfall.3
max_month <- tanz_rainfall%>%
  group_by(Region)%>%
  dplyr::summarise(max(yearmon))
max_month_data <- tanz_data_all_2017%>%
  group_by(Region)%>%
  dplyr::summarise(max(as.yearmon(yearmon)))
rainfall <- tanz_rainfall
inc.df <- lapply(1:length(results), 
                           function(x){
                             inc_history <- data.frame(t(results[[x]]$history['inc', 101:1000, -1]))%>%
                               mutate(t=c(1:nrow(inc_history)))%>%
                               melt(id='t')%>%
                               dplyr::rename(time=t)%>%
                               group_by(time)%>%
                               dplyr::summarise(median=median(value),
                                                mean=mean(value),
                                                upper=quantile(value,probs=0.975),
                                                lower=quantile(value,probs=0.025))%>%
                               mutate(sites = names(results[x]),
                                      month = as.Date(data_list[[x]]$month))
                             return(inc_history)
                           })
seas_prop_single_3<-function(time,out_vect){
  if(time<=10) return(sum(out_vect[time:(time+2)])/sum(out_vect))
  else return(sum(c(out_vect[time:12],out_vect[1:(2-(12-time))]))/sum(out_vect))
  }
seas_prop_single_4<-function(time,out_vect){
  if(time<=9) return(sum(out_vect[time:(time+3)])/sum(out_vect))
  else return(sum(c(out_vect[time:12],out_vect[1:(3-(12-time))]))/sum(out_vect))
  }
inc_vect_year <- inc.df[x]$Rainfall
get_max_single_3<-function(inc_vect_year){
  seas_fixed_year=sapply(1:12,seas_prop_single_3,out_vect=inc_vect_year)
  return(c(max(seas_fixed_year),which.max(seas_fixed_year)))}
get_max_single_4<-function(inc_vect_year){
  seas_fixed_year=sapply(1:12,seas_prop_single_4,out_vect=inc_vect_year)
  return(c(max(seas_fixed_year),which.max(seas_fixed_year)))}
##Rolling proportions
level <- 'Region'
prop_summary <- function(results,data_list,rainfall,level=c('Region','District','Council')){
  rainfall$sites <- gsub(' Region','',rainfall$Region)
  rainfall$month <- as.Date(as.yearmon(rainfall$yearmon))
  inc.df <- lapply(1:length(results), 
                   function(x){
                     inc_history <- data.frame(t(results[[x]]$history['inc', 101:1000, -1]))%>%
                       mutate(t=c(1:nrow(inc_history)))%>%
                       melt(id='t')%>%
                       dplyr::rename(time=t)%>%
                       group_by(time)%>%
                       dplyr::summarise(median=median(value),
                                        mean=mean(value),
                                        upper=quantile(value,probs=0.975),
                                        lower=quantile(value,probs=0.025))%>%
                       mutate(sites = names(results[x]),
                              month = as.Date(data_list[[x]]$month))
                     inc_plus_rainfall <- left_join(tail(inc_history,13)[-13,],rainfall,by=c('month','sites'))%>%
                       mutate(rainfall_norm = Rainfall/max(Rainfall),
                              rainfall_maxinc = Rainfall * (max(upper)/max(Rainfall)))
                     return(inc_plus_rainfall)
                   })

  prop.max.df <- bind_rows(lapply(1:length(inc.df), 
                              function(x){
                                rainfall_prop_vals <- get_max_single_3(inc.df[[x]]$Rainfall)
                                inc_prop_vals <- get_max_single_4(inc.df[[x]]$median)
                                month_vect <- month(inc.df[[x]]$month,label=TRUE,abbr=TRUE)
                                data.frame(sites = names(inc.df[x]),
                                           rainfall_prop_max = rainfall_prop_vals[1],
                                           rainfall_prop_month = month_vect[as.integer(rainfall_prop_vals[2])],
                                           inc_prop_max = inc_prop_vals[1],
                                           inc_prop_month = month_vect[as.integer(inc_prop_vals[2])])
                             }))
  prop.time.df <- bind_rows(lapply(1:length(inc.df), 
                                   function(x){
                                     data.frame(sites = names(inc.df[x]),
                                                month = month(inc.df[[x]]$month,label=TRUE,abbr=TRUE),
                                                rainfall_prop_vals = sapply(1:12,seas_prop_single_3,out_vect=inc.df[[x]]$Rainfall),
                                                inc_prop_vals = sapply(1:12,seas_prop_single_4,out_vect=inc.df[[x]]$median))
                                   }))

  inc.rainfall.prop <- ggplot(prop.time.df)+
    geom_line(aes(x=month,y=rainfall_prop_vals,group=sites),col="#2A788EFF",size=0.8)+
    geom_line(aes(x=month,y=inc_prop_vals,group=sites),size=0.8,col="#5DC863FF")+
    facet_geo(~ sites, grid = province_grid%>%
                               select(row,col,code,name))+
    ylab("Proportion of total year")+xlab("Starting month")+
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))+
    scale_x_discrete(breaks = c('Jan','Apr','Jul','Oct'))

scatter <- ggplot(prop.max.df)+
  geom_point(aes(x=rainfall_prop_max,y=inc_prop_max))+
  labs(x='Maximum Rainfall Proportion',y='Maximum Incidence Proportion')+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))
  
}
windows(height=20,width=20)
windows(10,10)
inc.rainfall.prop
scatter

###tanz_all_2015to2022_results
diagnostic_plots <- create_diag_figs_bulk(tanz_all_2015to2022_results)
windows(10,10)
diagnostic_plots$ar.plot
diagnostic_plots$ess.plot
diagnostic_plots$prior.trace
diagnostic_plots$posterior.trace
diagnostic_plots$likelihood.trace
diagnostic_plots$init_EIR.density
diagnostic_plots$EIR_SD.density

ANC_rainfall_region<-read_csv("./tanz/processed_inputs/ANC_data_Rainfall_ad1_2022.csv")%>%
  mutate(yearmon=as.yearmon(yearmon))
province_grid<-read_csv("./tanz/processed_inputs/Province_grid.csv")
province_grid$name=province_grid$code=province_grid$NAME_1

##Tanzania
tanzania_summary_plots <- create_summary_plots(results = tanz_all_2015to2022_results,
                                               data_list = tanz_data_list_15to22,
                                               rainfall = ANC_rainfall_region,
                                               level = 'Region',
                                               date_limits = c(as.Date('2016-01-01'),NA))
windows(10,10)
obs <- tanzania_summary_plots$obs_prev_plot
ggsave('tanz/figures/obs_prev_tanz_300523.pdf',plot = obs,width = 8,height=8)
prev <- tanzania_summary_plots$est_prev_plot
ggsave('tanz/figures/est_prev_tanz_300523.pdf',plot = prev,width = 8,height=8)

inc <- tanzania_summary_plots$inc.rainfall
ggsave('tanz/figures/est_inc_tanz_300523.pdf',plot = inc,width = 8,height=8)
inc_norm <- tanzania_summary_plots$inc.rainfall.3
ggsave('tanz/figures/est_incnorm_tanz_300523.pdf',plot = inc_norm,width = 8,height=8)

##National-level incidence estimates
tz_census <- read_excel('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Tanzania/tz_province_census_2022.xlsx')
tz_main_pop <- tz_census[tz_census$Region=='TanzaniaMainland',]$pop_2022_census
tz_pop_prop <- tz_census %>%
  filter(Region != 'Tanzania' & Region != 'TanzaniaMainland')%>%
  mutate(prop = pop_2022_census/tz_main_pop,
         Region = recode(Region, DaresSalaam = 'Dar Es Salaam',.default = Region))
  

results <- tanz_all_2015to2022_results#[c('Lindi','Mtwara','Ruvuma')]
rainfall <- ANC_rainfall_region
rainfall$sites <- sapply(1:nrow(rainfall), function(x) gsub(' Region','',rainfall[x,level]))
rainfall$month <- as.Date(as.yearmon(rainfall$yearmon))

inc.rainfall.df <- bind_rows(lapply(1:length(results), 
                                    function(x){
                                      history.df <- as.data.frame(t(results[[x]]$history['clininc_all', 101:1000, -1]))
                                      inc_history <- history.df%>%
                                        mutate(t=c(1:nrow(history.df)))%>%
                                        melt(id='t')%>%
                                        dplyr::rename(time=t)%>%
                                        group_by(time)%>%
                                        dplyr::summarise(inc.median=median(value*30),
                                                         inc.mean=mean(value*30),
                                                         inc.upper=quantile(value*30,probs=0.975),
                                                         inc.lower=quantile(value*30,probs=0.025))%>%
                                        mutate(sites = names(results[x]),
                                               month = as.Date(data_list[[x]]$month))
                                      inc_plus_rainfall <- left_join(inc_history,rainfall,by=c('month','sites'))%>%
                                        mutate(rainfall_norm = Rainfall/max(Rainfall),
                                               rainfall_maxinc = Rainfall * (max(inc.median)/max(Rainfall)))
                                      history.df.prev <- as.data.frame(t(results[[x]]$history['prev_05', 101:1000, -1]))
                                      prev_history <- history.df.prev%>%
                                        mutate(t=c(1:nrow(history.df.prev)))%>%
                                        melt(id='t')%>%
                                        dplyr::rename(time=t)%>%
                                        group_by(time)%>%
                                        dplyr::summarise(prev.median=median(value),
                                                         prev.mean=mean(value),
                                                         prev.upper=quantile(value,probs=0.975),
                                                         prev.lower=quantile(value,probs=0.025))%>%
                                        mutate(sites = names(results[x]),
                                               month = as.Date(data_list[[x]]$month))
                                      all <- left_join(inc_plus_rainfall,prev_history,by=c('month','sites'))%>%
                                        mutate(prev_maxinc = prev.median * (max(inc.median)),
                                               upper_maxinc = prev.upper *max(inc.median),
                                               lower_maxinc = prev.lower *max(inc.median),
                                               rainfall_maxprev = Rainfall * (max(prev.median)/max(Rainfall)),
                                               inc_maxinc = inc.median/ max(inc.median),
                                               upper_maxinc = inc.upper /max(inc.median),
                                               lower_maxinc = inc.lower /max(inc.median),
                                               rainfall_maxprev = Rainfall * (max(prev.median)/max(Rainfall)))
                                      return(all)
                                    }))
inc.rainfall.df <- addCIs(inc.rainfall.df,inc.rainfall.df$positive,inc.rainfall.df$tested)
est_cases <- left_join(inc.rainfall.df,tz_pop_prop,by='Region')%>%
  mutate(est_cases.median = inc.median*pop_2022_census,
         est_cases.mean = inc.mean*pop_2022_census,
         est_cases.upper = inc.upper*pop_2022_census,
         est_cases.lower = inc.lower*pop_2022_census
         )
est_cases_4PW <- est_cases %>%
  select(inc.median,inc.mean,inc.lower,inc.upper,sites,yearmon,positive,tested,total,
         pop_2012_census,pop_2022_census)
saveRDS(est_cases_4PW,'./tanz/est_cases_tz.rds')
windows(10,10)
scaler <- 400000
est_cases.rainfall <- ggplot(est_cases)+
  geom_line(aes(x=month,y=rainfall_norm*scaler),col="#2A788EFF",size=0.8)+
  geom_ribbon(aes(x=month,ymin=est_cases.lower,ymax=est_cases.upper),alpha=0.4,fill="#D95F02")+
  geom_line(aes(x=month,y=est_cases.median),size=0.8,col="#D95F02")+
  scale_x_date(date_labels = "'%y",limits = date_limits)+
  facet_geo(~ sites, grid = province_grid%>%
              select(row,col,code,name))+
  ylab("Estimated case counts/normalised rainfall")+xlab("Year")+
  coord_cartesian(ylim=c(0,scaler))

num_tested_plot <- ggplot(est_cases)+
  geom_line(aes(x=month,y=tested),size=0.8,col="#999999")+
  scale_x_date(date_labels = "'%y",limits = date_limits)+
  scale_y_continuous(limits = c(0,NA))+
  facet_geo(~ sites, grid = province_grid%>%
              select(row,col,code,name))+
  ylab("Number tested at first ANC")+xlab("Year")
rate_tested_plot <- ggplot(est_cases)+
  geom_line(aes(x=month,y=(tested*100*12/pop_2022_census)),size=0.8,col="#999999")+
  scale_x_date(date_labels = "'%y",limits = date_limits)+
  scale_y_continuous(limits = c(0,NA))+
  facet_geo(~ sites, grid = province_grid%>%
              select(row,col,code,name))+
  ylab("Percent of population tested\nper year")+xlab("Year")


  coord_cartesian(ylim=c(0,scaler))

est_cases_total <- est_cases %>%
  group_by(month)%>%
  summarise(est_cases.median=sum(est_cases.median),
            est_cases.mean=sum(est_cases.mean),
            est_cases.lower=sum(est_cases.lower),
            est_cases.upper=sum(est_cases.upper))
est_cases_total_plot <- ggplot(est_cases_total)+
  geom_ribbon(aes(x=month,ymin=est_cases.lower,ymax=est_cases.upper),alpha=0.4,fill="#D95F02")+
  geom_line(aes(x=month,y=est_cases.median),size=0.8,col="#D95F02")+
  scale_x_date(date_labels = "'%y",limits = date_limits)+
  ylab("Estimated case counts")+xlab("Year")+
  coord_cartesian(ylim=c(0,3000000))

est_cases_total_annual <- est_cases_total %>%
  mutate(year=year(month))%>%
  group_by(year)%>%
  summarise(est_cases.median=sum(est_cases.median),
            est_cases.mean=sum(est_cases.mean),
            est_cases.lower=sum(est_cases.lower),
            est_cases.upper=sum(est_cases.upper))
est_cases_total_annual_plot <- ggplot(est_cases_total_annual)+
  geom_ribbon(aes(x=as.Date(paste0(year,'-01-01')),ymin=est_cases.lower,ymax=est_cases.upper),alpha=0.4,fill="#D95F02")+
  geom_line(aes(x=as.Date(paste0(year,'-01-01')),y=est_cases.median),size=0.8,col="#D95F02")+
  scale_x_date(date_labels = "'%y",limits = date_limits)+
  ylab("Estimated case counts")+xlab("Year")+
  coord_cartesian(ylim=c(0,NA))

##Lindi and Mtwara illustration
lindi_mtwara_ruvuma<-Full_data%>%
  mutate(yearmon=as.yearmon(paste0(Full_data$Year,"-",Full_data$Month),"%Y-%m"))%>%
  filter(!is.na(ANC_test)&ANC_re_ad!=0)%>%
  group_by(yearmon,Region)%>%
  dplyr::summarise(count=n(),
                   positive=sum(ANC_pos,na.rm = TRUE),
                   tested=sum(ANC_test,na.rm = TRUE),
                   total=sum(Total_re_ad,na.rm = TRUE),
                   first_visit_total=sum(Before_12wk,na.rm = TRUE)+sum(After_12wk,na.rm = TRUE),
                   first_visit_b12wk=sum(Before_12wk,na.rm = TRUE),
                   first_visit_a12wk=sum(After_12wk,na.rm = TRUE))%>%
  filter(Region == 'Lindi Region' | Region == 'Mtwara Region' | Region == 'Ruvuma Region') %>%
  mutate(period = ifelse(yearmon>=as.yearmon('2018-01-01')&yearmon<=as.yearmon('2019-12-31'),'2018-2019',
                         ifelse(yearmon>=as.yearmon('2020-01-01')&yearmon<=as.yearmon('2021-12-31'),'2020-2021',NA)))
windows(5,10)
ylim.prim <- c(0, 5000)   # school prev limits
ylim.sec <- c(0, 500)    # ANC prev limits

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]
RColorBrewer::display.brewer.pal(8,'Paired')
RColorBrewer::brewer.pal(8,'Paired')
colors <- c('Women tested' = "#33A02C", 'Women attended (all visits)' = "#B2DF8A",
            'Facilities reporting' = "#1F78B4")
colors_2 <- c('Women tested' = "#33A02C", 'Women attended (all visits)' = "#B2DF8A",
              'Women attended (first visits)' = "#B2DF8A",
              'Facilities reporting' = "#1F78B4",'Normalised rainfall' = "#1582AD", 'Observed prevalence' = '#6D6A67',
              'Fitted prevalence' = "#EFBB12", 'Estimated incidence' = "#CE5126")
obs_tested_plot_dist <- ggplot(lindi_mtwara_ruvuma)+
  geom_line(aes(x=as.Date(as.yearmon(yearmon)),y=tested,color='Women tested'),size=1)+
  # geom_line(aes(x=as.Date(as.yearmon(yearmon)),y=total,color='Women attended (all visits)'),size=1)+
  geom_line(aes(x=as.Date(as.yearmon(yearmon)),y=first_visit_total,color='Women attended (first visits)'),size=1)+
  geom_line(aes(x=as.Date(as.yearmon(yearmon)),y=count*b+a,color='Facilities reporting'),size=1)+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_wrap(~ Region,nrow = 3)+
  scale_color_manual(values=colors_2)+
  scale_x_date(date_labels = "'%y")+
  scale_y_continuous("# Women Attended and Tested", sec.axis = sec_axis(~ (. - a)/b, name = "# Facilities Reporting"),limits = c(0,NA)) +
  labs(x='Date')+
  theme(legend.title = element_blank())

windows(12,12)
obs_testedprop_plot_dist <- ggplot(lindi_mtwara)+
  geom_line(aes(x=as.Date(as.yearmon(yearmon)),y=tested/first_visit_total),size=1)+
  geom_hline(yintercept = 1, linetype='dashed',color='#666666')+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_wrap(~ Region,nrow = 3)+
  scale_x_date(date_labels = "'%y")+
  scale_y_continuous("Proportion of women tested") +
  labs(x='Date')

toomanytests <- Full_data%>%
  mutate(yearmon=as.yearmon(paste0(Full_data$Year,"-",Full_data$Month),"%Y-%m"))%>%
  filter(!is.na(ANC_test)&ANC_re_ad!=0)%>%
  filter(ANC_test>(Before_12wk+After_12wk))
  
table(toomanytests$yearmon)
obs_prev_plot <- ggplot(inc.rainfall.df[!(inc.rainfall.df$sites=='Mbeya'&inc.rainfall.df$yearmon=='Oct 2020'),])+
  geom_line(aes(x=month,y=rainfall_norm*0.3,color='Normalised rainfall'),size=0.8)+
  # geom_ribbon(aes(x=month,ymin=lower,ymax=upper),alpha=0.2)+
  # geom_line(aes(x=month,y=median),size=1)+
  # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  geom_point(aes(x=month,y=mean,color='Observed prevalence'),pch = 19)+
  geom_errorbar(aes(x=month,ymin=lower,ymax=upper,color='Observed prevalence'),width = 0)+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_wrap(~ sites,nrow=3)+
  scale_x_date(date_labels = "'%y",limits = date_limits)+
  scale_color_manual(values=colors_2)+
  labs(x='Date',y='ANC Prevalence')+
  theme(legend.title = element_blank())

est_prev_plot <- ggplot(inc.rainfall.df[!(inc.rainfall.df$sites=='Mbeya'&inc.rainfall.df$yearmon=='Oct 2020'),])+
  geom_line(aes(x=month,y=rainfall_norm*0.3,color='Normalised rainfall'),col="#1582AD",size=0.8)+
  geom_point(aes(x=month,y=mean,color='Observed prevalence'),pch = 19)+
  geom_errorbar(aes(x=month,ymin=lower,ymax=upper,color='Observed prevalence'),width = 0)+
  geom_ribbon(aes(x=month,ymin=prev.lower,ymax=prev.upper),alpha=0.2,fill="#EFBB12")+
  geom_line(aes(x=month,y=prev.median,color='Fitted prevalence'),size=1)+
  # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_wrap(~ sites,nrow=3)+
  scale_x_date(date_labels = "'%y",limits = date_limits)+
  scale_color_manual(values=colors_2)+
  labs(x='Date',y='ANC Prevalence')+
  theme(legend.title = element_blank())

scaler <- 0.0075*30
inc.rainfall <- ggplot(inc.rainfall.df)+
  geom_line(aes(x=month,y=rainfall_norm*scaler,color='Normalised rainfall'),size=0.8)+
  geom_ribbon(aes(x=month,ymin=inc.lower,ymax=inc.upper),alpha=0.4,fill="#CE5126")+
  geom_line(aes(x=month,y=inc.median,color='Estimated incidence'),size=0.8)+
  scale_x_date(date_labels = "'%y",limits = date_limits)+
  scale_color_manual(values=colors_2)+
  facet_wrap(~ sites,nrow=3)+
  ylab("Estimated incidence/normalised rainfall")+xlab("Year")+
  theme(legend.title = element_blank())

windows(12,8)
obs_tested_plot_dist + obs_prev_plot + est_prev_plot + inc.rainfall + plot_layout(ncol=4,guides = 'collect') &
  guides(color=guide_legend(nrow=3, byrow=TRUE))
dash <- obs_prev_plot + est_prev_plot + inc.rainfall + plot_layout(ncol=3,guides = 'collect') &
  guides(color=guide_legend(nrow=3, byrow=TRUE))
ggsave('./tanz/figures/lindi_mtwara_ruvuma_model_300523.pdf',plot=dash,width = 12,height=8)
ggarrange(obs_tested_plot_dist, obs_prev_plot, est_prev_plot, inc.rainfall, ncol=4, nrow=1, common.legend = TRUE, legend="right")
