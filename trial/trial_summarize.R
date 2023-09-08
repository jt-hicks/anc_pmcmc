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
library(matrixStats)
library(epitools)
library(coda)

theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))

create_diag_figs <- function(result,country,district){
  print('acceptance rate')
  ar <- 1 - coda::rejectionRate(as.mcmc(result$mcmc))
  print(ar)
  print('effective size')
  ess <- coda::effectiveSize(as.mcmc(result$mcmc))
  print(ess)
  
  title <- paste0('Diagnostic plots for seasonal model - ',district,', ',country)
  pars_list <- c('log_prior','log_likelihood','log_posterior','volatility','log_init_EIR')
  trace_plots <- lapply(pars_list, function(x){
    bayesplot::mcmc_trace(result$mcmc[51:1000,],pars = x) + 
      ggtitle(paste0(x,' / AR: ',round(ar[x],3),' / ESS: ',round(ess[x],1)))+
      theme(title = element_text(size=6),
            axis.title.y = element_blank())
  })
  dense_plots <- lapply(pars_list, function(x){
    bayesplot::mcmc_dens(result$mcmc[51:1000,],pars = x) + 
      ggtitle(paste0(x,' / AR: ',round(ar[x],2),' / ESS: ',round(ess[x],1)))+
      theme(title = element_text(size=6),
            axis.title.x = element_blank())
  })
  diag <- (trace_plots[[1]]+dense_plots[[1]])/
    (trace_plots[[2]]+dense_plots[[2]])/
    (trace_plots[[3]]+dense_plots[[3]])/
    (trace_plots[[4]]+dense_plots[[4]])/
    (trace_plots[[5]]+dense_plots[[5]]) +
    plot_layout(guides = "collect") + plot_annotation(title = title)
  
  
  return(diag)
}

diag_wa_pg_seas_plots <- lapply(1:4,function(i) create_diag_figs(wa_pg_seas_result_list[[i]],country = '',district = names(WA_pg_data_list[i])))
diag_wa_pg_seas_2_plots <- lapply(1:4,function(i) create_diag_figs(wa_pg_seas_2_result_list[[i]],country = '',district = names(WA_pg_data_list[i])))
diag_wa_pg_std_plots <- lapply(1:4,function(i) create_diag_figs(wa_pg_std_result_list[[i]],country = '',district = names(WA_pg_data_list[i])))
diag_wa_pg_seas_padded_plots <- lapply(1:4,function(i) create_diag_figs(wa_pg_seas_padded_result_list[[i]],country = '',district = names(WA_pg_data_list[i])))
diag_gamb_pg_seas_plots <- create_diag_figs(wa_pg_gamb_seas$result(),country = '',district = 'Gambia')
diag_wa_pg_bulk_seas_090623_plots <- lapply(1:4,function(i) create_diag_figs(wa_pg_bulk_seas_090623_result_list[[i]],country = '',district = names(WA_pg_data_list[i])))
diag_wa_all_bulk_seas_090623_plots <- lapply(1:4,function(i) create_diag_figs(wa_all_bulk_seas_090623_result_list[[i]],country = '',district = names(WA_pg_data_list[i])))
diag_wa_pgsg_bulk_sifter_betaa_plots <- lapply(1:4,function(i) create_diag_figs(wa_pgsg_bulk_sifter_betaa_results[[i]],country = '',district = names(WA_pg_data_list[i])))
diag_wa_pgsg_bulk_sifter_betaa_2407_plots <- lapply(1:4,function(i) create_diag_figs(wa_pgsg_bulk_sifter_betaa_2407_results[[i]],country = '',district = names(WA_pg_data_list[i])))


str(wa_pg_bulk_seas_090623_result_list[[1]]$mcmc)
windows(15,15)
diag_wa_pg_seas_plots[1]
diag_wa_pg_seas_plots[2]
diag_wa_pg_seas_plots[3]
diag_wa_pg_seas_plots[4]
diag_wa_pg_seas_2_plots[1]
diag_wa_pg_seas_2_plots[2]
diag_wa_pg_seas_2_plots[3]
diag_wa_pg_seas_2_plots[4]
diag_wa_pg_std_plots[1]
diag_wa_pg_std_plots[2]
diag_wa_pg_std_plots[3]
diag_wa_pg_std_plots[4]
diag_wa_pg_seas_padded_plots[1]
diag_wa_pg_seas_padded_plots[2]
diag_wa_pg_seas_padded_plots[3]
diag_wa_pg_seas_padded_plots[4]
diag_wa_pg_bulk_seas_090623_plots[1]
diag_wa_pg_bulk_seas_090623_plots[2]
diag_wa_pg_bulk_seas_090623_plots[3]
diag_wa_pg_bulk_seas_090623_plots[4]
diag_wa_all_bulk_seas_090623_plots[1]
diag_wa_all_bulk_seas_090623_plots[2]
diag_wa_all_bulk_seas_090623_plots[3]
diag_wa_all_bulk_seas_090623_plots[4]
windows(15,15)
diag_wa_pgsg_bulk_sifter_betaa_plots[1]
windows(15,15)
diag_wa_pgsg_bulk_sifter_betaa_plots[2]
windows(15,15)
diag_wa_pgsg_bulk_sifter_betaa_plots[3]
windows(15,15)
diag_wa_pgsg_bulk_sifter_betaa_plots[4]
diag_wa_pgsg_bulk_sifter_betaa_plots
diag_wa_pgsg_bulk_sifter_betaa_max125_plots <- lapply(1:4,function(i) create_diag_figs(wa_pgsg_bulk_sifter_betaa_max125_results[[i]],country = '',district = names(WA_pg_data_list[i])))
for(i in 1:4){
  ggsave(filename=gsub(' ','',paste0('./trial/figures/Diagnostic/diag_wa_pgsg_bulk_sifter_betaa_max125_',names(wa_pgsg_bulk_sifter_betaa_max125_results[i]),'.pdf')),
         plot=diag_wa_pgsg_bulk_sifter_betaa_max125_plots[[i]],
         width=10,height=10,units='in')
}
for(i in 1:4){
  ggsave(filename=gsub(' ','',paste0('./trial/figures/Diagnostic/diag_wa_pgsg_bulk_sifter_betaa_2407_',names(wa_pgsg_bulk_sifter_betaa_2407_results[i]),'.pdf')),
         plot=diag_wa_pgsg_bulk_sifter_betaa_2407_plots[[i]],
         width=10,height=10,units='in')
}
diag_ea_pgmg_bulk_sifter_eir_max2000_plots <- lapply(1:2,function(i) create_diag_figs(ea_pgmg_bulk_sifter_eir_max2000_results[[i]],country = '',district = names(EA_pg_data_list[i])))
for(i in 1:2){
  ggsave(filename=gsub(' ','',paste0('./trial/figures/Diagnostic/diag_ea_pgmg_bulk_sifter_eir_max2000_',names(ea_pgmg_bulk_sifter_eir_max2000_results[i]),'.pdf')),
         plot=diag_ea_pgmg_bulk_sifter_eir_max2000_plots[[i]],
         width=10,height=10,units='in')
}

diag_ea_pgmg_bulk_sifter_betaa_max125_plots <- lapply(1:2,function(i) create_diag_figs(ea_pgmg_bulk_sifter_betaa_max125_results[[i]],country = '',district = names(EA_pg_data_list[i])))
for(i in 1:2){
  ggsave(filename=gsub(' ','',paste0('./trial/figures/Diagnostic/diag_ea_pgmg_bulk_sifter_betaa_max125_',names(ea_pgmg_bulk_sifter_betaa_max125_results[i]),'.pdf')),
         plot=diag_ea_pgmg_bulk_sifter_betaa_max125_plots[[i]],
         width=10,height=10,units='in')
}
diag_ea_pgmg_bulk_sifter_betaa_2407_plots <- lapply(1:2,function(i) create_diag_figs(ea_pgmg_bulk_sifter_betaa_2407_results[[i]],country = '',district = names(EA_pg_data_list[i])))
for(i in 1:2){
  ggsave(filename=gsub(' ','',paste0('./trial/figures/Diagnostic/diag_ea_pgmg_bulk_sifter_betaa_2407_',names(ea_pgmg_bulk_sifter_betaa_2407_results[i]),'.pdf')),
         plot=diag_ea_pgmg_bulk_sifter_betaa_2407_plots[[i]],
         width=10,height=10,units='in')
}

##Create dashboard by country
results <- wa_pg_bulk_seas_090623_result_list
prev_pg <- WA_pg_data_list
prev_mg <- WA_sg_data_list

create_dashboard_plots <- function(results,prev_pg=WA_pg_data_list,prev_mg=WA_sg_data_list){
  coefs_pg_df <- as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS'))
  coefs_mg_df <- as.data.frame(readRDS('./nnp/Corr/sg_corr_sample.RDS'))
  
  ## fn to return prevalence from log_odds
  get_prev_from_log_odds<-function(log_odds){
    return(exp(log_odds)/(1+exp(log_odds)))
  }
  
  ## fn to return odds from prevalence
  get_odds_from_prev<-function(prev){
    return(prev/(1-prev))
  }
  
  districts <- c("Burkina Faso","Gambia","Ghana","Mali" )
  district_labels <- c("Burkina Faso","Gambia","Ghana","Mali")
  colors <- c(`Burkina Faso` = "#1B9E77", Gambia = "#999999", Ghana = "#D95F02", Mali = "#377EB8")
  
  
  df_prev <- data.frame(time = integer(),
                        median = numeric(),
                        mean = numeric(),
                        upper = numeric(),
                        lower = numeric(),
                        site = character(),
                        month = character())
  df_prev_sample <- data.frame(time = integer(),
                               value = numeric(),
                               variable = character(),
                               site = character(),
                               month = character())
  
  for(i in 1:4){
    prev_history <- data.frame(t(results[[i]]$history['prev_05', 51:1000, -1]))
    dates <- WA_pg_data_list[[i]]$month
    
    long_prev_sum <- prev_history%>%
      mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      mutate(logodds_child = log(get_odds_from_prev(value)),
             prev_pg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
             }))),
             prev_mg = rowMedians(as.matrix(sapply(1:nrow(coefs_mg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_mg_df$gradient[x]*(logodds_child-coefs_mg_df$av_lo_child[x])+coefs_mg_df$intercept[x])
             }))))%>%
      group_by(time)%>%
      summarise(median_pg=median(prev_pg),
                mean_pg=mean(prev_pg),
                upper_pg=quantile(prev_pg,probs=0.975),
                lower_pg=quantile(prev_pg,probs=0.025),
                median_mg=median(prev_mg),
                mean_mg=mean(prev_mg),
                upper_mg=quantile(prev_mg,probs=0.975),
                lower_mg=quantile(prev_mg,probs=0.025))%>%
      mutate(site = districts[[i]],
             month = as.yearmon(dates))
    df_prev <- rbind(df_prev,long_prev_sum)
    
    prev_sample <- prev_history[, sample(ncol(prev_history), 100)] %>%
      mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      mutate(logodds_child = log(get_odds_from_prev(value)),
             prev_pg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
             }))),
             prev_mg = rowMedians(as.matrix(sapply(1:nrow(coefs_mg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_mg_df$gradient[x]*(logodds_child-coefs_mg_df$av_lo_child[x])+coefs_mg_df$intercept[x])
             }))))%>%
      mutate(site = districts[[i]],
             month = as.yearmon(rep(dates,100)))
    df_prev_sample <- rbind(df_prev_sample,prev_sample)
    
    
  }
  
  df_data_pg <- data.frame(t=numeric(),
                           tested=integer(),
                           positive=numeric(),
                           mean=numeric(),
                           upper=numeric(),
                           lower=numeric(),
                           site=character())
  
  for(i in 1:4){
    dates <- WA_pg_data_list[[i]]$month
    data_cis <- addCIs(prev_pg[[i]],prev_pg[[i]]$positive,prev_pg[[i]]$tested)%>%
      mutate(site=names(prev_pg[i]),
             month = dates)
    df_data_pg <- rbind(df_data_pg,data_cis)
  }
  
  df_data_mg <- data.frame(t=numeric(),
                           tested=integer(),
                           positive=numeric(),
                           mean=numeric(),
                           upper=numeric(),
                           lower=numeric(),
                           site=character())
  
  for(i in 1:4){
    dates <- WA_sg_data_list[[i]]$month
    data_cis <- addCIs(prev_mg[[i]],prev_mg[[i]]$positive,prev_mg[[i]]$tested)%>%
      mutate(site=names(prev_pg[i]),
             month = dates)
    df_data_mg <- rbind(df_data_mg,data_cis)
  }
  
  # df_data_all <- data.frame(t=numeric(),
  #                          tested=integer(),
  #                          positive=numeric(),
  #                          mean=numeric(),
  #                          upper=numeric(),
  #                          lower=numeric(),
  #                          district=character())
  
  # for(i in start:(start+number-1)){
  #   data_cis <- addCIs(prev_all[[i]],prev_all[[i]]$positive,prev_all[[i]]$tested)%>%
  #     mutate(district=names(prev_all[i]),
  #            month = dates_list[[country]])
  #   df_data_all <- rbind(df_data_all,data_cis)
  # }
  
  df_inc <- data.frame(time = integer(),
                       median = numeric(),
                       mean = numeric(),
                       upper = numeric(),
                       lower = numeric(),
                       site = character(),
                       month = character())
  df_inc_sample <- data.frame(time = integer(),
                              value = numeric(),
                              variable = character(),
                              site = character(),
                              month = character())
  for(i in 1:4){
    inc_history <- data.frame(t(results[[i]]$history['clininc_all', 51:1000, -1]))
    dates <- WA_pg_data_list[[i]]$month
    
    long_inc_sum <- inc_history%>%
      dplyr::mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value)*10000*30,
                       mean=mean(value)*10000*30,
                       upper=quantile(value,probs=0.975)*10000*30,
                       lower=quantile(value,probs=0.025)*10000*30)%>%
      mutate(site = districts[[i]],
             month = as.yearmon(dates))
    df_inc <- rbind(df_inc,long_inc_sum)
    
    inc_sample <- inc_history[, sample(ncol(inc_history), 100)] %>%
      mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value*10000*30,
                    site = districts[[i]],
                    month = as.yearmon(rep(dates,100)))
    df_inc_sample <- rbind(df_inc_sample,inc_sample)
  }

  df_eir <- data.frame(time = integer(),
                       median = numeric(),
                       mean = numeric(),
                       upper = numeric(),
                       lower = numeric(),
                       site = character(),
                       month = character(),
                       rel_eir = numeric())
  df_eir_sample <- data.frame(time = integer(),
                              value = numeric(),
                              variable = character(),
                              site = character(),
                              month = character(),
                              rel_eir = numeric())
  for(i in 1:4){
    eir_history <- data.frame(t(results[[i]]$history['EIR', 51:1000, -1]))
    dates <- WA_pg_data_list[[i]]$month
    
    long_eir_sum <- eir_history%>%
      dplyr::mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value),
                       mean=mean(value),
                       upper=quantile(value,probs=0.975),
                       lower=quantile(value,probs=0.025))%>%
      mutate(site = districts[[i]],
             month = as.yearmon(dates),
             rel_eir = median/max(median))
    df_eir <- rbind(df_eir,long_eir_sum)
    
    eir_sample <- eir_history[, sample(ncol(eir_history), 100)] %>%
      mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value,
                    site = districts[[i]],
                    month = as.yearmon(rep(dates,100)),
                    rel_eir = value/max(long_eir_sum$median))
    df_eir_sample <- rbind(df_eir_sample,eir_sample)
  }
  df_betaa <- data.frame(time = integer(),
                       median = numeric(),
                       mean = numeric(),
                       upper = numeric(),
                       lower = numeric(),
                       site = character(),
                       month = character(),
                       rel_betaa = numeric())
  df_betaa_sample <- data.frame(time = integer(),
                              value = numeric(),
                              variable = character(),
                              site = character(),
                              month = character(),
                              rel_betaa = numeric())
  for(i in 1:4){
    betaa_history <- data.frame(t(results[[i]]$history['betaa', 51:1000, -1]))
    dates <- WA_pg_data_list[[i]]$month
    
    long_betaa_sum <- betaa_history%>%
      dplyr::mutate(t=c(1:nrow(betaa_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value),
                       mean=mean(value),
                       upper=quantile(value,probs=0.975),
                       lower=quantile(value,probs=0.025))%>%
      mutate(site = districts[[i]],
             month = as.yearmon(dates),
             rel_betaa = median/max(median))
    df_betaa <- rbind(df_betaa,long_betaa_sum)
    
    betaa_sample <- betaa_history[, sample(ncol(betaa_history), 100)] %>%
      mutate(t=c(1:nrow(betaa_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value,
                    site = districts[[i]],
                    month = as.yearmon(rep(dates,100)),
                    rel_betaa = value/max(long_betaa_sum$median))
    df_betaa_sample <- rbind(df_betaa_sample,betaa_sample)
  }
  max_eir <- df_eir_sample %>%
    group_by(site) %>%
    summarise(max_eir = ifelse(max(value)==1000,800,max(value)*0.5))

  rainfall <- read_csv('./trial/Data/WA_ISTp_rainfall.csv') %>%
    group_by(Country) %>%
    filter(between(as.yearmon(Month), min(df_data_pg$month), max(df_data_pg$month)))%>%
    left_join(max_eir, by = join_by(Country == site))%>%
    mutate(rel_rainfall = Rainfall/max(Rainfall))%>%
    dplyr::rename(site = Country)
  
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  # est_inc_plot <- ggplot()+
  #   geom_line(data=df_inc_sample,aes(x=as.Date(month),y=value,color=site,group=variable),alpha=0.1)+
  #   # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  #   geom_line(data=df_inc,aes(x=as.Date(month),y=median,color=site,group=site),size=1)+
  #   # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
  #   scale_color_manual(values=colors)+
  #   scale_fill_manual(values=colors)+
  #   scale_x_date(date_labels = "%b %Y")+
  #   # coord_cartesian(ylim=c(0, 1000))+
  #   coord_cartesian(ylim=c(0, max(df_inc$median)*2))+
  #   # scale_y_continuous(limits=c(0,500))+
  #   # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  #   # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
  #   # labs(x='Date',y='EIR')+
  #   facet_grid(factor(site,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
  #   labs(title = 'Estimated Incidence\nper 10,000 person-mos.')+
  #   theme(legend.title = element_blank(),
  #         legend.position = 'none',
  #         axis.title.x = element_blank(),
  #         axis.title.y = element_blank(),
  #         axis.text.x=element_text(angle=45, hjust=1, vjust=1),
  #         axis.ticks.x = element_line(size = 0.5),
  #         axis.ticks.length = unit(3, "pt")
  #   )
  est_eir_plot <- ggplot()+
    # geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
    geom_line(data=df_eir_sample,aes(x=as.Date(month),y=value,color=site,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_eir,aes(x=as.Date(month),y=median,color=site,group=site),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # coord_cartesian(ylim=c(0, 1000))+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.1,1,10,100,1000),labels=c(.1,1,10,100,1000))+
    coord_cartesian(ylim=c(.1, 800))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    facet_grid(factor(site,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.,scales = 'free_y')+
    labs(title = 'EIR\nInfectious bites/person/year')+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  rel_eir_plot <- ggplot()+
    geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
    # geom_line(data=df_eir_sample,aes(x=as.Date(month),y=rel_eir,color=site,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_eir,aes(x=as.Date(month),y=rel_eir,color=site,group=site),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # coord_cartesian(ylim=c(0, 1000))+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.1,1,10,100,1000),labels=c(.1,1,10,100,1000))+
    coord_cartesian(ylim=c(0, 1.2))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    facet_grid(factor(site,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    labs(title = 'Relative EIR and\nRelative rainfall')+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  # print('est_inc_plot')
  obs_prev_plot_mg <- ggplot()+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=site,group=site),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=site,group=site),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data_mg,aes(x=as.Date(month),y=mean,color=site,group=site),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_mg,aes(x=as.Date(month),ymin=lower,ymax=upper,color=site,group=site),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    facet_grid(factor(site,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    scale_y_continuous(limits = c(0,1))+
    coord_cartesian(ylim=c(0, 1),xlim = range(as.Date(df_data_mg$month)))+
    labs(title = 'ANC Prevalence\nSecundigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  # print('obs_prev_plot_mg')
  # # obs_prev_plot_all <- annotations[[country]]+
  # #   # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  # #   # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
  # #   # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # #   geom_point(data=df_data_all,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
  # #   geom_errorbar(data=df_data_all,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
  # #   scale_color_manual(values=colors)+
  # #   scale_fill_manual(values=colors)+
  # #   scale_x_date(date_labels = "%b %Y")+
  # #   facet_grid(district~.)+
  # #   scale_y_continuous(limits = c(0,1))+
  # #   coord_cartesian(ylim=c(0, 1))+
  # #   labs(title = 'ANC Prevalence\nAll gravidities')+
  # #   theme(legend.title = element_blank(),
  # #         axis.title.x = element_blank(),
  # #         axis.title.y = element_blank(),
  # #         axis.text.x=element_text(angle=45, hjust=1, vjust=1),
  # #         axis.ticks.x = element_line(linewidth = 0.5),
  # #         axis.ticks.length = unit(3, "pt"),
  # #         legend.position = 'none'
  # #   )
  est_prev_pg_plot <- ggplot()+
    geom_line(data=df_prev_sample,aes(x=as.Date(month),y=prev_pg,color=site,group=variable),alpha=0.2)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=site,group=site),alpha=0.2)+
    geom_line(data=df_prev,aes(x=as.Date(month),y=median_pg,color=site,group=site),linewidth=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data_pg,aes(x=as.Date(month),y=mean,color=site,group=site),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_pg,aes(x=as.Date(month),ymin=lower,ymax=upper,color=site,group=site),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits=c(0,1))+
    coord_cartesian(xlim = range(as.Date(df_data_pg$month)))+
    facet_grid(factor(site,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    labs(title = 'ANC Prevalence\nPrimigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  est_prev_mg_plot <- ggplot()+
    geom_line(data=df_prev_sample,aes(x=as.Date(month),y=prev_mg,color=site,group=variable),alpha=0.2)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=site,group=site),alpha=0.2)+
    geom_line(data=df_prev,aes(x=as.Date(month),y=median_mg,color=site,group=site),linewidth=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data_mg,aes(x=as.Date(month),y=mean,color=site,group=site),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_mg,aes(x=as.Date(month),ymin=lower,ymax=upper,color=site,group=site),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits=c(0,1))+
    coord_cartesian(xlim = range(as.Date(df_data_pg$month)))+
    facet_grid(factor(site,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    labs(title = 'ANC Prevalence\nSecundigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  # print('est_prev_plot')
  sample_size_plot <- ggplot()+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=site,group=site),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=site,group=site),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_line(data=df_data_pg,aes(x=as.Date(month),y=tested,color=site,group=site))+
    geom_line(data=df_data_mg,aes(x=as.Date(month),y=tested,color=site,group=site),linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    facet_grid(factor(site,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    labs(title = 'Number Tested\nPrimigrav (solid)\nSecundigrav (dashed)')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  # sample_size_plot+est_prev_plot+obs_prev_plot_mg+est_eir_plot+est_inc_plot+ plot_layout(guides = "collect",ncol=5)
  sample_size_plot+est_prev_pg_plot+est_prev_mg_plot+est_eir_plot+rel_eir_plot+ plot_layout(guides = "collect",ncol=5)
}

create_dashboard_plots_all <- function(results,prev=WA_all_data_list){
  coefs_pg_df <- as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS'))
  ## fn to return prevalence from log_odds
  get_prev_from_log_odds<-function(log_odds){
    return(exp(log_odds)/(1+exp(log_odds)))
  }
  
  ## fn to return odds from prevalence
  get_odds_from_prev<-function(prev){
    return(prev/(1-prev))
  }
  
  districts <- c("Burkina Faso","Gambia","Ghana","Mali" )
  district_labels <- c("Burkina Faso","Gambia","Ghana","Mali")
  colors <- c(`Burkina Faso` = "#1B9E77", Gambia = "#999999", Ghana = "#D95F02", Mali = "#377EB8")
  
  
  df_prev <- data.frame(time = integer(),
                        median = numeric(),
                        mean = numeric(),
                        upper = numeric(),
                        lower = numeric(),
                        site = character(),
                        month = character())
  df_prev_sample <- data.frame(time = integer(),
                               value = numeric(),
                               variable = character(),
                               site = character(),
                               month = character())
  
  for(i in 1:4){
    prev_history <- data.frame(t(results[[i]]$history['prev', 51:1000, -1]))
    dates <- prev[[i]]$month
    
    long_prev_sum <- prev_history%>%
      mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      mutate(logodds_child = log(get_odds_from_prev(value)),
             prev_pg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
             }))))%>%
      group_by(time)%>%
      summarise(median=median(prev_pg),
                mean=mean(prev_pg),
                upper=quantile(prev_pg,probs=0.975),
                lower=quantile(prev_pg,probs=0.025))%>%
      mutate(site = districts[[i]],
             month = as.yearmon(dates))
    df_prev <- rbind(df_prev,long_prev_sum)
    
    prev_sample <- prev_history[, sample(ncol(prev_history), 100)] %>%
      mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      mutate(logodds_child = log(get_odds_from_prev(value)),
             prev_pg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
             }))))%>%
      mutate(site = districts[[i]],
             month = as.yearmon(rep(dates,100)))
    df_prev_sample <- rbind(df_prev_sample,prev_sample)
  }
  
  df_data_all <- data.frame(t=numeric(),
                           tested=integer(),
                           positive=numeric(),
                           mean=numeric(),
                           upper=numeric(),
                           lower=numeric(),
                           site=character())
  
  for(i in 1:4){
    dates <- prev[[i]]$month
    data_cis <- addCIs(prev[[i]],prev[[i]]$positive,prev[[i]]$tested)%>%
      mutate(site=names(prev[i]),
             month = dates)
    df_data_all <- rbind(df_data_all,data_cis)
  }
  
  # df_data_all <- data.frame(t=numeric(),
  #                          tested=integer(),
  #                          positive=numeric(),
  #                          mean=numeric(),
  #                          upper=numeric(),
  #                          lower=numeric(),
  #                          district=character())
  
  # for(i in start:(start+number-1)){
  #   data_cis <- addCIs(prev_all[[i]],prev_all[[i]]$positive,prev_all[[i]]$tested)%>%
  #     mutate(district=names(prev_all[i]),
  #            month = dates_list[[country]])
  #   df_data_all <- rbind(df_data_all,data_cis)
  # }
  
  df_inc <- data.frame(time = integer(),
                       median = numeric(),
                       mean = numeric(),
                       upper = numeric(),
                       lower = numeric(),
                       site = character(),
                       month = character())
  df_inc_sample <- data.frame(time = integer(),
                              value = numeric(),
                              variable = character(),
                              site = character(),
                              month = character())
  for(i in 1:4){
    inc_history <- data.frame(t(results[[i]]$history['clininc_all', 51:1000, -1]))
    dates <- WA_pg_data_list[[i]]$month
    
    long_inc_sum <- inc_history%>%
      dplyr::mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value)*10000*30,
                       mean=mean(value)*10000*30,
                       upper=quantile(value,probs=0.975)*10000*30,
                       lower=quantile(value,probs=0.025)*10000*30)%>%
      mutate(site = districts[[i]],
             month = as.yearmon(dates))
    df_inc <- rbind(df_inc,long_inc_sum)
    
    inc_sample <- inc_history[, sample(ncol(inc_history), 100)] %>%
      mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value*10000*30,
                    site = districts[[i]],
                    month = as.yearmon(rep(dates,100)))
    df_inc_sample <- rbind(df_inc_sample,inc_sample)
  }
  
  df_eir <- data.frame(time = integer(),
                       median = numeric(),
                       mean = numeric(),
                       upper = numeric(),
                       lower = numeric(),
                       site = character(),
                       month = character(),
                       rel_eir = numeric())
  df_eir_sample <- data.frame(time = integer(),
                              value = numeric(),
                              variable = character(),
                              site = character(),
                              month = character(),
                              rel_eir = numeric())
  for(i in 1:4){
    eir_history <- data.frame(t(results[[i]]$history['EIR', 51:1000, -1]))
    dates <- WA_pg_data_list[[i]]$month
    
    long_eir_sum <- eir_history%>%
      dplyr::mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value),
                       mean=mean(value),
                       upper=quantile(value,probs=0.975),
                       lower=quantile(value,probs=0.025))%>%
      mutate(site = districts[[i]],
             month = as.yearmon(dates),
             rel_eir = median/max(median))
    df_eir <- rbind(df_eir,long_eir_sum)
    
    eir_sample <- eir_history[, sample(ncol(eir_history), 100)] %>%
      mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value,
                    site = districts[[i]],
                    month = as.yearmon(rep(dates,100)),
                    rel_eir = value/max(long_eir_sum$median))
    df_eir_sample <- rbind(df_eir_sample,eir_sample)
  }
  max_eir <- df_eir_sample %>%
    group_by(site) %>%
    summarise(max_eir = ifelse(max(value)==1000,800,max(value)*0.5))
  
  rainfall <- read_csv('./trial/Data/WA_ISTp_rainfall.csv') %>%
    group_by(Country) %>%
    filter(between(as.yearmon(Month), min(df_data_pg$month), max(df_data_pg$month)))%>%
    left_join(max_eir, by = join_by(Country == site))%>%
    mutate(rel_rainfall = Rainfall/max(Rainfall))%>%
    dplyr::rename(site = Country)
  
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  # est_inc_plot <- ggplot()+
  #   geom_line(data=df_inc_sample,aes(x=as.Date(month),y=value,color=site,group=variable),alpha=0.1)+
  #   # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  #   geom_line(data=df_inc,aes(x=as.Date(month),y=median,color=site,group=site),size=1)+
  #   # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
  #   scale_color_manual(values=colors)+
  #   scale_fill_manual(values=colors)+
  #   scale_x_date(date_labels = "%b %Y")+
  #   # coord_cartesian(ylim=c(0, 1000))+
  #   coord_cartesian(ylim=c(0, max(df_inc$median)*2))+
  #   # scale_y_continuous(limits=c(0,500))+
  #   # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  #   # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
  #   # labs(x='Date',y='EIR')+
  #   facet_grid(factor(site,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
  #   labs(title = 'Estimated Incidence\nper 10,000 person-mos.')+
  #   theme(legend.title = element_blank(),
  #         legend.position = 'none',
  #         axis.title.x = element_blank(),
  #         axis.title.y = element_blank(),
  #         axis.text.x=element_text(angle=45, hjust=1, vjust=1),
  #         axis.ticks.x = element_line(size = 0.5),
  #         axis.ticks.length = unit(3, "pt")
  #   )
  est_eir_plot <- ggplot()+
    # geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
    geom_line(data=df_eir_sample,aes(x=as.Date(month),y=value,color=site,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_eir,aes(x=as.Date(month),y=median,color=site,group=site),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # coord_cartesian(ylim=c(0, 1000))+
    # scale_y_continuous(limits=c(0,500))+
    scale_y_log10(breaks=c(.1,1,8,80,800),labels=c(.1,1,8,80,800))+
    coord_cartesian(ylim=c(.1, 800))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    facet_grid(factor(site,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.,scales = 'free_y')+
    labs(title = 'EIR - log scale\nInfectious bites/person/year')+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  rel_eir_plot <- ggplot()+
    geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
    # geom_line(data=df_eir_sample,aes(x=as.Date(month),y=rel_eir,color=site,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_eir,aes(x=as.Date(month),y=rel_eir,color=site,group=site),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # coord_cartesian(ylim=c(0, 1000))+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.1,1,10,100,1000),labels=c(.1,1,10,100,1000))+
    coord_cartesian(ylim=c(0, 1.2))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    facet_grid(factor(site,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    labs(title = 'Relative EIR and\nRelative rainfall')+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  # print('est_inc_plot')
  # print('obs_prev_plot_mg')
  # # obs_prev_plot_all <- annotations[[country]]+
  # #   # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  # #   # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
  # #   # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # #   geom_point(data=df_data_all,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
  # #   geom_errorbar(data=df_data_all,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
  # #   scale_color_manual(values=colors)+
  # #   scale_fill_manual(values=colors)+
  # #   scale_x_date(date_labels = "%b %Y")+
  # #   facet_grid(district~.)+
  # #   scale_y_continuous(limits = c(0,1))+
  # #   coord_cartesian(ylim=c(0, 1))+
  # #   labs(title = 'ANC Prevalence\nAll gravidities')+
  # #   theme(legend.title = element_blank(),
  # #         axis.title.x = element_blank(),
  # #         axis.title.y = element_blank(),
  # #         axis.text.x=element_text(angle=45, hjust=1, vjust=1),
  # #         axis.ticks.x = element_line(linewidth = 0.5),
  # #         axis.ticks.length = unit(3, "pt"),
  # #         legend.position = 'none'
  # #   )
  est_prev_plot <- ggplot()+
    geom_line(data=df_prev_sample,aes(x=as.Date(month),y=prev_pg,color=site,group=variable),alpha=0.2)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=site,group=site),alpha=0.2)+
    geom_line(data=df_prev,aes(x=as.Date(month),y=median,color=site,group=site),linewidth=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data_all,aes(x=as.Date(month),y=mean,color=site,group=site),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_all,aes(x=as.Date(month),ymin=lower,ymax=upper,color=site,group=site),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits=c(0,1))+
    coord_cartesian(xlim = range(as.Date(df_data_all$month)))+
    facet_grid(factor(site,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    labs(title = 'ANC Prevalence\nPrimi- and secundigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  # print('est_prev_plot')
  sample_size_plot <- ggplot()+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=site,group=site),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=site,group=site),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_line(data=df_data_all,aes(x=as.Date(month),y=tested,color=site,group=site))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    facet_grid(factor(site,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    labs(title = 'Number Tested')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  sample_size_plot+est_prev_plot+est_eir_plot+rel_eir_plot+ plot_layout(guides = "collect",ncol=4)
}


site <- 'Burkina Faso'

create_dashboard_plots_single <- function(results,prev_pg=WA_pg_data_list,prev_mg=WA_sg_data_list,site,burnin=50,chain_length=1000,sample=100){
  ## fn to return prevalence from log_odds
  get_prev_from_log_odds<-function(log_odds){
    return(exp(log_odds)/(1+exp(log_odds)))
  }
  
  ## fn to return odds from prevalence
  get_odds_from_prev<-function(prev){
    return(prev/(1-prev))
  }
  
  site_code <- c(`Burkina Faso` = 1, Gambia = 2, Ghana = 3, Mali = 4)
  i <- site_code[[site]]
  colors <- c(`Burkina Faso` = "#1B9E77", Gambia = "#999999", Ghana = "#D95F02", Mali = "#377EB8")

  prev_history <- data.frame(t(results$history['prev', (burnin+1):chain_length, -1]))
  dates <- prev_pg[[i]]$month
  
  sample_index <- sample(ncol(prev_history),sample)

  if(!is.null(results$seas_history)){
    
    prehistory_seas <- data.frame(bind_rows(results$seas_history[sample_index+burnin],.id='id'))
    prehistory_seas$date <- as.Date(min(as.Date(dates))-(max(prehistory_seas$t)-prehistory_seas$t))

  }
  coefs_pg_df <- as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS'))
  
  df_prev <- prev_history%>%
    mutate(t=c(1:nrow(prev_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(logodds_child = log(get_odds_from_prev(value)),
           prev_pg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
             prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
           }))))%>%
    group_by(time)%>%
    summarise(median=median(prev_pg),
              mean=mean(prev_pg),
              upper=quantile(prev_pg,probs=0.975),
              lower=quantile(prev_pg,probs=0.025))%>%
    mutate(site = site,
           month = as.yearmon(dates))
  # df_prev <- prev_history%>%
  #   mutate(t=c(1:nrow(prev_history)))%>%
  #   melt(id='t')%>%
  #   rename(time=t)%>%
  #   group_by(time)%>%
  #   summarise(median=median(value),
  #             mean=mean(value),
  #             upper=quantile(value,probs=0.975),
  #             lower=quantile(value,probs=0.025))%>%
  #   mutate(site = site,
  #          month = as.yearmon(dates))
  
  df_prev_sample <- prev_history %>%
    select(all_of(sample_index))%>%
    mutate(t=c(1:nrow(prev_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    mutate(logodds_child = log(get_odds_from_prev(value)),
           prev_pg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
             prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
           }))))%>%
    mutate(site = site,
           month = as.yearmon(rep(dates,sample)))
  # df_prev_sample <- prev_history %>%
  #   select(all_of(sample_index))%>%
  #   mutate(t=c(1:nrow(prev_history)))%>%
  #   melt(id='t')%>%
  #   rename(time=t)%>%
  #   mutate(site = site,
  #          month = as.yearmon(rep(dates,sample)))
  # 
  dates_pg <- prev_pg[[i]]$month
  df_data_pg <- addCIs(prev_pg[[i]],prev_pg[[i]]$positive,prev_pg[[i]]$tested)%>%
    mutate(site=site,
           month = dates_pg)
  
  
  dates_mg <- prev_mg[[i]]$month
  df_data_mg <- addCIs(prev_mg[[i]],prev_mg[[i]]$positive,prev_mg[[i]]$tested)%>%
    mutate(site=site,
           month = dates_mg)
  
  
  inc_history <- data.frame(t(results$history['clininc_all', (burnin+1):chain_length, -1]))
  
  df_inc <- inc_history%>%
    dplyr::mutate(t=c(1:nrow(inc_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    group_by(time)%>%
    dplyr::summarise(median=median(value)*10000*30,
                     mean=mean(value)*10000*30,
                     upper=quantile(value,probs=0.975)*10000*30,
                     lower=quantile(value,probs=0.025)*10000*30)%>%
    mutate(site = site,
           month = as.yearmon(dates))
  
  df_inc_sample <- inc_history %>%
    select(all_of(sample_index))%>%
    mutate(t=c(1:nrow(inc_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    dplyr::mutate(value = value*10000*30,
                  site = site,
                  month = as.yearmon(rep(dates,sample)))
  
  
  eir_history <- data.frame(t(results$history['EIR', (burnin+1):chain_length, -1]))
  
  df_eir <- eir_history%>%
    dplyr::mutate(t=c(1:nrow(eir_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    group_by(time)%>%
    dplyr::summarise(median=median(value),
                     mean=mean(value),
                     upper=quantile(value,probs=0.975),
                     lower=quantile(value,probs=0.025))%>%
    mutate(site = site,
           month = as.yearmon(dates))
  
  df_eir_sample <- eir_history %>%
    select(all_of(sample_index))%>%
    mutate(t=c(1:nrow(eir_history)))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    dplyr::mutate(value = value,
                  site = site,
                  month = as.yearmon(rep(dates,sample)))
  
  if(!is.null(results$seas_history)){
    plot_base_prev <- ggplot()+
      geom_line(data=prehistory_seas,aes(x=date,y=prev,group=id,color=site),alpha=0.1)
    plot_base_inc <- ggplot()+
      geom_line(data=prehistory_seas,aes(x=date,y=inc*10000*30,group=id,color=site),alpha=0.1)
    plot_base_eir <- ggplot()+
      geom_line(data=prehistory_seas,aes(x=date,y=EIR.1,group=id,color=site),alpha=0.1)
  }else{
    plot_base_prev <- ggplot()
    plot_base_inc <- ggplot()
    plot_base_eir <- ggplot()
  }
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  est_inc_plot <- plot_base_inc+
    geom_line(data=df_inc_sample,aes(x=as.Date(month),y=value,color=site,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_inc,aes(x=as.Date(month),y=median,color=site,group=site),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # coord_cartesian(ylim=c(0, 1000))+
    # coord_cartesian(ylim=c(0, max(df_inc$median)*2))+
    coord_cartesian(xlim = c(min(as.Date(dates))-2*365,NA), ylim=c(0, max(df_inc$median)*2))+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    labs(title = 'Estimated Incidence\nper 10,000 person-months')+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  est_eir_plot <- plot_base_eir+
    geom_line(data=df_eir_sample,aes(x=as.Date(month),y=value,color=site,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_eir,aes(x=as.Date(month),y=median,color=site,group=site),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # coord_cartesian(ylim=c(0, 1000))+
    coord_cartesian(xlim = c(min(as.Date(dates))-2*365,NA),ylim=c(0, 1500))+
    # coord_cartesian(ylim=c(0, 1500))+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    labs(title = 'EIR\n(infectious bites per person per year)')+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  # print('est_inc_plot')
  obs_prev_plot_mg <- ggplot()+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=site,group=site),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=site,group=site),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data_mg,aes(x=as.Date(month),y=mean,color=site,group=site),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_mg,aes(x=as.Date(month),ymin=lower,ymax=upper,color=site,group=site),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits = c(0,1))+
    coord_cartesian(ylim=c(0, 1),xlim = c(min(as.Date(dates))-2*365,NA))+
    # coord_cartesian(ylim=c(0, 1))+
    labs(title = 'ANC Prevalence\nSecundigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  # print('obs_prev_plot_mg')
  # # obs_prev_plot_all <- annotations[[country]]+
  # #   # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  # #   # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
  # #   # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # #   geom_point(data=df_data_all,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
  # #   geom_errorbar(data=df_data_all,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
  # #   scale_color_manual(values=colors)+
  # #   scale_fill_manual(values=colors)+
  # #   scale_x_date(date_labels = "%b %Y")+
  # #   facet_grid(district~.)+
  # #   scale_y_continuous(limits = c(0,1))+
  # #   coord_cartesian(ylim=c(0, 1))+
  # #   labs(title = 'ANC Prevalence\nAll gravidities')+
  # #   theme(legend.title = element_blank(),
  # #         axis.title.x = element_blank(),
  # #         axis.title.y = element_blank(),
  # #         axis.text.x=element_text(angle=45, hjust=1, vjust=1),
  # #         axis.ticks.x = element_line(linewidth = 0.5),
  # #         axis.ticks.length = unit(3, "pt"),
  # #         legend.position = 'none'
  # #   )
  est_prev_plot <- plot_base_prev+
    # geom_line(data=df_prev_sample,aes(x=as.Date(month),y=value,color=site,group=variable),alpha=0.2)+
    geom_line(data=df_prev_sample,aes(x=as.Date(month),y=prev_pg,color=site,group=variable),alpha=0.2)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=site,group=site),alpha=0.2)+
    geom_line(data=df_prev,aes(x=as.Date(month),y=median,color=site,group=site),linewidth=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data_pg,aes(x=as.Date(month),y=mean,color=site,group=site),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_pg,aes(x=as.Date(month),ymin=lower,ymax=upper,color=site,group=site),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits=c(0,1))+
    coord_cartesian(xlim = c(min(as.Date(dates))-2*365,NA))+
    labs(title = 'ANC Prevalence\nPrimigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  # print('est_prev_plot')
  sample_size_plot <- ggplot()+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=site,group=site),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=site,group=site),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_line(data=df_data_pg,aes(x=as.Date(month),y=tested,color=site,group=site))+
    geom_line(data=df_data_mg,aes(x=as.Date(month),y=tested,color=site,group=site),linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    facet_grid(site~.)+
    labs(title = 'Number Tested - Primigrav (solid)\nSecundigrav (dashed)')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  all <- sample_size_plot+est_prev_plot+obs_prev_plot_mg+est_eir_plot+est_inc_plot+ plot_layout(guides = "collect")
  estimates <- est_prev_plot+est_eir_plot+est_inc_plot+ plot_layout(guides = "collect")
  return(list(all=all,estimates=estimates))
}

ggsave('./trial/figures/primi_fitting_dash_060623.pdf',width=15,height=10)
  # obs_prev_plot_mg

windows(30,15)
wa_seas_2_dash <- create_dashboard_plots(results=wa_pg_seas_2_result_list)
wa_seas_2_dash

wa_std_dash <- create_dashboard_plots(results=wa_pg_std_result_list)
windows(30,15)
wa_std_dash

wa_seas_padded_dash <- create_dashboard_plots(results=wa_pg_seas_padded_result_list)
windows(30,15)
wa_seas_padded_dash

wa_seas_padded_dash_new_gamb
wa_pg_seas_padded_gamb_result_list <- wa_pg_seas_padded_result_list
wa_pg_seas_padded_gamb_result_list[[2]] <- wa_pg_gamb_seas$result()
wa_seas_padded_gamb_dash <- create_dashboard_plots(results=wa_pg_seas_padded_gamb_result_list)
windows(30,15)
wa_seas_padded_gamb_dash

wa_pg_bfall_seas_pc_dash <- create_dashboard_plots_single(results=wa_pg_bfall_seas_pc,burnin = 0,chain_length = 20,sample=20,site = 'Burkina Faso')
windows(30,15)
wa_pg_bfall_seas_pc_dash

wa_pg_bfall_seas_16_dash <- create_dashboard_plots_single(results=wa_pg_bfall_seas_16$result(),prev_pg = list(data_bf_all,data_bf_all),burnin = 50,chain_length = 250,sample=100,site = 'Burkina Faso')
windows(30,15)
wa_pg_bfall_seas_16_dash$estimates
results <- wa_pg_bfall_seas_16$result()
prehistory <- results$seas_history
results$history

wa_pg_bf_seas_16_u5prev_dash <- create_dashboard_plots_single(results=wa_pg_bf_seas_16_u5prev$result(),burnin = 50,chain_length = 500,sample=100,site = 'Burkina Faso')
windows(30,15)
wa_pg_bf_seas_16_u5prev_dash$estimates

wa_pg_bfall_seas_16_u5prev_dash <- create_dashboard_plots_single(results=wa_pg_bfall_seas_16_u5prev$result(),prev_pg = list(data_bf_all,data_bf_all),burnin = 50,chain_length = 500,sample=100,site = 'Burkina Faso')
wa_pg_bfall_seas_16_u5prev_dash$estimates
wa_pg_bf_seas_16_u5prev_dash$estimates
wa_pg_bfall_seas_16_sampsz_dash <- create_dashboard_plots_single(results=wa_pg_bfall_seas_16_sampsz$result(),prev_pg = list(data_bf_all_standardized_cis,data_bf_all_standardized_cis),burnin = 50,chain_length = 250,sample=100,site = 'Burkina Faso')
windows(30,15)
wa_pg_bfall_seas_16_sampsz_dash$estimates

wa_pg_bf_seas_16_dash <- create_dashboard_plots_single(results=wa_pg_bf_seas_16$result(),burnin = 50,chain_length = 500,sample=100,site = 'Burkina Faso')
windows(30,15)
wa_pg_bf_seas_16_dash

wa_pg_bulk_seas_090623_dash <- create_dashboard_plots(results=wa_pg_bulk_seas_090623_result_list)
windows(30,15)
wa_pg_bulk_seas_090623_dash

wa_all_bulk_seas_090623_dash <- create_dashboard_plots_all(results=wa_all_bulk_seas_090623_result_list, prev=WA_all_data_list)
windows(30,15)
wa_all_bulk_seas_090623_dash

wa_pgsg_bulk_sifter_betaa_dash <- create_dashboard_plots(results=wa_pgsg_bulk_sifter_betaa_results)

##Sifter
source('./trial/create_plots_trial.R')
wa_pgsg_bulk_sifter_eir_results
wa_pgsg_bulk_sifter_eir_dash <- create_dashboard_plots_trial(wa_pgsg_bulk_sifter_eir_results,
                                                               prev_pg=WA_pg_data_list,
                                                               prev_mg=WA_sg_data_list,
                                                               coefs_pg_df = as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS')),
                                                               coefs_mg_df = as.data.frame(readRDS('./nnp/Corr/sg_corr_sample.RDS')),
                                                               start_pf_time = 30*4,
                                                               param = 'EIR')
wa_pgsg_bulk_sifter_eir_dash$rel_dash
ggsave(paste0('Q:/anc_pmcmc/trial/figures/wa_pgsg_bulk_sifter_eir_dash_rel.pdf'), plot = wa_pgsg_bulk_sifter_eir_dash$rel_dash, width = 12, height = 6)

wa_pgsg_bulk_sifter_betaa_dash <- create_dashboard_plots_trial(wa_pgsg_bulk_sifter_betaa_results,
                                                               prev_pg=WA_pg_data_list,
                                                               prev_mg=WA_sg_data_list,
                                                               coefs_pg_df = as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS')),
                                                               coefs_mg_df = as.data.frame(readRDS('./nnp/Corr/sg_corr_sample.RDS')),
                                                               start_pf_time = 30*4,
                                                               param = 'betaa')
windows(30,15)
wa_pgsg_bulk_sifter_betaa_dash$rel_dash
ggsave(paste0('Q:/anc_pmcmc/trial/figures/wa_pgsg_bulk_sifter_betaa_dash_rel_max62.pdf'), plot = wa_pgsg_bulk_sifter_betaa_dash$rel_dash, width = 12, height = 6)

wa_pgsg_bulk_sifter_betaa_max125_dash <- create_dashboard_plots_trial(wa_pgsg_bulk_sifter_betaa_max125_results,
                                                               prev_pg=WA_pg_data_list,
                                                               prev_mg=WA_sg_data_list,
                                                               coefs_pg_df = as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS')),
                                                               coefs_mg_df = as.data.frame(readRDS('./nnp/Corr/sg_corr_sample.RDS')),
                                                               start_pf_time = 30*4,
                                                               param = 'betaa')
windows(30,15)
wa_pgsg_bulk_sifter_betaa_max125_dash$rel_dash
ggsave(paste0('Q:/anc_pmcmc/trial/figures/wa_pgsg_bulk_sifter_betaa_max125_rel.pdf'), plot = wa_pgsg_bulk_sifter_betaa_max125_dash$rel_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/trial/figures/wa_pgsg_bulk_sifter_betaa_max125_full.pdf'), plot = wa_pgsg_bulk_sifter_betaa_max125_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/trial/figures/wa_pgsg_bulk_sifter_betaa_max125_mcmc.pdf'), plot = wa_pgsg_bulk_sifter_betaa_max125_dash$mcmc_dash, width = 6, height = 6)

wa_pgsg_bulk_sifter_betaa_2407_dash <- create_dashboard_plots_trial(wa_pgsg_bulk_sifter_betaa_2407_results,
                                                                      prev_pg=WA_pg_data_list,
                                                                      prev_mg=WA_sg_data_list,
                                                                      coefs_pg_df = as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS')),
                                                                      coefs_mg_df = as.data.frame(readRDS('./nnp/Corr/sg_corr_sample.RDS')),
                                                                      start_pf_time = 30*4,
                                                                      param = 'betaa')
windows(30,15)
wa_pgsg_bulk_sifter_betaa_2407_dash$rel_dash
wa_pgsg_bulk_sifter_betaa_2407_dash$overlay_plot

ggsave(paste0('Q:/anc_pmcmc/trial/figures/wa_pgsg_bulk_sifter_betaa_2407_rel.pdf'), plot = wa_pgsg_bulk_sifter_betaa_2407_dash$rel_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/trial/figures/wa_pgsg_bulk_sifter_betaa_2407_full.pdf'), plot = wa_pgsg_bulk_sifter_betaa_2407_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/trial/figures/wa_pgsg_bulk_sifter_betaa_2407_mcmc.pdf'), plot = wa_pgsg_bulk_sifter_betaa_2407_dash$mcmc_dash, width = 6, height = 6)
ggsave(paste0('Q:/anc_pmcmc/trial/figures/wa_pgsg_bulk_sifter_betaa_2407_overlay.pdf'), plot = wa_pgsg_bulk_sifter_betaa_2407_dash$overlay_plot, width = 6, height = 12)

source('./trial/create_plots_trial_ea.R')
ea_pgmg_bulk_sifter_eir_max2000_dash <- create_dashboard_plots_trial_ea(ea_pgmg_bulk_sifter_eir_max2000_results,
                                                                      prev_pg=EA_pg_data_list,
                                                                      prev_mg=EA_mg_data_list,
                                                                      coefs_pg_df = as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS')),
                                                                      coefs_mg_df = as.data.frame(readRDS('./nnp/Corr/mg_corr_sample.RDS')),
                                                                      start_pf_time = 30*4,
                                                                      param = 'EIR')
windows(30,15)
ea_pgmg_bulk_sifter_eir_max2000_dash$rel_dash
ggsave(paste0('Q:/anc_pmcmc/trial/figures/ea_pgmg_bulk_sifter_eir_max2000.pdf'), plot = ea_pgmg_bulk_sifter_eir_max2000_dash$rel_dash, width = 12, height = 6)

ea_pgmg_bulk_sifter_betaa_max125_dash <- create_dashboard_plots_trial_ea(ea_pgmg_bulk_sifter_betaa_max125_results,
                                                                        prev_pg=EA_pg_data_list,
                                                                        prev_mg=EA_mg_data_list,
                                                                        coefs_pg_df = as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS')),
                                                                        coefs_mg_df = as.data.frame(readRDS('./nnp/Corr/mg_corr_sample.RDS')),
                                                                        start_pf_time = 30*4,
                                                                        param = 'betaa')
windows(30,15)
ea_pgmg_bulk_sifter_betaa_max125_dash$rel_dash
ggsave(paste0('Q:/anc_pmcmc/trial/figures/ea_pgmg_bulk_sifter_betaa_max125.pdf'), plot = ea_pgmg_bulk_sifter_betaa_max125_dash$rel_dash, width = 12, height = 6)
