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
source('shared/addCIs_inc.R')
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
result <- nnp_pgmg_bulk_sifter_betaa_max125_results[[1]]
country <- 'Burkina Faso'
district <- 'Banfora'
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

diag_std_mzbf_plots <- lapply(1:6,function(i) create_diag_figs(std_mzbf_result_list[[i]],country = country[i],district = district_list[i]))
diag_seas_mzbf_plots <- lapply(1:6,function(i) create_diag_figs(seas_mzbf_result_list[[i]],country = country[i],district = district_list[i]))
diag_std_pgmg_plots <- lapply(1:10,function(i) create_diag_figs(nnp_mgcorr_bulk_std_results[[i]],country = country[i],district = district_list[i]))
diag_seas_pg_plots <- lapply(1:10,function(i) create_diag_figs(nnp_pgcorr_bulk_seas_results[[i]],country = country[i],district = district_list[i]))
diag_std_pg_plots <- lapply(1:10,function(i) create_diag_figs(nnp_pgcorr_bulk_std_results[[i]],country = country[i],district = district_list[i]))
diag_seas_pgmg_plots <- lapply(1:10,function(i) create_diag_figs(nnp_mgcorr_bulk_seas_results_update[[i]],country = country[i],district = district_list[i]))

diag_pgmg_sifter_eir_plots <- lapply(1:10,function(i) create_diag_figs(nnp_pgmg_bulk_sifter_eir_results[[i]],country = country[i],district = district_list[i]))

diag_pgmg_sifter_betaa_plots<- lapply(1:10,function(i) create_diag_figs(nnp_pgmg_bulk_sifter_betaa_results[[i]],country = country[i],district = district_list[i]))

diag_pgmg_sifter_betaa_max125_plots<- lapply(1:10,function(i) create_diag_figs(nnp_pgmg_bulk_sifter_betaa_max125_results[[i]],country = country[i],district = district_list[i]))

diag_pgmg_sifter_betaa_2407_plots<- lapply(1:10,function(i) create_diag_figs(nnp_pgmg_bulk_sifter_betaa_2407_results[[i]],country = country[i],district = district_list[i]))

diag_pgmg_sifter_betaa_plots[1]
diag_pgmg_sifter_betaa_plots[2]

create_diag_figs(nnp_pgmg_bulk_sifter_betaa_2407.2$tasks[[1]]$result(),country=country[2],district=district_list[2])

country <- c('Burkina Faso','Burkina Faso','Burkina Faso',
             'Mozambique','Mozambique','Mozambique',
             'Nigeria','Nigeria','Nigeria','Nigeria')
district_list <- c('Banfora','Gaoua','Orodara','Changara','Chemba','Guro','Asa','Ejigbo','Ife North','Moro')

ggsave("test.pdf", width = 10, height = 7)
for(i in c(1:6)){
  windows(10,7)
  print(diag_std_mzbf_plots[i])
  ggsave(paste0('Q:/anc_pmcmc/nnp/figures/diag/',district_list[i],'_std_310223.pdf'), width = 10, height = 7)
}
for(i in c(1:6)){
  windows(10,7)
  print(diag_seas_mzbf_plots[i])
  ggsave(paste0('Q:/anc_pmcmc/nnp/figures/diag/',district_list[i],'_seas_310223.pdf'), width = 10, height = 7)
}
for(i in c(1:6)){
  windows(10,7)
  print(diag_seas_mzbf_plots[i])
  ggsave(paste0('Q:/anc_pmcmc/nnp/figures/diag/',district_list[i],'_seas_310223.pdf'), width = 10, height = 7)
}
for(i in c(1:10)){
  windows(10,7)
  print(diag_std_pgmg_plots[i])
  ggsave(paste0('Q:/anc_pmcmc/nnp/figures/diag/',district_list[i],'_pgmg_std_220223.pdf'), width = 10, height = 7)
}
for(i in c(1:10)){
  windows(10,7)
  print(diag_seas_pg_plots[i])
  ggsave(paste0('Q:/anc_pmcmc/nnp/figures/diag/',district_list[i],'_pg_seas_270223.pdf'), width = 10, height = 7)
}
for(i in c(1:10)){
  windows(10,7)
  print(diag_std_pg_plots[i])
  ggsave(paste0('Q:/anc_pmcmc/nnp/figures/diag/',district_list[i],'_pg_std_280223.pdf'), width = 10, height = 7)
}
for(i in c(1:10)){
  windows(10,7)
  print(diag_seas_pgmg_plots[i])
  ggsave(paste0('Q:/anc_pmcmc/nnp/figures/diag/',district_list[i],'_pgmg_seas_020323.pdf'), width = 10, height = 7)
}
for(i in c(1:10)){
  ggsave(paste0('Q:/anc_pmcmc/nnp/figures/diag/',district_list[i],'_pgmg_sifter_eir.pdf'),plot=diag_pgmg_sifter_eir_plots[[i]], width = 10, height = 7)
}
for(i in c(1:10)){
  ggsave(paste0('Q:/anc_pmcmc/nnp/figures/diag/',district_list[i],'_pgmg_sifter_betaa.pdf'),plot=diag_pgmg_sifter_betaa_plots[[i]], width = 10, height = 7)
}
for(i in c(1:10)){
  ggsave(paste0('Q:/anc_pmcmc/nnp/figures/diag/',district_list[i],'_pgmg_sifter_betaa.pdf'),plot=diag_pgmg_sifter_betaa_plots[[i]], width = 10, height = 7)
}
for(i in c(1:10)){
  ggsave(paste0('Q:/anc_pmcmc/nnp/figures/diag/',district_list[i],'_pgmg_max125_betaa.pdf'),plot=diag_pgmg_sifter_betaa_max125_plots[[i]], width = 7, height = 7)
}
for(i in c(1:10)){
  ggsave(paste0('Q:/anc_pmcmc/nnp/figures/diag/',district_list[i],'_pgmg_2407_betaa.pdf'),plot=diag_pgmg_sifter_betaa_2407_plots[[i]], width = 7, height = 7)
}
diag_pgmg_sifter_betaa_max125_plots[[1]] + diag_pgmg_sifter_betaa_2407_plots[[1]] + plot_layout(ncol=2)
nnp_mgcorr_bulk_seas_results_update
summary_init_EIR <- bind_rows(lapply(1:10,function(x,result){
  exp(quantile(result[[x]]$mcmc[501:1000,'log_init_EIR'], probs = c(0.025,0.5,0.975)))
},result=nnp_mgcorr_bulk_seas_results_update))
names(nnp_pg_list) <- district_list
names(nnp_mg_list) <- district_list

nnp_pg_list <- list(Banfora = data_raw_bf_pg_banfora,Gaoua = data_raw_bf_pg_gaoua,Orodara = data_raw_bf_pg_orodara,
                 Changara = data_raw_mz_pg_changara, Chemba = data_raw_mz_pg_chemba,Guro = data_raw_mz_pg_guro,
                 Asa = data_raw_ng_pg_asa,Ejigbo = data_raw_ng_pg_ejigbo,`Ife North` = data_raw_ng_pg_ifenorth, Moro = data_raw_ng_pg_moro)
nnp_mg_list <- list(Banfora = data_raw_bf_mg_banfora,Gaoua = data_raw_bf_mg_gaoua,Orodara = data_raw_bf_mg_orodara,
                    Changara = data_raw_mz_mg_changara, Chemba = data_raw_mz_mg_chemba,Guro = data_raw_mz_mg_guro,
                    Asa = data_raw_ng_mg_asa,Ejigbo = data_raw_ng_mg_ejigbo,`Ife North` = data_raw_ng_mg_ifenorth, Moro = data_raw_ng_mg_moro)
nnp_all_list <- list(Banfora = BF_all_banfora,Gaoua = BF_all_gaoua,Orodara = BF_all_orodara,
                    Changara = MZ_all_changara, Chemba = MZ_all_chemba,Guro = MZ_all_guro,
                    Asa = NG_all_asa,Ejigbo = NG_all_ejigbo,`Ife North` = NG_all_ifenorth, Moro = NG_all_moro)

create_prev_plots <- function(results,data_list=nnp_list,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                 MZ = c('Changara','Chemba','Guro'),
                 NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-6-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2022-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2022-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   distrit = character())
  
  for(i in start:(start+number-1)){
    prev_history <- data.frame(t(results[[i]]$history['prev', 51:1000, -1]))
    long_prev_sum <- prev_history%>% mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      group_by(time)%>%
      summarise(median=median(value),
                mean=mean(value),
                upper=quantile(value,probs=0.975),
                lower=quantile(value,probs=0.025))%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df <- rbind(df,long_prev_sum)
    
  }
  df_data <- data.frame(t=numeric(),
                        tested=integer(),
                        positive=numeric(),
                        mean=numeric(),
                        upper=numeric(),
                        lower=numeric(),
                        district=character())
  
  for(i in start:(start+number-1)){
    data_cis <- addCIs(data_list[[i]],data_list[[i]]$positive,data_list[[i]]$tested)%>%
      mutate(district=names(data_list[i]),
             month = dates_list[[country]])
    df_data <- rbind(df_data,data_cis)
  }
  
  annotations <- list(
    BF = ggplot(df)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-31'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-31'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    MZ = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-1-1'), xmax = as.Date('2022-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    NG = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-7-1'), xmax = as.Date('2022-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")
  )
  prev_plot <- annotations[[country]]+
    geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    labs(x='Date',y='ANC Prevalence')+
    theme(legend.title = element_blank(),
          # axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"),
          )
  prev_plot
}

##Prev plots
BF_prev_pg <- create_prev_plots(results = std_mzbf_result_list,data_list=nnp_pg_list,country='BF')
MZ_prev_pg <- create_prev_plots(results = std_mzbf_result_list,data_list=nnp_pg_list,country='MZ')
NG_prev_pg <- create_prev_plots(results = nnp_pg_result_list,data_list=nnp_pg_list,country='NG')

windows(10,7)
BF_prev_pg
MZ_prev_pg
NG_prev_pg


##Inc and EIR plots
create_inc_plots <- function(results,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-6-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2022-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2022-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   distrit = character())
  
  for(i in start:(start+number-1)){
    inc_history <- data.frame(t(results[[i]]$history['inc', 51:1000, -1]))
    long_inc_sum <- inc_history%>% mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      group_by(time)%>%
      summarise(median=median(value)*10000*30,
                mean=mean(value)*10000*30,
                upper=quantile(value,probs=0.975)*10000*30,
                lower=quantile(value,probs=0.025)*10000*30)%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df <- rbind(df,long_inc_sum)
    
  }

  df_eir <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   distrit = character())
  
  for(i in start:(start+number-1)){
    eir_history <- data.frame(t(results[[i]]$history['EIR', 51:1000, -1]))
    long_eir_sum <- eir_history%>% mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      group_by(time)%>%
      summarise(median=median(value),
                mean=mean(value),
                upper=quantile(value,probs=0.975),
                lower=quantile(value,probs=0.025))%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df_eir <- rbind(df_eir,long_eir_sum)
  }
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  annotations <- list(
    BF = ggplot(df)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits = c(0,3750)),
    MZ = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-1-1'), xmax = as.Date('2022-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    scale_y_continuous(limits=c(0,400)),
    NG = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-7-1'), xmax = as.Date('2022-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,3750))

  )
  inc_plot <- annotations[[country]]+
    geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    labs(x='Date',y='Estimated Incidence (<5 year old)\nper 10,000 person-months')+
    # labs(x='Date',y='EIR')+
    theme(legend.title = element_blank(),
          # axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"))
  # +
  #   coord_cartesian(ylim=c(0, 40))
  inc_plot
}

BF_inc_pg <- create_inc_plots(results = std_mzbf_result_list,country='BF')
MZ_inc_pg <- create_inc_plots(results = std_mzbf_result_list,country='MZ')
NG_inc_pg <- create_inc_plots(results = nnp_pg_bulk_std_all,country='NG')

BF_inc_pg <- create_inc_plots(results = seas_mzbf_result_list,country='BF')
MZ_inc_pg <- create_inc_plots(results = seas_mzbf_result_list,country='MZ')

windows(7,5)
BF_inc_pg
MZ_inc_pg
NG_inc_pg

create_relinc_plots <- function(results,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-5-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2021-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2022-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   distrit = character())
  
  for(i in start:(start+number-1)){
    inc_history <- data.frame(t(results[[i]]$history['inc', 51:1000, -1]))
    long_inc_sum <- inc_history%>% mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      group_by(time)%>%
      summarise(median=median(value)*10000,
                mean=mean(value)*10000,
                upper=quantile(value,probs=0.975)*10000,
                lower=quantile(value,probs=0.025)*10000)%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    long_relinc_sum <- long_inc_sum %>%
      mutate(median=median/max(long_inc_sum$median),
             lower=lower/max(long_inc_sum$upper),
             upper=upper/max(long_inc_sum$upper),
      )
    df <- rbind(df,long_relinc_sum)
    
  }
  
  annotations <- list(
    BF = ggplot(df)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-311'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-31'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    MZ = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    NG = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-7-1'), xmax = as.Date('2022-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")
  )
  inc_plot <- annotations[[country]]+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    labs(x='Date',y='Normalised Median Incidence\n(<5 year old)')+
    # labs(x='Date',y='EIR')+
    theme(legend.title = element_blank(),
          # axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"))
  # +
  #   coord_cartesian(ylim=c(0, 40))
  inc_plot
}
create_relinc_plots(nnp_result_list,country='BF')
create_relinc_plots(nnp_result_list,country='MZ')
create_relinc_plots(nnp_result_list,country='NG')

##EIR plots
create_eir_plots <- function(results,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-6-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2022-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2022-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df_eir <- data.frame(time = integer(),
                       median = numeric(),
                       mean = numeric(),
                       upper = numeric(),
                       lower = numeric(),
                       distrit = character())
  
  for(i in start:(start+number-1)){
    eir_history <- data.frame(t(results[[i]]$history['EIR', 51:1000, -1]))
    long_eir_sum <- eir_history%>% dplyr::mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value),
                mean=mean(value),
                upper=quantile(value,probs=0.975),
                lower=quantile(value,probs=0.025))%>%
      dplyr::mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df_eir <- rbind(df_eir,long_eir_sum)
  }
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  annotations <- list(
    BF = ggplot(df_eir)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-1'), ymin = 0, ymax = 1000,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-1'), ymin = 0, ymax = 1000,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits = c(0,105)),
    MZ = ggplot(df_eir)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-1-1'), xmax = as.Date('2022-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    scale_y_continuous(limits=c(0,400)),
    NG = ggplot(df_eir)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-1'), ymin = 0, ymax = 1000,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,45))
    
  )
  eir_plot <- annotations[[country]]+
    geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # scale_y_continuous(limits=c(0,500))+
    scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    labs(x='Date',y='EIR')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          # axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"))
  # +
  #   coord_cartesian(ylim=c(0, 40))
  eir_plot
}

BF_eir_pg <- create_eir_plots(results = std_mzbf_result_list,country='BF')
MZ_eir_pg <- create_eir_plots(results = std_mzbf_result_list,country='MZ')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/esteir_mz_std_310223.pdf'), width = 7, height = 5)
NG_eir_pg <- create_eir_plots(results = nnp_pg_result_list,country='NG')
MZ_eir_pg_seas <- create_eir_plots(results = seas_mzbf_result_list,country='MZ')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/esteir_mz_seas_310223.pdf'), width = 7, height = 5)

# ##Compare posterior distributions
# # Each country separately
# # Grab posterior distributions for each run type
# # Melt dataframe
# # Add variable, country name, and run type as new variable
# # rbind 4 different run types
# # Plot as violin plot (maybe ridge plot), variables as facets, run types as colors
# 
# create_comp_plots <- function(aa,u5,pg,pgcorr,country=c('BF','MZ','NG')){
#   district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
#                         MZ = c('Changara','Chemba','Guro'),
#                         NG = c('Asa','Ejigbo','Ife North','Moro'))
#   dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-5-1'),by='months'),
#                      MZ = seq(as.Date('2020-12-1'),as.Date('2021-9-1'),by='months'),
#                      NG = seq(as.Date('2020-11-1'),as.Date('2022-12-1'),by='months'))
#   colors <- c(`All ages` = "#1B9E77", `Under 5 years old` = "#999999", `Primigrav` = "#D95F02", `Primigrav Correlation` = "#377EB8")
#   
#   start_list <- c(BF = 1, MZ = 4, NG = 7)
#   number_list <- c(BF = 3, MZ = 3, NG = 4)
#   
#   districts <- district_list[[country]]
#   start <- start_list[[country]]
#   number <- number_list[[country]]
# 
#   df <- data.frame(value = numeric(),
#                    variable = character(),
#                    run = character(),
#                    distrit = character())
#   aa[[1]]$mcmc
# ##All ages results  
#   for(i in start:(start+number-1)){
#     mcmc <- data.frame(aa[[i]]$mcmc[51:1000,])
#     long_mcmc <- mcmc%>% mutate(step=c(1:nrow(mcmc)))%>%
#       melt(id='step')%>%
#       mutate(district = districts[[i-start+1]],
#              run = 'All ages')
#     df <- rbind(df,long_mcmc)
#     
#   }
#   ##Under 5 results  
#   for(i in start:(start+number-1)){
#     mcmc <- data.frame(u5[[i]]$mcmc[51:1000,])
#     long_mcmc <- mcmc%>% mutate(step=c(1:nrow(mcmc)))%>%
#       melt(id='step')%>%
#       mutate(district = districts[[i-start+1]],
#              run = 'Under 5 years old')
#     df <- rbind(df,long_mcmc)
#   }
#   ##Primigravs 
#   for(i in start:(start+number-1)){
#     mcmc <- data.frame(pg[[i]]$mcmc[51:1000,])
#     long_mcmc <- mcmc%>% mutate(step=c(1:nrow(mcmc)))%>%
#       melt(id='step')%>%
#       mutate(district = districts[[i-start+1]],
#              run = 'Primigrav')
#     df <- rbind(df,long_mcmc)
#   }
#   ##Primigravs 
#   for(i in start:(start+number-1)){
#     mcmc <- data.frame(pgcorr[[i]]$mcmc[51:1000,])
#     long_mcmc <- mcmc%>% mutate(step=c(1:nrow(mcmc)))%>%
#       melt(id='step')%>%
#       mutate(district = districts[[i-start+1]],
#              run = 'Primigrav Correlation')
#     df <- rbind(df,long_mcmc)
#   }
#   
#   df_out <- df%>%
#     mutate(run=factor(run,levels=c('All ages','Under 5 years old','Primigrav','Primigrav Correlation')),
#            value = ifelse(variable=='log_init_EIR',exp(value),value))
#   levels(df_out$variable)[levels(df_out$variable)=='log_init_EIR'] <- 'init_EIR'
# 
#   comp_prior_plot <- ggplot(data=df_out[df_out$variable=='log_prior',])+
#     geom_violin(aes(x=district,y=value,fill=run))+
#     scale_fill_manual(values=colors)+
#     labs(y='Log Prior')+
#     facet_wrap(~district,scales='free_x')+
#     theme(axis.title.x = element_blank())
#   comp_likelihood_plot <- ggplot(data=df_out[df_out$variable=='log_likelihood',])+
#     geom_violin(aes(x=district,y=value,fill=run))+
#     scale_fill_manual(values=colors)+
#     labs(y='Log Likelihood')+
#     facet_wrap(~district,scales='free_x')+
#     theme(axis.title.x = element_blank())
#   comp_posterior_plot <- ggplot(data=df_out[df_out$variable=='log_posterior',])+
#     geom_violin(aes(x=district,y=value,fill=run))+
#     scale_fill_manual(values=colors)+
#     labs(y='Log Posterior')+
#     facet_wrap(~district,scales='free_x')+
#     theme(axis.title.x = element_blank())
#   
#   probs_plot <- comp_prior_plot/comp_likelihood_plot/comp_posterior_plot +  plot_layout(guides = "collect")
#   
#   comp_eirsd_plot <- ggplot(data=df_out[df_out$variable=='EIR_SD',])+
#     geom_violin(aes(x=district,y=value,fill=run))+
#     scale_fill_manual(values=colors)+
#     labs(y='EIR Volatility')+
#     facet_wrap(~district,scales='free_x')+
#     theme(axis.title.x = element_blank())
#   comp_eir_plot <- ggplot(data=df_out[df_out$variable=='init_EIR',])+
#     geom_violin(aes(x=district,y=value,fill=run))+
#     scale_fill_manual(values=colors)+
#     labs(y='Initial EIR')+
#     facet_wrap(~district,scales='free_x')+
#     scale_y_log10()+
#     theme(axis.title.x = element_blank())
#   params_plot <- comp_eirsd_plot/comp_eir_plot +  plot_layout(guides = "collect")
# }
# #Burkina Faso
# windows(10,15)
# bf_comp_plot <- create_comp_plots(aa=nnp_result_list_aa,u5=nnp_result_list,
#                                   pg=nnp_pg_result_list,pgcorr=nnp_pgcorr_result_list,
#                                   country = 'BF')
##For presentation
create_obsprev_plots <- function(results,data_list=nnp_list,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-6-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2022-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2022-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   distrit = character())
  
  for(i in start:(start+number-1)){
    prev_history <- data.frame(t(results[[i]]$history['prev', 51:1000, -1]))
    long_prev_sum <- prev_history%>% mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      group_by(time)%>%
      summarise(median=median(value),
                mean=mean(value),
                upper=quantile(value,probs=0.975),
                lower=quantile(value,probs=0.025))%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df <- rbind(df,long_prev_sum)
    
  }
  df_data <- data.frame(t=numeric(),
                        tested=integer(),
                        positive=numeric(),
                        mean=numeric(),
                        upper=numeric(),
                        lower=numeric(),
                        district=character())
  
  for(i in start:(start+number-1)){
    data_cis <- addCIs(data_list[[i]],data_list[[i]]$positive,data_list[[i]]$tested)%>%
      mutate(district=names(data_list[i]),
             month = dates_list[[country]])
    df_data <- rbind(df_data,data_cis)
  }
  
  annotations <- list(
    BF = ggplot(df)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-31'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-31'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    MZ = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-1-1'), xmax = as.Date('2022-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    NG = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-7-1'), xmax = as.Date('2022-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")
    
  )
  prev_plot <- annotations[[country]]+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    labs(y='ANC Prevalence')+
    facet_grid(~district)+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5), 
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  prev_plot
}

BF_obsprev_pg <- create_obsprev_plots(results = std_mzbf_result_list,data_list=nnp_pg_list,country='BF')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/obsprev_bf_std_310223.pdf'), width = 7, height = 5)
MZ_obsprev_pg <- create_obsprev_plots(results = std_mzbf_result_list,data_list=nnp_pg_list,country='MZ')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/obsprev_mz_std_310223.pdf'), width = 7, height = 5)
NG_obsprev_pg <- create_obsprev_plots(results = std_mzbf_result_list,data_list=nnp_pg_list,country='NG')


windows(7,5)
BF_obsprev_pg
MZ_obsprev_pg
NG_obsprev_pg

create_estprev_plots <- function(results,data_list=nnp_list,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-6-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2022-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2022-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   district = character(),
                   month = character())
  df_sample <- data.frame(time = integer(),
                          value = numeric(),
                          variable = character(),
                          district = character(),
                          month = character())
  for(i in start:(start+number-1)){
    prev_history <- data.frame(t(results[[i]]$history['prev', 51:1000, -1]))
      
    
    long_prev_sum <- prev_history%>%
      mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      group_by(time)%>%
      summarise(median=median(value),
                mean=mean(value),
                upper=quantile(value,probs=0.975),
                lower=quantile(value,probs=0.025))%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df <- rbind(df,long_prev_sum)
    
    prev_sample <- prev_history[, sample(ncol(prev_history), 100)] %>%
      mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      mutate(district = districts[[i-start+1]],
             month = rep(dates_list[[country]],100))
    df_sample <- rbind(df_sample,prev_sample)
    
    
  }
  
  df_data <- data.frame(t=numeric(),
                        tested=integer(),
                        positive=numeric(),
                        mean=numeric(),
                        upper=numeric(),
                        lower=numeric(),
                        district=character())
  
  for(i in start:(start+number-1)){
    data_cis <- addCIs(data_list[[i]],data_list[[i]]$positive,data_list[[i]]$tested)%>%
      mutate(district=names(data_list[i]),
             month = dates_list[[country]])
    df_data <- rbind(df_data,data_cis)
  }
  
  annotations <- list(
    BF = ggplot(df)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-31'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-31'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    MZ = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-1-1'), xmax = as.Date('2022-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    NG = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-7-1'), xmax = as.Date('2022-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")
  )
  prev_plot <- annotations[[country]]+
    geom_line(data=df_sample,aes(x=month,y=value,color=district,group=variable),alpha=0.2)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(aes(x=month,y=median,color=district,group=district),linewidth=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    labs(y='ANC Prevalence')+
    facet_grid(~district)+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5), 
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  prev_plot
}
windows(7,5)

BF_estprev_pg_std <- create_estprev_plots(results = std_mzbf_result_list,data_list=nnp_pg_list,country='BF')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/estprev_bf_std_310223.pdf'), width = 7, height = 5)
MZ_estprev_pg_std <- create_estprev_plots(results = std_mzbf_result_list,data_list=nnp_pg_list,country='MZ')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/estprev_mz_std_310223.pdf'), width = 7, height = 5)
NG_estprev_pg <- create_estprev_plots(results = nnp_pg_result_list,data_list=nnp_pg_list,country='NG')
BF_estprev_pg
MZ_estprev_pg
BF_estprev_pg_seas <- create_estprev_plots(results = seas_mzbf_result_list,data_list=nnp_pg_list,country='BF')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/estprev_bf_seas_310223.pdf'), width = 7, height = 5)
MZ_estprev_pg_seas <- create_estprev_plots(results = seas_mzbf_result_list,data_list=nnp_pg_list,country='MZ')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/estprev_mz_seas_310223.pdf'), width = 7, height = 5)
NG_estprev_pg_seas <- create_estprev_plots(results = seas_all_result_list,data_list=nnp_pg_list,country='NG')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/estprev_ng_seas_310223.pdf'), width = 7, height = 5)

NG_estprev_pg / plot_spacer()
windows(6,7)
BF_estprev_pg/BF_inc_pg
MZ_estprev_pg/MZ_inc_pg
NG_estprev_pg/NG_inc_pg

##Create plots for BF and NG HMIS data

bf_hmis <- readxl::read_excel('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Burkina Faso/Routine HMIS/Routine Data All Ages.xlsx')
bf_hmis <- bf_hmis %>%
  rename(district = 'Distrist')%>%
  mutate(population_2017 = case_when(
    district=='Nouna' ~ 368395-(15001+54109),
    district=='Tougan' ~ 287801-(11720+42273),
    district=='Banfora' ~ 392498-(14131+50817),
    district=='Orodara' ~ 257258-(8262+31725),
    district=='Gaoua' ~ 252385-(9927+35702)
  ))
bf_hmis_u5 <- readxl::read_excel('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Burkina Faso/Routine HMIS/Routine Data Under 5.xlsx')
names(bf_hmis_u5) <- c('district','date','mal_severe','mal_simple','mal_simple_act','positive','tested')
bf_hmis_u5 <- bf_hmis_u5 %>%
  mutate(population = dplyr::case_when(
    district=='Nouna' ~ 15001+54109,
    district=='Tougan' ~ 11720+42273,
    district=='Banfora' ~ 14131+50817,
    district=='Orodara' ~ 8262+31725,
    district=='Gaoua' ~ 9927+35702),
    date = as.yearmon(date)
  )

bf_hmis <- readxl::read_excel('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Burkina Faso/Routine HMIS/Routine Data All Ages.xlsx')
bf_hmis <- bf_hmis %>%
  rename(district = 'Distrist')%>%
  mutate(population_2017 = case_when(
    district=='Nouna' ~ 368395-(15001+54109),
    district=='Tougan' ~ 287801-(11720+42273),
    district=='Banfora' ~ 392498-(14131+50817),
    district=='Orodara' ~ 257258-(8262+31725),
    district=='Gaoua' ~ 252385-(9927+35702)
  ))
mz_hmis <- readxl::read_excel('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Mozambique/Routine HMIS/Mozambique_Routine_Sept2022.xlsx')
ng_hmis <- readxl::read_excel('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl020323/Nigeria/Routine data/NNP Nigeria HMIS Data 2019-2022 - LGA rev 2023.03.02.xlsx')

bf_hmis_sum <- bf_hmis %>%
  group_by(Month,Year,`Net type`,district)%>%
  dplyr::summarise(confirmed = sum(Confirmed),
            population = sum(Population))%>%
  mutate(inc = confirmed*10000/population)

bf_hmis_sum$date = as.yearmon(paste(bf_hmis_sum$Year, bf_hmis_sum$Month), "%Y %b")
bf_hmis$date = as.yearmon(paste(bf_hmis$Year, bf_hmis$Month), "%Y %b")
bf_hmis_all <- left_join(bf_hmis,bf_hmis_u5,by=c('date','district'),suffix=c('.allpop','.u5'))
bf_hmis_all$positive_o5 <- bf_hmis_all$`RDT +`-bf_hmis_all$positive
bf_hmis_forplot <- bf_hmis_all[bf_hmis_all$district %in% c('Banfora','Gaoua','Orodara')&bf_hmis_all$date>=as.yearmon('Sep 2020')&
                                 bf_hmis_all$date<=as.yearmon('Jun 2022')&!is.na(bf_hmis_all$Confirmed),]
bf_hmis_forplot <- addCIs_inc(bf_hmis_forplot,bf_hmis_forplot$positive_o5,bf_hmis_forplot$population_2017)
bf_hmis_forplot <- bf_hmis_forplot %>%
  mutate(date_ex = as.Date(date, frac = 0),
         prop_pos = positive_o5/(Tested-tested))%>%
  dplyr::rename(inc = mean)

mz_hmis_filtered <- mz_hmis[mz_hmis$District %in% c('Changara','Chemba','Guro')&!is.na(mz_hmis$Positive),]%>%
  group_by(District, Month, Year, `Net type`)%>%
  dplyr::summarise(Tested=sum(Tested),
                   Positive=sum(Positive),
                   Population=sum(u5.Population))
mz_hmis_forplot <- addCIs_inc(mz_hmis_filtered,mz_hmis_filtered$Positive,mz_hmis_filtered$Population)%>%
  mutate(date_ex=as.Date(as.yearmon(paste(Year, Month), "%Y %b")),
         prop_pos = Positive/Tested)%>%
  filter(date_ex>=as.Date('2020-12-1')&date_ex<=as.Date('2022-09-1'))%>%
  dplyr::rename(district = District,
                inc = mean)

ng_hmis_filtered <- ng_hmis[!is.na(ng_hmis$year),]
ng_hmis_forplot <- addCIs_inc(ng_hmis_filtered,ng_hmis_filtered$rdt_positive,ng_hmis_filtered$pop_est_micro)%>%
  mutate(date=as.Date(period),
         prop_pos = rdt_positive/rdt_tested)%>%
  filter(date>=as.Date('2020-11-1')&date<=as.Date('2022-12-1'))%>%
  dplyr::rename(date_ex = date,
                district = lga,
                inc = mean)

colors_nets <- c(IG2 = "#1B9E77", Stdr = "#999999", PBO = "#D95F02")
windows(10,4)
ggplot(bf_hmis_sum)+
  geom_line(aes(x=date,y=inc,col=`Net type`),size=1)+
  scale_y_continuous(limits=c(0,900))+
  scale_color_manual(values=colors_nets)

create_inc_plots_threads <- function(results,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-6-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2022-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2022-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   district = character(),
                   month = character())
  df_sample <- data.frame(time = integer(),
                          value = numeric(),
                          variable = character(),
                          district = character(),
                          month = character())
  for(i in start:(start+number-1)){
    inc_history <- data.frame(t(results[[i]]$history['inc', 51:1000, -1]))
    
    long_inc_sum <- inc_history%>%
      dplyr::mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value)*10000*30,
                mean=mean(value)*10000*30,
                upper=quantile(value,probs=0.975)*10000*30,
                lower=quantile(value,probs=0.025)*10000*30)%>%
      dplyr::mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df <- rbind(df,long_inc_sum)
    
    inc_sample <- inc_history[, sample(ncol(inc_history), 100)] %>%
      mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value*10000*30,
             district = districts[[i-start+1]],
             month = rep(dates_list[[country]],100))
    df_sample <- rbind(df_sample,inc_sample)
  }
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  annotations <- list(
    BF = ggplot(df)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-1'), ymin = 0, ymax = 6000,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-1'), ymin = 0, ymax = 6000,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits = c(0,6000))+
      coord_cartesian(ylim=c(0, 4000)),
    MZ = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = 9000,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-1-1'), xmax = as.Date('2022-6-30'), ymin = 0, ymax = 9000,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,9000))+
      coord_cartesian(ylim=c(0, 4000)),
    NG = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-7-1'), xmax = as.Date('2022-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,3750))+
      coord_cartesian(ylim=c(0, 1000))
    
  )
  inc_plot <- annotations[[country]]+
    geom_line(data=df_sample,aes(x=month,y=value,color=district,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    labs(x='Date',y='Estimated Incidence\nper 10,000 person-months')+
    # labs(x='Date',y='EIR')+
    facet_grid(~district)+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"))
  # +
  #   coord_cartesian(ylim=c(0, 40))
  inc_plot
}
windows(7,5)

BF_inc_pg_threads <- create_inc_plots_threads(results = std_mzbf_result_list,country='BF')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/estinc_bf_std_310223.pdf'), width = 7, height = 5)
MZ_inc_pg_threads <- create_inc_plots_threads(results = std_mzbf_result_list,country='MZ')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/estinc_mz_std_310223.pdf'), width = 7, height = 5)
NG_inc_pg_threads <- create_inc_plots_threads(results = nnp_pg_result_list,country='NG')
MZ_inc_pg_threads_seas <- create_inc_plots_threads(results = seas_mzbf_result_list,country='MZ')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/estinc_mz_seas_310223.pdf'), width = 7, height = 5)

BF_inc_pg_threads_seas <- create_inc_plots_threads(results = seas_all_result_list,country='BF')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/estinc_bf_seas_310223.pdf'), width = 7, height = 5)
MZ_inc_pg_threads_seas <- create_inc_plots_threads(results = seas_all_result_list,country='MZ')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/estinc_mz_seas_310223.pdf'), width = 7, height = 5)
NG_inc_pg_threads_seas <- create_inc_plots_threads(results = seas_all_result_list,country='NG')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/estinc_ng_seas_310223.pdf'), width = 7, height = 5)

create_obsinc_plots <- function(results,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-6-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2022-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2022-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  annotations <- list(
    BF = ggplot(results)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")
      ,#scale_y_continuous(limits = c(0,3750)),
      #coord_cartesian(ylim=c(0, 2000)),
    MZ = ggplot(results)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = 3000,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-1-1'), xmax = as.Date('2022-6-30'), ymin = 0, ymax = 3000,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,3000)),
    NG = ggplot(results)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-1'), ymin = 0, ymax = 1500,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,1500))+
      coord_cartesian(ylim=c(0, 1000))
    
  )
  inc_plot <- annotations[[country]]+
    geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(aes(x=date_ex,y=inc,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    labs(x='Date',y='Clinical Incidence\nper 10,000 person-months')+
    # labs(x='Date',y='EIR')+
    facet_grid(~district)+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"))
  inc_plot
}
source('shared/addCIs_inc.R')
bf_hmis$date = as.yearmon(paste(bf_hmis$Year, bf_hmis$Month), "%Y %b")
bf_hmis_forplot <- bf_hmis[bf_hmis$district %in% c('Banfora','Gaoua','Orodara')&bf_hmis$date>=as.yearmon('Sep 2020')&
                     bf_hmis$date<=as.yearmon('Jun 2022')&!is.na(bf_hmis$Confirmed),]
bf_hmis_forplot <- addCIs_inc(bf_hmis_forplot,bf_hmis_forplot$Confirmed,bf_hmis_forplot$Population)
bf_hmis_forplot <- bf_hmis_forplot %>%
  mutate(date_ex = as.Date(date, frac = 0),
         prop_pos = Confirmed/Tested)%>%
  dplyr::rename(inc = mean)
BF_obsinc_plot <- create_obsinc_plots(bf_hmis_forplot,'BF')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/obsinc_bf_std_310223.pdf'), width = 7, height = 5)

##BF <5yo
bf_hmis_u5$date <- as.yearmon(bf_hmis_u5$date)
bf_hmis_forplot <- bf_hmis_u5[bf_hmis_u5$district %in% c('Banfora','Gaoua','Orodara')&bf_hmis_u5$date>=as.yearmon('Sep 2020')&
                             bf_hmis_u5$date<=as.yearmon('Jun 2022')&!is.na(bf_hmis_u5$positive),]
bf_hmis_forplot <- addCIs_inc(bf_hmis_forplot,bf_hmis_forplot$positive,bf_hmis_forplot$population)
bf_hmis_forplot <- bf_hmis_forplot %>%
  mutate(date_ex = as.Date(date, frac = 0),
         prop_pos = positive/tested)%>%
  dplyr::rename(inc = mean)
BF_obsinc_plot <- create_obsinc_plots(bf_hmis_forplot,'BF')

mz_hmis_filtered <- mz_hmis[mz_hmis$District %in% c('Changara','Chemba','Guro')&!is.na(mz_hmis$Positive),]%>%
  group_by(District, Month, Year, `Net type`)%>%
  dplyr::summarise(Tested=sum(Tested),
                   Positive=sum(Positive),
                   Population=sum(u5.Population))
mz_hmis_forplot <- addCIs_inc(mz_hmis_filtered,mz_hmis_filtered$Positive,mz_hmis_filtered$Population)%>%
  mutate(date_ex=as.Date(as.yearmon(paste(Year, Month), "%Y %b")),
         prop_pos = Positive/Tested)%>%
  filter(date_ex>=as.Date('2020-12-1')&date_ex<=as.Date('2022-09-1'))%>%
  dplyr::rename(district = District,
         inc = mean)

##BF >5yo

NG_obsinc_plot <- create_obsinc_plots(ng_hmis_forplot,'NG')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/obsinc_ng_std_310223.pdf'), width = 7, height = 5)

pois.daly(bf_hmis_forplot$Confirmed[1],bf_hmis_forplot$Population[1])
windows(10,10)
BF_inc_pg_threads/BF_obsinc_plot + plot_layout(guides = "collect")
NG_inc_pg_threads / NG_obsinc_plot + plot_layout(guides = "collect")

BF_CS_all_grouped_site <- readRDS('nnp/data/BF_CS_all_grouped_site_0822.rds')
MZ_CS_all_grouped_site <- readRDS('nnp/data/MZ_CS_all_grouped_site_180123.rds')
NG_CS_all_grouped_site <- readRDS('nnp/data/NG_CS_all_grouped_site.rds')

cs_data_list <- list(BF = BF_CS_all_grouped_site, MZ = MZ_CS_all_grouped_site, NG = NG_CS_all_grouped_site)
##Create dashboard by country
create_dashboard_plots <- function(results,prev_pg,prev_mg,prev_all=NULL,incidence,cs_data_list,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  district_labels <- c('Banfora - IG2','Gaoua - Standard','Orodara - PBO',
                        'Changara - PBO','Chemba - Standard','Guro - IG2',
                        'Asa - IG2','Ejigbo - Standard','Ife North - PBO','Moro - RG')
  names(district_labels) <- unlist(district_list)
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-6-1'),by='months'),
                     MZ = seq(as.Date('2020-07-1'),as.Date('2022-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2022-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  cs_data <- cs_data_list[[country]]%>%
    rename(district=site)%>%
    mutate(month=as.Date(month))
  
  df_prev <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   district = character(),
                   month = character())
  df_prev_sample <- data.frame(time = integer(),
                          value = numeric(),
                          variable = character(),
                          district = character(),
                          month = character())
  for(i in start:(start+number-1)){
    prev_history <- data.frame(t(results[[i]]$history['prev', 51:1000, -1]))
    
    
    long_prev_sum <- prev_history%>%
      mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      group_by(time)%>%
      summarise(median=median(value),
                mean=mean(value),
                upper=quantile(value,probs=0.975),
                lower=quantile(value,probs=0.025))%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df_prev <- rbind(df_prev,long_prev_sum)
    
    prev_sample <- prev_history[, sample(ncol(prev_history), 100)] %>%
      mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      mutate(district = districts[[i-start+1]],
             month = rep(dates_list[[country]],100))
    df_prev_sample <- rbind(df_prev_sample,prev_sample)
    
    
  }
  
  df_data_pg <- data.frame(t=numeric(),
                        tested=integer(),
                        positive=numeric(),
                        mean=numeric(),
                        upper=numeric(),
                        lower=numeric(),
                        district=character())
  
  for(i in start:(start+number-1)){
    data_cis <- addCIs(prev_pg[[i]],prev_pg[[i]]$positive,prev_pg[[i]]$tested)%>%
      mutate(district=names(prev_pg[i]),
             month = dates_list[[country]])
    df_data_pg <- rbind(df_data_pg,data_cis)
  }

  df_data_mg <- data.frame(t=numeric(),
                           tested=integer(),
                           positive=numeric(),
                           mean=numeric(),
                           upper=numeric(),
                           lower=numeric(),
                           district=character())
  
  for(i in start:(start+number-1)){
    data_cis <- addCIs(prev_mg[[i]],prev_mg[[i]]$positive,prev_mg[[i]]$tested)%>%
      mutate(district=names(prev_mg[i]),
             month = dates_list[[country]])
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
                   district = character(),
                   month = character())
  df_inc_sample <- data.frame(time = integer(),
                          value = numeric(),
                          variable = character(),
                          district = character(),
                          month = character())
  for(i in start:(start+number-1)){
    inc_history <- data.frame(t(results[[i]]$history['inc', 51:1000, -1]))
    
    long_inc_sum <- inc_history%>%
      dplyr::mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value)*10000*30,
                       mean=mean(value)*10000*30,
                       upper=quantile(value,probs=0.975)*10000*30,
                       lower=quantile(value,probs=0.025)*10000*30)%>%
      dplyr::mutate(district = districts[[i-start+1]],
                    month = dates_list[[country]])
    df_inc <- rbind(df_inc,long_inc_sum)
    
    inc_sample <- inc_history[, sample(ncol(inc_history), 100)] %>%
      mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value*10000*30,
                    district = districts[[i-start+1]],
                    month = rep(dates_list[[country]],100))
    df_inc_sample <- rbind(df_inc_sample,inc_sample)
  }
  
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  annotations <- list(
    BF = ggplot()+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-1'), ymin = 0, ymax = Inf,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-1'), ymin = 0, ymax = Inf,alpha = .1,fill = "#999999"),
    MZ = ggplot()+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = Inf,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-1-1'), xmax = as.Date('2022-6-30'), ymin = 0, ymax = Inf,alpha = .1,fill = "#999999"),
    NG = ggplot()+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-1'), ymin = 0, ymax = Inf,alpha = .1,fill = "#999999")
  )
  est_inc_plot <- annotations[[country]]+
    geom_line(data=df_inc_sample,aes(x=month,y=value,color=district,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_inc,aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # coord_cartesian(ylim=c(0, 1000))+
    coord_cartesian(ylim=c(0, ifelse(country=='NG',1000,max(incidence$inc)*2)))+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    labs(title = 'Estimated Incidence\nper 10,000 person-months')+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
          )
  # print('est_inc_plot')
  obs_prev_plot_mg <- annotations[[country]]+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),pch= 2, size=1.5,color="black")+
    geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color='black')+
    geom_point(data=df_data_mg,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_mg,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    scale_y_continuous(limits = c(0,1))+
    coord_cartesian(ylim=c(0, 1),xlim = range(df_data_mg$month))+
    labs(title = 'ANC Prevalence\nMultigrav')+
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
  est_prev_plot <- annotations[[country]]+
    geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),size=1.5,pch=2,color="black")+
    geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color="black")+
    geom_line(data=df_prev_sample,aes(x=month,y=value,color=district,group=variable),alpha=0.2)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_prev,aes(x=month,y=median,color=district,group=district),linewidth=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data_pg,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_pg,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits=c(0,1))+
    coord_cartesian(xlim = range(df_data_pg$month))+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
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
  obs_inc_plot <- annotations[[country]]+
    # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
    geom_line(data=incidence,aes(x=date_ex,y=inc,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    labs(title='Observed Incidence\nper 10,000 person-months')+
    coord_cartesian(ylim=c(0, ifelse(country=='NG',1000,max(incidence$inc)*2)))+
    # coord_cartesian(ylim=c(0, max(incidence$inc)*2))+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt"))
  # print('obs_inc_plot')
  obs_pos_plot <- annotations[[country]]+
    # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
    geom_line(data=incidence,aes(x=date_ex,y=prop_pos,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits=c(0,1))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    labs(title = 'Test positivity\nproportion')+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt"))
  # print('obs_pos_plot')
  # # obs_test_plot <- annotations[[country]]+
  # #   # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
  # #   geom_line(aes(x=date_ex,y=inc,color=district,group=district),size=1)+
  # #   # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
  # #   scale_color_manual(values=colors)+
  # #   scale_fill_manual(values=colors)+
  # #   scale_x_date(date_labels = "%b %Y")+
  # #   # scale_y_continuous(limits=c(0,500))+
  # #   # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  # #   # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
  # #   labs(x='Date',y='Clinical Incidence\nper 10,000 person-months')+
  # #   # labs(x='Date',y='EIR')+
  # #   facet_grid(~district)+
  # #   theme(legend.title = element_blank(),
  # #         axis.title.x = element_blank(),
  # #         axis.text.x=element_text(angle=45, hjust=1, vjust=1),
  # #         axis.ticks.x = element_line(size = 0.5), 
  # #         axis.ticks.length = unit(3, "pt"))
  est_prev_plot+obs_prev_plot_mg+est_inc_plot+obs_inc_plot+obs_pos_plot+ plot_layout(guides = "collect",ncol=5)
  # obs_prev_plot_mg
}
windows(30,15)
bf_dash <- create_dashboard_plots(results=seas_all_result_list,
                       prev_pg=nnp_pg_list,
                       prev_mg=nnp_mg_list,
                       incidence=bf_hmis_forplot,
                       cs_data_list = cs_data_list,
                       country='BF')
bf_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_bf_seas_060223.pdf'),plot=bf_dash, width = 12, height = 6)
nnp_pgmg_bulk_sifter_betaa_results[[1]][["history"]][,1,1]
mz_dash <- create_dashboard_plots(results=seas_all_result_list,
                                  prev_pg=nnp_pg_list,
                                  prev_mg=nnp_mg_list,
                                  cs_data_list = cs_data_list,
                                  incidence=mz_hmis_forplot,
                                  country='MZ')
mz_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_seas_060223.pdf'),plot = mz_dash, width = 12, height = 6)
ng_dash <- create_dashboard_plots(results=seas_all_result_list,
                                  prev_pg=nnp_pg_list,
                                  prev_mg=nnp_mg_list,
                                  cs_data_list = cs_data_list,
                                  incidence=ng_hmis_forplot,
                                  country='NG')
ng_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_seas_060223.pdf'), plot = ng_dash, width = 12, height = 6)


ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_seas_310223.pdf'), width = 7, height = 5)
results <- nnp_mgcorr_bulk_std_results
results <- nnp_pgcorr_bulk_seas_results
coefs_pg_df <- as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS'))
coefs_mg_df <- as.data.frame(readRDS('./nnp/Corr/mg_corr_sample.RDS'))
country <- 'NG'
create_dashboard_pgmg_plots <- function(results,prev_pg,prev_mg,prev_all=NULL,coefs_pg_df,coefs_mg_df,incidence,cs_data_list,country=c('BF','MZ','NG')){
  ## fn to return prevalence from log_odds
  get_prev_from_log_odds<-function(log_odds){
    return(exp(log_odds)/(1+exp(log_odds)))
  }
  
  ## fn to return odds from prevalence
  get_odds_from_prev<-function(prev){
    return(prev/(1-prev))
  }
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  admin_key <- data.frame(district=c('Banfora','Gaoua','Orodara',
                                     'Changara','Chemba','Guro',
                                     'Asa','Ejigbo','Ife North','Moro'),
                          admin=c('Cascades','Sud-Ouest','Haut-Bassins',
                                  'Tete','Sofala','Manica',
                                  'Kwara','Osun','Osun','Kwara'))
  district_labels <- c('Banfora - IG2','Gaoua - Standard','Orodara - PBO',
                       'Changara - PBO','Chemba - Standard','Guro - IG2',
                       'Asa - IG2','Ejigbo - Standard','Ife North - PBO','Moro - RG')
  names(district_labels) <- unlist(district_list)
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-6-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2022-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2022-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  cs_data <- cs_data_list[[country]]%>%
    dplyr::rename(district=site)%>%
    mutate(month=as.Date(month))
  
  df_prev <- data.frame(time = integer(),
                        median = numeric(),
                        mean = numeric(),
                        upper = numeric(),
                        lower = numeric(),
                        district = character(),
                        month = character())
  df_prev_sample <- data.frame(time = integer(),
                               value = numeric(),
                               variable = character(),
                               district = character(),
                               month = character())
 for(i in start:(start+number-1)){
    prev_history <- data.frame(t(results[[i]]$history['prev', 501:1000, -1]))
    
    
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
      dplyr::summarise(median_pg=median(prev_pg),
                mean_pg=mean(prev_pg),
                upper_pg=quantile(prev_pg,probs=0.975),
                lower_pg=quantile(prev_pg,probs=0.025),
                median_mg=median(prev_mg),
                mean_mg=mean(prev_mg),
                upper_mg=quantile(prev_mg,probs=0.975),
                lower_mg=quantile(prev_mg,probs=0.025))%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df_prev <- rbind(df_prev,long_prev_sum)
    
    prev_sample <- prev_history[, sample(ncol(prev_history), 100)] %>%
      mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      mutate(logodds_child = log(get_odds_from_prev(value)),
             prev_pg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
             }))),
             prev_mg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_mg_df$gradient[x]*(logodds_child-coefs_mg_df$av_lo_child[x])+coefs_mg_df$intercept[x])
             }))))%>%
      mutate(district = districts[[i-start+1]],
             month = rep(dates_list[[country]],100))
    df_prev_sample <- rbind(df_prev_sample,prev_sample)
    
    
  }
  
  df_data_pg <- data.frame(t=numeric(),
                           tested=integer(),
                           positive=numeric(),
                           mean=numeric(),
                           upper=numeric(),
                           lower=numeric(),
                           district=character())
  
  for(i in start:(start+number-1)){
    data_cis <- addCIs(prev_pg[[i]],prev_pg[[i]]$positive,prev_pg[[i]]$tested)%>%
      mutate(district=names(prev_pg[i]),
             month = dates_list[[country]])
    df_data_pg <- rbind(df_data_pg,data_cis)
  }
  
  df_data_mg <- data.frame(t=numeric(),
                           tested=integer(),
                           positive=numeric(),
                           mean=numeric(),
                           upper=numeric(),
                           lower=numeric(),
                           district=character())
  
  for(i in start:(start+number-1)){
    data_cis <- addCIs(prev_mg[[i]],prev_mg[[i]]$positive,prev_mg[[i]]$tested)%>%
      mutate(district=names(prev_mg[i]),
             month = dates_list[[country]])
    df_data_mg <- rbind(df_data_mg,data_cis)
  }
  
  df_data_pg$grav <- 'pg'
  df_data_mg$grav <- 'mg'
  df_data_all <- rbind(df_data_pg,df_data_mg)

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
                       district = character(),
                       month = character(),
                       rel_inc = numeric())
  df_inc_sample <- data.frame(time = integer(),
                              value = numeric(),
                              variable = character(),
                              district = character(),
                              month = character())
  for(i in start:(start+number-1)){
    inc_history <- data.frame(t(results[[i]]$history['inc', 501:1000, -1]))
    
    long_inc_sum <- inc_history%>%
      dplyr::mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value)*10000*30,
                       mean=mean(value)*10000*30,
                       upper=quantile(value,probs=0.975)*10000*30,
                       lower=quantile(value,probs=0.025)*10000*30)%>%
      dplyr::mutate(district = districts[[i-start+1]],
                    month = dates_list[[country]])%>%
      group_by(district)%>%
      mutate(rel_inc = median/max(median))
    df_inc <- rbind(df_inc,long_inc_sum)
    
    inc_sample <- inc_history[, sample(ncol(inc_history), 100)] %>%
      mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value*10000*30,
                    district = districts[[i-start+1]],
                    month = rep(dates_list[[country]],100))
    df_inc_sample <- rbind(df_inc_sample,inc_sample)
  }

  df_eir <- data.frame(time = integer(),
                       median = numeric(),
                       mean = numeric(),
                       upper = numeric(),
                       lower = numeric(),
                       district = character(),
                       month = character(),
                       rel_eir = numeric())
  df_eir_sample <- data.frame(time = integer(),
                              value = numeric(),
                              variable = character(),
                              district = character(),
                              month = character(),
                              rel_eir = numeric())
  for(i in start:(start+number-1)){
    eir_history <- data.frame(t(results[[i]]$history['EIR', 501:1000, -1]))

    long_eir_sum <- eir_history%>%
      dplyr::mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value),
                       mean=mean(value),
                       upper=quantile(value,probs=0.975),
                       lower=quantile(value,probs=0.025))%>%
      dplyr::mutate(district = districts[[i-start+1]],
                    month = dates_list[[country]],
                    rel_eir = median/max(median))
    df_eir <- rbind(df_eir,long_eir_sum)

    eir_sample <- eir_history[, sample(ncol(eir_history), 100)] %>%
      mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value,
                    district = districts[[i-start+1]],
                    month = rep(dates_list[[country]],100))
    df_eir_sample <- rbind(df_eir_sample,eir_sample)
  }

  rainfall <- read_csv('./nnp/data/NNP_rainfall.csv') %>%
    left_join(admin_key,by=join_by(NAME_1==admin))%>%
    filter(district %in% district_list[[country]])%>%
    group_by(district) %>%
    filter(between(as.yearmon(Month), as.yearmon(min(dates_list[[country]])), as.yearmon(max(dates_list[[country]]))))%>%
    mutate(rel_rainfall = Rainfall/max(Rainfall))
  
  incidence <- incidence%>%
    filter(district %in% district_list[[country]])%>%
    group_by(district)%>%
    mutate(rel_inc = (inc/max(inc)))
  
  annotations <- list(
    BF = ggplot()+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-1'), ymin = 0, ymax = Inf,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-1'), ymin = 0, ymax = Inf,alpha = .1,fill = "#999999"),
    MZ = ggplot()+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = Inf,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-1-1'), xmax = as.Date('2022-6-30'), ymin = 0, ymax = Inf,alpha = .1,fill = "#999999"),
    NG = ggplot()+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-1'), ymin = 0, ymax = Inf,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-7-1'), xmax = as.Date('2022-11-1'), ymin = 0, ymax = Inf,alpha = .1,fill = "#999999")
  )
  est_eir_plot <- annotations[[country]]+
    # geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
    geom_line(data=df_eir_sample,aes(x=as.Date(month),y=value,color=district,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_eir,aes(x=as.Date(month),y=median,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # coord_cartesian(ylim=c(0, 1000))+
    # scale_y_continuous(limits=c(0,500))+
    scale_y_log10(breaks=c(.01,.1,1,10,100,1000),labels=c(.01,.1,1,10,100,1000))+
    coord_cartesian(ylim=c(.01, 800))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    labs(title = 'EIR - log scale')+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  rel_eir_plot <- annotations[[country]]+
    geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
    # geom_line(data=df_eir_sample,aes(x=as.Date(month),y=rel_eir,color=site,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_eir,aes(x=as.Date(month),y=rel_eir,color=district,group=district),size=1)+
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
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    labs(title = 'Relative EIR and\nRainfall')+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  
  est_inc_plot <- annotations[[country]]+
    geom_line(data=df_inc_sample,aes(x=month,y=value,color=district,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_inc,aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # coord_cartesian(ylim=c(0, 1000))+
    coord_cartesian(ylim=c(0, ifelse(country=='NG',1000,max(incidence$inc)*2)))+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    labs(title = 'Estimated Incidence\nper 10,000 person-months')+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  rel_est_inc_plot <- annotations[[country]]+
    geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
    # geom_line(data=df_inc_sample,aes(x=month,y=value,color=district,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_inc,aes(x=month,y=rel_inc,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # coord_cartesian(ylim=c(0, 1000))+
    # coord_cartesian(ylim=c(0, ifelse(country=='NG',1000,max(incidence$inc)*2)))+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    labs(title = 'Rel. Est. Cases\nand Rainfall')+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  # print('est_inc_plot')
  obs_prev_plot_mg <- annotations[[country]]+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data_mg,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_mg,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),pch= 2, size=1.5,color="black")+
    geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color='black')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    scale_y_continuous(limits = c(0,1))+
    coord_cartesian(ylim=c(0, 1),xlim = range(df_data_mg$month))+
    labs(title = 'ANC Prevalence\nMultigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  est_prev_plot_mg <- annotations[[country]]+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_line(data=df_prev_sample,aes(x=month,y=prev_mg,color=district,group=variable),alpha=0.2)+
    geom_line(data=df_prev,aes(x=month,y=median_mg,color=district,group=district),linewidth=1)+
    geom_point(data=df_data_mg,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_mg,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),pch= 2, size=1.5,color="black")+
    geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color='black')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    scale_y_continuous(limits = c(0,1))+
    coord_cartesian(ylim=c(0, 1),xlim = range(df_data_mg$month))+
    labs(title = 'ANC Prevalence\nMultigrav')+
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
  est_prev_plot_pg <- annotations[[country]]+
    geom_line(data=df_prev_sample,aes(x=month,y=prev_pg,color=district,group=variable),alpha=0.2)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_prev,aes(x=month,y=median_pg,color=district,group=district),linewidth=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data_pg,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_pg,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),size=1.5,pch=2,color="black")+
    geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color="black")+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits=c(0,1))+
    coord_cartesian(xlim = range(df_data_pg$month))+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    labs(title = 'ANC Prevalence\nPrimigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  obs_prev_plot_pg <- annotations[[country]]+
    geom_point(data=df_data_pg,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_pg,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),size=1.5,pch=2,color="black")+
    geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color="black")+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits=c(0,1))+
    coord_cartesian(xlim = range(df_data_pg$month))+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
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
  obs_inc_plot <- annotations[[country]]+
    # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
    geom_line(data=incidence,aes(x=date_ex,y=inc,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    labs(title='Observed Incidence\nper 10,000 person-months')+
    coord_cartesian(ylim=c(0, ifelse(country=='NG',1000,max(incidence$inc)*2)))+
    # coord_cartesian(ylim=c(0, max(incidence$inc)*2))+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt"))
  obs_rel_inc_plot <- annotations[[country]]+
    geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
    # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
    geom_line(data=incidence,aes(x=date_ex,y=rel_inc,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    labs(title='Rel. Obs. Cases\nand Rainfall')+
    # coord_cartesian(ylim=c(0, ifelse(country=='NG',1000,max(incidence$inc)*2)))+
    # coord_cartesian(ylim=c(0, max(incidence$inc)*2))+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt"))
  
    # print('obs_inc_plot')
  obs_pos_plot <- annotations[[country]]+
    # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
    geom_line(data=incidence,aes(x=date_ex,y=prop_pos,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits=c(0,1))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    labs(title = 'Test positivity\nproportion')+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt"))
  sample_size_plot <- annotations[[country]]+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=site,group=site),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=site,group=site),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_line(data=df_data_all,aes(x=month,y=tested,color=district,linetype=factor(grav,levels=c('pg','mg'))))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
    labs(title = 'Number Tested')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  # print('obs_pos_plot')
  # # obs_test_plot <- annotations[[country]]+
  # #   # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
  # #   geom_line(aes(x=date_ex,y=inc,color=district,group=district),size=1)+
  # #   # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
  # #   scale_color_manual(values=colors)+
  # #   scale_fill_manual(values=colors)+
  # #   scale_x_date(date_labels = "%b %Y")+
  # #   # scale_y_continuous(limits=c(0,500))+
  # #   # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  # #   # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
  # #   labs(x='Date',y='Clinical Incidence\nper 10,000 person-months')+
  # #   # labs(x='Date',y='EIR')+
  # #   facet_grid(~district)+
  # #   theme(legend.title = element_blank(),
  # #         axis.title.x = element_blank(),
  # #         axis.text.x=element_text(angle=45, hjust=1, vjust=1),
  # #         axis.ticks.x = element_line(size = 0.5), 
  # #         axis.ticks.length = unit(3, "pt"))
  full_dash <- est_prev_plot_pg+est_prev_plot_mg+est_inc_plot+obs_inc_plot+obs_pos_plot+ plot_layout(guides = "collect",ncol=5)
  obs_data_dash <- obs_prev_plot_pg+obs_prev_plot_mg+est_inc_plot+obs_inc_plot+ plot_layout(guides = "collect",ncol=4)
  pres_dash <- est_prev_plot_pg+est_prev_plot_mg+est_inc_plot+obs_inc_plot+ plot_layout(guides = "collect",ncol=4)
  rel_dash <- sample_size_plot+est_prev_plot_pg+est_prev_plot_mg+est_eir_plot+rel_est_inc_plot+obs_rel_inc_plot+plot_layout(guides = "collect",ncol=6)
  return(list(full_dash = full_dash,
              obs_data_dash = obs_data_dash,
              pres_dash = pres_dash,
              rel_dash = rel_dash))
  # obs_prev_plot_mg
}
coefs_pg_df <- as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS'))
coefs_mg_df <- as.data.frame(readRDS('./nnp/Corr/mg_corr_sample.RDS'))


windows(30,15)
bf_pgmg_dash <- create_dashboard_pgmg_plots(results=nnp_mgcorr_bulk_std_results,
                                  prev_pg=nnp_pg_list,
                                  prev_mg=nnp_mg_list,
                                  coefs_pg_df = coefs_pg_df,
                                  coefs_mg_df = coefs_mg_df,
                                  incidence=bf_hmis_forplot,
                                  cs_data_list = cs_data_list,
                                  country='BF')
bf_pgmg_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_bf_pgmg_std_220223.pdf'),plot=bf_pgmg_dash$full_dash, width = 12, height = 6)
bf_pgmg_dash$obs_data_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_bf_pgmg_std_220223_obs.pdf'),plot=bf_pgmg_dash$obs_data_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_bf_pgmg_std_220223_pres.pdf'),plot=bf_pgmg_dash$pres_dash, width = 12, height = 6)

mz_pgmg_dash <- create_dashboard_pgmg_plots(results=nnp_mgcorr_bulk_std_results,
                                            prev_pg=nnp_pg_list,
                                            prev_mg=nnp_mg_list,
                                            coefs_pg_df = coefs_pg_df,
                                            coefs_mg_df = coefs_mg_df,
                                            cs_data_list = cs_data_list,
                                            incidence=mz_hmis_forplot,
                                            country='MZ')
mz_pgmg_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_pgmg_std_220223.pdf'),plot = mz_pgmg_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_pgmg_std_220223_obs.pdf'),plot = mz_pgmg_dash$obs_data_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_pgmg_std_220223_pres.pdf'),plot = mz_pgmg_dash$pres_dash, width = 12, height = 6)

ng_pgmg_dash <- create_dashboard_pgmg_plots(results=nnp_mgcorr_bulk_std_results,
                                  prev_pg=nnp_pg_list,
                                  prev_mg=nnp_mg_list,
                                  coefs_pg_df = coefs_pg_df,
                                  coefs_mg_df = coefs_mg_df,
                                  cs_data_list = cs_data_list,
                                  incidence=ng_hmis_forplot,
                                  country='NG')
ng_pgmg_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pgmg_std_220223.pdf'), plot = ng_pgmg_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pgmg_std_220223_obs.pdf'), plot = ng_pgmg_dash$obs_data_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pgmg_std_220223_pres.pdf'), plot = ng_pgmg_dash$pres_dash, width = 12, height = 6)

##PG only - Seasonal model

windows(30,15)
bf_pg_seas_dash <- create_dashboard_pgmg_plots(results=nnp_pgcorr_bulk_seas_results,
                                            prev_pg=nnp_pg_list,
                                            prev_mg=nnp_mg_list,
                                            coefs_pg_df = coefs_pg_df,
                                            coefs_mg_df = coefs_mg_df,
                                            incidence=bf_hmis_forplot,
                                            cs_data_list = cs_data_list,
                                            country='BF')
bf_pg_seas_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_bf_pg_seas_270223.pdf'),plot=bf_pg_seas_dash$full_dash, width = 12, height = 6)

mz_pg_seas_dash <- create_dashboard_pgmg_plots(results=nnp_pgcorr_bulk_seas_results,
                                            prev_pg=nnp_pg_list,
                                            prev_mg=nnp_mg_list,
                                            coefs_pg_df = coefs_pg_df,
                                            coefs_mg_df = coefs_mg_df,
                                            cs_data_list = cs_data_list,
                                            incidence=mz_hmis_forplot,
                                            country='MZ')
mz_pg_seas_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_pg_seas_270223.pdf'),plot = mz_pg_seas_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_pg_seas_270223_obs.pdf'),plot = mz_pg_seas_dash$obs_data_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_pg_seas_270223_pres.pdf'),plot = mz_pg_seas_dash$pres_dash, width = 12, height = 6)

ng_pg_seas_dash <- create_dashboard_pgmg_plots(results=nnp_pgcorr_bulk_seas_results_update,
                                            prev_pg=nnp_pg_list,
                                            prev_mg=nnp_mg_list,
                                            coefs_pg_df = coefs_pg_df,
                                            coefs_mg_df = coefs_mg_df,
                                            cs_data_list = cs_data_list,
                                            incidence=ng_hmis_forplot,
                                            country='NG')
ng_pg_seas_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pg_seas_270223.pdf'), plot = ng_pg_seas_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pg_seas_270223_obs.pdf'), plot = ng_pg_seas_dash$obs_data_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pg_seas_270223_pres.pdf'), plot = ng_pg_seas_dash$pres_dash, width = 12, height = 6)

##PG only - Standard model

windows(30,15)
bf_pg_std_dash <- create_dashboard_pgmg_plots(results=nnp_pgcorr_bulk_std_results,
                                               prev_pg=nnp_pg_list,
                                               prev_mg=nnp_mg_list,
                                               coefs_pg_df = coefs_pg_df,
                                               coefs_mg_df = coefs_mg_df,
                                               incidence=bf_hmis_forplot,
                                               cs_data_list = cs_data_list,
                                               country='BF')
bf_pg_std_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_bf_pg_std_280223.pdf'),plot=bf_pg_std_dash$full_dash, width = 12, height = 6)

mz_pg_std_dash <- create_dashboard_pgmg_plots(results=nnp_pgcorr_bulk_std_results,
                                               prev_pg=nnp_pg_list,
                                               prev_mg=nnp_mg_list,
                                               coefs_pg_df = coefs_pg_df,
                                               coefs_mg_df = coefs_mg_df,
                                               cs_data_list = cs_data_list,
                                               incidence=mz_hmis_forplot,
                                               country='MZ')
mz_pg_std_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_pg_std_280223.pdf'),plot = mz_pg_std_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_pg_std_280223_obs.pdf'),plot = mz_pg_std_dash$obs_data_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_pg_std_280223_pres.pdf'),plot = mz_pg_std_dash$pres_dash, width = 12, height = 6)

ng_pg_std_dash <- create_dashboard_pgmg_plots(results=nnp_pgcorr_bulk_std_results,
                                               prev_pg=nnp_pg_list,
                                               prev_mg=nnp_mg_list,
                                               coefs_pg_df = coefs_pg_df,
                                               coefs_mg_df = coefs_mg_df,
                                               cs_data_list = cs_data_list,
                                               incidence=ng_hmis_forplot,
                                               country='NG')
ng_pg_std_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pg_std_280223.pdf'), plot = ng_pg_std_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pg_std_280223_obs.pdf'), plot = ng_pg_std_dash$obs_data_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pg_std_280223_pres.pdf'), plot = ng_pg_std_dash$pres_dash, width = 12, height = 6)

windows(30,15)
bf_pgmg_seas_dash <- create_dashboard_pgmg_plots(results=nnp_mgcorr_bulk_seas_results_update,
                                               prev_pg=nnp_pg_list,
                                               prev_mg=nnp_mg_list,
                                               coefs_pg_df = coefs_pg_df,
                                               coefs_mg_df = coefs_mg_df,
                                               incidence=bf_hmis_forplot,
                                               cs_data_list = cs_data_list,
                                               country='BF')
bf_pgmg_seas_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_bf_pgmg_seas_300523.pdf'),plot=bf_pgmg_seas_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_bf_pgmg_seas_300523_pres.pdf'),plot=bf_pgmg_seas_dash$pres_dash, width = 12, height = 6)
bf_pgmg_seas_dash$rel_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_bf_pgmg_seas_300523_rel.pdf'),plot=bf_pgmg_seas_dash$rel_dash, width = 12, height = 6)

mz_pgmg_seas_dash <- create_dashboard_pgmg_plots(results=nnp_mgcorr_bulk_seas_results_update,
                                               prev_pg=nnp_pg_list,
                                               prev_mg=nnp_mg_list,
                                               coefs_pg_df = coefs_pg_df,
                                               coefs_mg_df = coefs_mg_df,
                                               cs_data_list = cs_data_list,
                                               incidence=mz_hmis_forplot,
                                               country='MZ')
mz_pgmg_seas_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_pgmg_seas_300523.pdf'),plot = mz_pgmg_seas_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_pgmg_seas_300523_pres.pdf'),plot = mz_pgmg_seas_dash$pres_dash, width = 12, height = 6)
mz_pgmg_seas_dash$rel_dash

ng_pgmg_seas_dash <- create_dashboard_pgmg_plots(results=nnp_mgcorr_bulk_seas_results_update,
                                              prev_pg=nnp_pg_list,
                                              prev_mg=nnp_mg_list,
                                              coefs_pg_df = coefs_pg_df,
                                              coefs_mg_df = coefs_mg_df,
                                              cs_data_list = cs_data_list,
                                              incidence=ng_hmis_forplot,
                                              country='NG')
ng_pgmg_seas_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pgmg_seas_080323.pdf'), plot = ng_pgmg_seas_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pgmg_seas_080323_pres.pdf'), plot = ng_pgmg_seas_dash$pres_dash, width = 12, height = 6)
ng_pgmg_seas_dash$rel_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pgmg_seas_080323_rel.pdf'), plot = ng_pgmg_seas_dash$rel_dash, width = 12, height = 6)

source('./nnp/create_plots_nnp.R')
windows(30,15)
names(nnp_pgmg_bulk_sifter_eir_results) <- names(nnp_pg_list)
bf_pgmg_seas_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_eir_results,
                                                prev_pg=nnp_pg_list,
                                                prev_mg=nnp_mg_list,
                                                coefs_pg_df = coefs_pg_df,
                                                coefs_mg_df = coefs_mg_df,
                                                incidence=bf_hmis_forplot,
                                                cs_data_list = cs_data_list,
                                                country='BF')
bf_pgmg_seas_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_bf_pgmg_sifter_eir.pdf'), plot = bf_pgmg_seas_dash$full_dash, width = 12, height = 6)
bf_pgmg_seas_dash$rel_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/rel_dash_bf_pgmg_sifter_eir.pdf'), plot = bf_pgmg_seas_dash$rel_dash, width = 12, height = 6)
bf_pgmg_seas_dash$mcmc_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/mcmc_dash_bf_pgmg_sifter_eir.pdf'), plot = bf_pgmg_seas_dash$full_dash, width = 12, height = 6)

bf_pgmg_seas_betaa_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_betaa_results,
                                                prev_pg=nnp_pg_list,
                                                prev_mg=nnp_mg_list,
                                                coefs_pg_df = coefs_pg_df,
                                                coefs_mg_df = coefs_mg_df,
                                                incidence=bf_hmis_forplot,
                                                cs_data_list = cs_data_list,
                                                country='BF',
                                                param='betaa')
bf_pgmg_seas_betaa_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_bf_pgmg_sifter_betaa.pdf'), plot = bf_pgmg_seas_betaa_dash$full_dash, width = 12, height = 6)
bf_pgmg_seas_betaa_dash$rel_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/rel_dash_bf_pgmg_sifter_betaa.pdf'), plot = bf_pgmg_seas_betaa_dash$rel_dash, width = 12, height = 6)
bf_pgmg_seas_betaa_dash$mcmc_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/mcmc_dash_bf_pgmg_sifter_betaa.pdf'), plot = bf_pgmg_seas_betaa_dash$full_dash, width = 12, height = 6)


mz_pgmg_seas_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_eir_results,
                                                                                 prev_pg=nnp_pg_list,
                                                                                 prev_mg=nnp_mg_list,
                                                                                 coefs_pg_df = coefs_pg_df,
                                                                                 coefs_mg_df = coefs_mg_df,
                                                                                 incidence=mz_hmis_forplot,
                                                                                 cs_data_list = cs_data_list,
                                                                                 country='MZ')
mz_pgmg_seas_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_pgmg_sifter_eir.pdf'), plot = mz_pgmg_seas_dash$full_dash, width = 12, height = 6)
mz_pgmg_seas_dash$rel_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/rel_dash_mz_pgmg_sifter_eir.pdf'), plot = mz_pgmg_seas_dash$rel_dash, width = 12, height = 6)
mz_pgmg_seas_dash$mcmc_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/mcmc_dash_mz_pgmg_sifter_eir.pdf'), plot = mz_pgmg_seas_dash$mcmc_dash, width = 12, height = 6)


mz_pgmg_seas_betaa_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_betaa_results,
                                                      prev_pg=nnp_pg_list,
                                                      prev_mg=nnp_mg_list,
                                                      coefs_pg_df = coefs_pg_df,
                                                      coefs_mg_df = coefs_mg_df,
                                                      incidence=mz_hmis_forplot,
                                                      cs_data_list = cs_data_list,
                                                      country='MZ',
                                                      param='betaa')
mz_pgmg_seas_betaa_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_mz_pgmg_sifter_betaa.pdf'), plot = mz_pgmg_seas_betaa_dash$full_dash, width = 12, height = 6)
mz_pgmg_seas_betaa_dash$rel_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/rel_dash_mz_pgmg_sifter_betaa.pdf'), plot = mz_pgmg_seas_betaa_dash$rel_dash, width = 12, height = 6)
mz_pgmg_seas_betaa_dash$mcmc_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/mcmc_dash_mz_pgmg_sifter_betaa.pdf'), plot = mz_pgmg_seas_betaa_dash$mcmc_dash, width = 12, height = 6)

ng_pgmg_seas_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_eir_results,
                                                prev_pg=nnp_pg_list,
                                                prev_mg=nnp_mg_list,
                                                coefs_pg_df = coefs_pg_df,
                                                coefs_mg_df = coefs_mg_df,
                                                incidence=ng_hmis_forplot,
                                                cs_data_list = cs_data_list,
                                                country='NG')
ng_pgmg_seas_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pgmg_sifter_eir.pdf'), plot = ng_pgmg_seas_dash$full_dash, width = 12, height = 6)
ng_pgmg_seas_dash$rel_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/rel_dash_ng_pgmg_sifter_eir.pdf'), plot = ng_pgmg_seas_dash$rel_dash, width = 12, height = 6)
ng_pgmg_seas_dash$mcmc_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/mcmc_dash_ng_pgmg_sifter_eir.pdf'), plot = ng_pgmg_seas_dash$mcmc_dash, width = 12, height = 6)

ng_pgmg_seas_betaa_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_betaa_results,
                                                      prev_pg=nnp_pg_list,
                                                      prev_mg=nnp_mg_list,
                                                      coefs_pg_df = coefs_pg_df,
                                                      coefs_mg_df = coefs_mg_df,
                                                      incidence=ng_hmis_forplot,
                                                      cs_data_list = cs_data_list,
                                                      country='NG',
                                                      param='betaa')
ng_pgmg_seas_betaa_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/dash_ng_pgmg_sifter_betaa.pdf'), plot = ng_pgmg_seas_betaa_dash$full_dash, width = 12, height = 6)
ng_pgmg_seas_betaa_dash$rel_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/rel_dash_ng_pgmg_sifter_betaa.pdf'), plot = ng_pgmg_seas_betaa_dash$rel_dash, width = 12, height = 6)
ng_pgmg_seas_betaa_dash$mcmc_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/mcmc_dash_ng_pgmg_sifter_betaa.pdf'), plot = ng_pgmg_seas_betaa_dash$mcmc_dash, width = 12, height = 6)
nnp_pgmg_bulk_sifter_betaa_results[[7]]

ratio <- nnp_pgmg_bulk_sifter_betaa_results[[1]]$history['EIR',600,]/nnp_pgmg_bulk_sifter_betaa_results[[1]]$history['betaa',600,]
plot(ratio)
ratio[1:15]
plot(ratio[200:400])

###tests###
bf_pgmg_seas_betaa_test1_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_betaa_test1_results,
                                                      prev_pg=nnp_pg_list,
                                                      prev_mg=nnp_mg_list,
                                                      coefs_pg_df = coefs_pg_df,
                                                      coefs_mg_df = coefs_mg_df,
                                                      incidence=bf_hmis_forplot,
                                                      cs_data_list = cs_data_list,
                                                      burnin = 0,
                                                      sample_size = 50,
                                                      chain_length = 50,
                                                      country='BF',
                                                      param='betaa')
windows(12,6)
names(nnp_pg_list) <- names(nnp_mg_list) <- names(nnp_pgmg_bulk_sifter_betaa_test1_results) <- district_list <- c('Banfora','Gaoua','Orodara','Changara','Chemba','Guro','Asa','Ejigbo','Ife North','Moro')

bf_pgmg_seas_betaa_test1_dash$rel_dash
bf_pgmg_seas_betaa_test1_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_betaa_test1_results,
                                                            prev_pg=nnp_pg_list,
                                                            prev_mg=nnp_mg_list,
                                                            coefs_pg_df = coefs_pg_df,
                                                            coefs_mg_df = coefs_mg_df,
                                                            incidence=bf_hmis_forplot,
                                                            cs_data_list = cs_data_list,
                                                            burnin = 0,
                                                            sample_size = 50,
                                                            chain_length = 50,
                                                            country='BF',
                                                            param='betaa')
windows(12,6)
mz_pgmg_seas_betaa_test1_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_betaa_test1_results,
                                                            prev_pg=nnp_pg_list,
                                                            prev_mg=nnp_mg_list,
                                                            coefs_pg_df = coefs_pg_df,
                                                            coefs_mg_df = coefs_mg_df,
                                                            incidence=mz_hmis_forplot,
                                                            cs_data_list = cs_data_list,
                                                            burnin = 0,
                                                            sample_size = 50,
                                                            chain_length = 50,
                                                            country='MZ',
                                                            param='betaa')
windows(12,6)
mz_pgmg_seas_betaa_test1_dash$rel_dash

ng_pgmg_seas_betaa_test1_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_betaa_test1_results,
                                                            prev_pg=nnp_pg_list,
                                                            prev_mg=nnp_mg_list,
                                                            coefs_pg_df = coefs_pg_df,
                                                            coefs_mg_df = coefs_mg_df,
                                                            incidence=ng_hmis_forplot,
                                                            cs_data_list = cs_data_list,
                                                            burnin = 0,
                                                            sample_size = 50,
                                                            chain_length = 50,
                                                            country='NG',
                                                            param='betaa')
windows(12,6)
ng_pgmg_seas_betaa_test1_dash$rel_dash

nnp_pgmg_bulk_sifter_betaa_max125_results

bf_pgmg_sifter_betaa_max125_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_betaa_max125_results,
                                                               prev_pg=nnp_pg_list,
                                                               prev_mg=nnp_mg_list,
                                                               coefs_pg_df = coefs_pg_df,
                                                               coefs_mg_df = coefs_mg_df,
                                                               incidence=bf_hmis_forplot,
                                                               cs_data_list = cs_data_list,
                                                               burnin = 100,
                                                               sample_size = 100,
                                                               chain_length = 1000,
                                                               country='BF',
                                                               param='betaa')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/rel_dash_bf_pgmg_sifter_betaa_max125.pdf'), plot = bf_pgmg_sifter_betaa_max125_dash$rel_dash, width = 12, height = 6)
bf_pgmg_sifter_betaa_max125_dash$mcmc_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/full_dash_bf_pgmg_sifter_betaa_max125.pdf'), plot = bf_pgmg_sifter_betaa_max125_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/mcmc_dash_bf_pgmg_sifter_betaa_max125.pdf'), plot = bf_pgmg_sifter_betaa_max125_dash$mcmc_dash, width = 6, height = 6)

mz_pgmg_sifter_betaa_max125_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_betaa_max125_results,
                                                               prev_pg=nnp_pg_list,
                                                               prev_mg=nnp_mg_list,
                                                               coefs_pg_df = coefs_pg_df,
                                                               coefs_mg_df = coefs_mg_df,
                                                               incidence=mz_hmis_forplot,
                                                               cs_data_list = cs_data_list,
                                                               burnin = 100,
                                                               sample_size = 100,
                                                               chain_length = 1000,
                                                               country='MZ',
                                                               param='betaa')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/rel_dash_mz_pgmg_sifter_betaa_max125.pdf'), plot = mz_pgmg_sifter_betaa_max125_dash$rel_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/full_dash_mz_pgmg_sifter_betaa_max125.pdf'), plot = mz_pgmg_sifter_betaa_max125_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/mcmc_dash_mz_pgmg_sifter_betaa_max125.pdf'), plot = mz_pgmg_sifter_betaa_max125_dash$mcmc_dash, width = 6, height = 6)
ng_pgmg_sifter_betaa_max125_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_betaa_max125_results,
                                                            prev_pg=nnp_pg_list,
                                                            prev_mg=nnp_mg_list,
                                                            coefs_pg_df = coefs_pg_df,
                                                            coefs_mg_df = coefs_mg_df,
                                                            incidence=ng_hmis_forplot,
                                                            cs_data_list = cs_data_list,
                                                            burnin = 100,
                                                            sample_size = 100,
                                                            chain_length = 1000,
                                                            country='NG',
                                                            param='betaa')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/rel_dash_ng_pgmg_sifter_betaa_2407.pdf'), plot = ng_pgmg_sifter_betaa_2407_dash$rel_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/full_dash_ng_pgmg_sifter_betaa_2407.pdf'), plot = ng_pgmg_sifter_betaa_2407_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/mcmc_dash_ng_pgmg_sifter_betaa_2407.pdf'), plot = ng_pgmg_sifter_betaa_2407_dash$mcmc_dash, width = 6, height = 6)
windows(12,6)
ng_pgmg_seas_betaa_test1_dash$rel_dash

names(nnp_pgmg_bulk_sifter_betaa_2407_results) <- district_list
bf_pgmg_sifter_betaa_2407_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_betaa_2407_results,
                                                               prev_pg=nnp_pg_list,
                                                               prev_mg=nnp_mg_list,
                                                               coefs_pg_df = coefs_pg_df,
                                                               coefs_mg_df = coefs_mg_df,
                                                               incidence=bf_hmis_forplot,
                                                               cs_data_list = cs_data_list,
                                                               burnin = 100,
                                                               sample_size = 100,
                                                               chain_length = 1000,
                                                               country='BF',
                                                               param='betaa')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/rel_dash_bf_pgmg_sifter_betaa_2407.pdf'), plot = bf_pgmg_sifter_betaa_2407_dash$rel_dash, width = 12, height = 6)
windows(12,6)
bf_pgmg_sifter_betaa_2407_dash$full_dash

bf_pgmg_sifter_betaa_2407_dash$mcmc_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/full_dash_bf_pgmg_sifter_betaa_2407.pdf'), plot = bf_pgmg_sifter_betaa_2407_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/mcmc_dash_bf_pgmg_sifter_betaa_2407.pdf'), plot = bf_pgmg_sifter_betaa_2407_dash$mcmc_dash, width = 6, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/overlay_bf_pgmg_sifter_betaa_2407.pdf'), plot = bf_pgmg_sifter_betaa_2407_dash$overlay_plot, width = 12, height = 6)

bf_pgmg_sifter_betaa_2407_dash$overlay_plot

mz_pgmg_sifter_betaa_2407_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_betaa_2407_results,
                                                               prev_pg=nnp_pg_list,
                                                               prev_mg=nnp_mg_list,
                                                               coefs_pg_df = coefs_pg_df,
                                                               coefs_mg_df = coefs_mg_df,
                                                               incidence=mz_hmis_forplot,
                                                               cs_data_list = cs_data_list,
                                                               burnin = 100,
                                                               sample_size = 100,
                                                               chain_length = 1000,
                                                               country='MZ',
                                                               param='betaa')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/rel_dash_mz_pgmg_sifter_betaa_2407.pdf'), plot = mz_pgmg_sifter_betaa_2407_dash$rel_dash, width = 12, height = 6)
mz_pgmg_sifter_betaa_2407_dash$full_dash
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/full_dash_mz_pgmg_sifter_betaa_2407.pdf'), plot = mz_pgmg_sifter_betaa_2407_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/mcmc_dash_mz_pgmg_sifter_betaa_2407.pdf'), plot = mz_pgmg_sifter_betaa_2407_dash$mcmc_dash, width = 6, height = 6)
ng_pgmg_sifter_betaa_2407_dash <- create_dashboard_plots_nnp(results=nnp_pgmg_bulk_sifter_betaa_2407_results,
                                                               prev_pg=nnp_pg_list,
                                                               prev_mg=nnp_mg_list,
                                                               coefs_pg_df = coefs_pg_df,
                                                               coefs_mg_df = coefs_mg_df,
                                                               incidence=ng_hmis_forplot,
                                                               cs_data_list = cs_data_list,
                                                               burnin = 100,
                                                               sample_size = 100,
                                                               chain_length = 1000,
                                                               country='NG',
                                                               param='betaa')
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/rel_dash_ng_pgmg_sifter_betaa_2407.pdf'), plot = ng_pgmg_sifter_betaa_2407_dash$rel_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/full_dash_ng_pgmg_sifter_betaa_2407.pdf'), plot = ng_pgmg_sifter_betaa_2407_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/nnp/figures/obs_prev/mcmc_dash_ng_pgmg_sifter_betaa_2407.pdf'), plot = ng_pgmg_sifter_betaa_2407_dash$mcmc_dash, width = 6, height = 6)
