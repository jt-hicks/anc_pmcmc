create_dashboard_plots_trial <- function(results,
                                       prev_pg,
                                       prev_mg,
                                       prev_all=NULL,
                                       coefs_pg_df,
                                       coefs_mg_df,
                                       start_pf_time = 30*4,
                                       param = 'EIR'){
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
  
  
  dates_list <- lapply(1:length(prev_pg),
                       function(x,start_pf_time){
                         start_obs <- min(zoo::as.Date(zoo::as.yearmon((prev_pg[[x]]$month))))#Month of first observation (in Date format)
                         time_origin <- as.Date(paste0(year(start_obs)-1,'-01-01')) #January 1 of year before observation (in Date format)
                         date <- zoo::as.Date(zoo::as.yearmon(prev_pg[[x]]$month), frac = 0.5) #Convert dates to middle of month
                         t <- as.integer(difftime(date,time_origin,units="days"))
                         
                         return(data.frame(date=date,t=t))
                       })
  
  start_times <- sapply(dates_list, function(x,start_pf_time){
    return(list(t=min(x$t)-start_pf_time,date=min(x$date)-start_pf_time))
  },start_pf_time=start_pf_time)
  
  sample_index <- sample(c(1:900), 100)
  
  mcmc.df <- bind_rows(lapply(1:4, 
                              function(x){
                                df <- results[[x]]$mcmc[101:1000,]
                                df$sites <- names(results[x])
                                df$step <- 1:nrow(df)
                                return(df)
                              }))
  
  df_prev <- bind_rows(lapply(1:4,function(i){
    prev_history <- data.frame(t(results[[i]]$history['prev_05', 101:1000, -1]))
    prev_history$t <- results[[i]]$times[-1]
    prev_history <- left_join(dates_list[[i]],prev_history,by=join_by(t==t))
    long_prev_sum <- prev_history%>%
      reshape2::melt(id=c('t','date'))%>%
      dplyr::rename(time=t)%>%
      mutate(logodds_child = log(get_odds_from_prev(value)),
             prev_pg = matrixStats::rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
             }))),
             prev_mg = matrixStats::rowMedians(as.matrix(sapply(1:nrow(coefs_mg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_mg_df$gradient[x]*(logodds_child-coefs_mg_df$av_lo_child[x])+coefs_mg_df$intercept[x])
             }))))%>%
      group_by(time,date)%>%
      dplyr::summarise(median_pg=median(prev_pg),
                       mean_pg=mean(prev_pg),
                       upper_pg=quantile(prev_pg,na.rm=TRUE,probs=0.975),
                       lower_pg=quantile(prev_pg,na.rm=TRUE,probs=0.025),
                       median_mg=median(prev_mg),
                       mean_mg=mean(prev_mg),
                       upper_mg=quantile(prev_mg,na.rm=TRUE,probs=0.975),
                       lower_mg=quantile(prev_mg,na.rm=TRUE,probs=0.025))%>%
      mutate(district = districts[[i]])
    long_prev_sum$std_prev_pg = (long_prev_sum$mean_pg - mean(long_prev_sum$mean_pg))/sd(long_prev_sum$mean_pg)
    long_prev_sum$std_prev_mg = (long_prev_sum$mean_mg - mean(long_prev_sum$mean_mg))/sd(long_prev_sum$mean_mg)
    return(long_prev_sum)
  }))
  df_prev_sample <- bind_rows(lapply(1:4,function(i){
    prev_history <- data.frame(t(results[[i]]$history['prev_05', 101:1000, -1]))
    prev_history <- prev_history[, sample_index]
    prev_history$t <- results[[i]]$times[-1]
    prev_history <- left_join(dates_list[[i]],prev_history,by=join_by(t==t))
    prev_sample <- prev_history%>%
      reshape2::melt(id=c('t','date'))%>%
      dplyr::rename(time=t)%>%
      mutate(logodds_child = log(get_odds_from_prev(value)),
             prev_pg = matrixStats::rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
             }))),
             prev_mg = matrixStats::rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_mg_df$gradient[x]*(logodds_child-coefs_mg_df$av_lo_child[x])+coefs_mg_df$intercept[x])
             }))))%>%
      mutate(district = districts[[i]])
    return(prev_sample)
  }))
  
  observed.pg <- bind_rows(lapply(1:4,
                                  function(x){
                                    df <- prev_pg[[x]]
                                    df$district <- names(prev_pg[x])
                                    return(df)
                                  }))
  observed.pg <- addCIs(observed.pg,observed.pg$positive,observed.pg$tested)
  observed.pg$date <- as.Date(observed.pg$month,frac = 0.5)

  observed.mg <- bind_rows(lapply(1:4,
                                  function(x){
                                    df <- prev_mg[[x]]
                                    df$district <- names(prev_mg[x])
                                    return(df)
                                  }))
  observed.mg <- addCIs(observed.mg,observed.mg$positive,observed.mg$tested)
  observed.mg$date <- as.Date(observed.mg$month,frac = 0.5)
  
  observed.pg$grav <- 'pg'
  observed.mg$grav <- 'mg'
  observed.all <- rbind(observed.pg,observed.mg)
  
  df_inc <- bind_rows(lapply(1:4,function(i){
    inc_history <- data.frame(t(results[[i]]$history['clininc_05', 101:1000, -1]))
    inc_history$t <- results[[i]]$times[-1]
    inc_history <- left_join(dates_list[[i]],inc_history,by=join_by(t==t))
    
    long_inc_sum <- inc_history%>%
      reshape2::melt(id=c('t','date'))%>%
      dplyr::rename(time=t)%>%
      group_by(time,date)%>%
      dplyr::summarise(inc.median=median(value)*30,
                       inc.mean=mean(value)*30,
                       inc.upper=quantile(value,probs=0.975)*30,
                       inc.lower=quantile(value,probs=0.025)*30)%>%
      mutate(month = zoo::as.Date(zoo::as.yearmon(date),frac=0.5))%>%
      mutate(district = districts[[i]])
    long_inc_sum <- left_join(dates_list[[i]],long_inc_sum,by=join_by(date==month))%>%
      mutate(month=as.yearmon(date))
    long_inc_sum$rel_inc = long_inc_sum$inc.median/max(long_inc_sum$inc.median)
    long_inc_sum$std_inc = (long_inc_sum$inc.median-mean(long_inc_sum$inc.median))/sd(long_inc_sum$inc.median)
    return(long_inc_sum)
  }))
  df_inc_sample <- bind_rows(lapply(1:4,function(i){
    inc_history <- data.frame(t(results[[i]]$history['clininc_05', 101:1000, -1]))
    inc_history <- inc_history[, sample_index]
    inc_history$t <- results[[i]]$times[-1]
    inc_history <- left_join(dates_list[[i]],inc_history,by=join_by(t==t))
    
    inc_sample <- inc_history %>%
      reshape2::melt(id=c('t','date'))%>%
      dplyr::rename(time=t)%>%
      mutate(month = zoo::as.Date(zoo::as.yearmon(date),frac=0.5))%>%
      mutate(district = districts[[i]],
             value = value*30)
    
    return(inc_sample)
  }))
  
  df_eir <- bind_rows(lapply(1:4,function(i){
    eir_history <- data.frame(t(results[[i]]$history['EIR', 101:1000, -1]))
    eir_history$t <- results[[i]]$times[-1]
    eir_history <- left_join(dates_list[[i]],eir_history,by=join_by(t==t))
    
    long_eir_sum <- eir_history%>%
      mutate(t=c(1:nrow(eir_history)))%>%
      reshape2::melt(id=c('t','date'))%>%
      dplyr::rename(time=t)%>%
      group_by(time,date)%>%
      dplyr::summarise(median=median(value),
                       mean=mean(value),
                       upper=quantile(value,probs=0.975),
                       lower=quantile(value,probs=0.025))%>%
      dplyr::mutate(district = districts[[i]])
    long_eir_sum$rel_eir = long_eir_sum$median/max(long_eir_sum$median)
    long_eir_sum$std_eir = (long_eir_sum$median-mean(long_eir_sum$median))/sd(long_eir_sum$median)
    return(long_eir_sum)
  }))
  df_eir_sample <- bind_rows(lapply(1:4,function(i){
    eir_history <- data.frame(t(results[[i]]$history['EIR', 101:1000, -1]))
    eir_history <- eir_history[, sample_index]
    eir_history$t <- results[[i]]$times[-1]
    eir_history <- left_join(dates_list[[i]],eir_history,by=join_by(t==t))
    eir_sample <- eir_history %>%
      reshape2::melt(id=c('t','date'))%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value,
                    district = districts[[i]])
    return(eir_sample)
  }))
  
  rainfall <- read_csv('./trial/Data/WA_ISTp_rainfall.csv') %>%
    group_by(Country) %>%
    filter(between(as.yearmon(Month), min(observed.pg$month), max(observed.pg$month)))%>%
    mutate(rel_rainfall = Rainfall/max(Rainfall),
           std_rainfall = (Rainfall - mean(Rainfall))/sd(Rainfall))%>%
    dplyr::rename(district = Country)
  
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  est_eir_plot <- ggplot()+
    # geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
    geom_line(data=df_eir_sample,aes(x=as.Date(date),y=value,color=district,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_eir,aes(x=as.Date(date),y=median,color=district,group=district),size=1)+
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
    facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.,scales = 'free_y')+
    labs(title = 'EIR - log scale')+
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
    geom_line(data=df_eir,aes(x=as.Date(date),y=rel_eir,color=district,group=district),size=1)+
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
    facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
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
  std_eir_plot <- ggplot()+
    geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=std_rainfall),color="black",size=1)+
    # geom_line(data=df_eir_sample,aes(x=as.Date(month),y=std_eir,color=site,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_eir,aes(x=as.Date(date),y=std_eir,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # coord_cartesian(ylim=c(0, 1000))+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.1,1,10,100,1000),labels=c(.1,1,10,100,1000))+
    # coord_cartesian(ylim=c(0, 1.2))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    labs(title = 'Standardized EIR and\nRainfall')+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  
  est_inc_plot <- ggplot()+
    geom_line(data=df_inc_sample,aes(x=date,y=value*100,color=district,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_inc,aes(x=date,y=inc.median*100,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    coord_cartesian(ylim=c(0, 150))+
    # coord_cartesian(ylim=c(0, ifelse(country=='NG',1000,max(incidence$inc)*2)))+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    labs(title = 'Clinical cases\nper month per 100 children under 5')+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  rel_est_inc_plot <- ggplot()+
    geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
    # geom_line(data=df_inc_sample,aes(x=month,y=value,color=district,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_inc,aes(x=date,y=rel_inc,color=district,group=district),size=1)+
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
    facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
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
  std_est_inc_plot <- ggplot()+
    geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=std_rainfall),color="black",size=1)+
    # geom_line(data=df_inc_sample,aes(x=month,y=value,color=district,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_inc,aes(x=date,y=std_inc,color=district,group=district),size=1)+
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
    facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    labs(title = 'Std. Est. Cases\nand Rainfall')+
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
  obs_prev_plot_mg <- ggplot()+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=observed.mg,aes(x=date,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=observed.mg,aes(x=date,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    # geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),pch= 2, size=1.5,color="black")+
    # geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color='black')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    scale_y_continuous(limits = c(0,1))+
    coord_cartesian(ylim=c(0, 1),xlim = range(observed.mg$date))+
    labs(title = 'ANC Prevalence\nSecundigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  est_prev_plot_mg <- ggplot()+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_line(data=df_prev_sample,aes(x=date,y=prev_mg,color=district,group=variable),alpha=0.2)+
    geom_line(data=df_prev,aes(x=date,y=median_mg,color=district,group=district),linewidth=1)+
    geom_point(data=observed.mg,aes(x=date,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=observed.mg,aes(x=date,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    # geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),pch= 2, size=1.5,color="black")+
    # geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color='black')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    scale_y_continuous(limits = c(0,1))+
    coord_cartesian(ylim=c(0, 1),xlim = range(observed.mg$date))+
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
  # # obs_prev_plot_all <- ggplot()+
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
  est_prev_plot_pg <- ggplot()+
    geom_line(data=df_prev_sample,aes(x=date,y=prev_pg,color=district,group=variable),alpha=0.2)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_prev,aes(x=date,y=median_pg,color=district,group=district),linewidth=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=observed.pg,aes(x=date,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=observed.pg,aes(x=date,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    # geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),size=1.5,pch=2,color="black")+
    # geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color="black")+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits=c(0,1))+
    coord_cartesian(xlim = range(observed.pg$date))+
    facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    labs(title = 'ANC Prevalence\nPrimigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  obs_prev_plot_pg <- ggplot()+
    geom_point(data=observed.pg,aes(x=date,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=observed.pg,aes(x=date,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    # geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),size=1.5,pch=2,color="black")+
    # geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color="black")+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits=c(0,1))+
    coord_cartesian(xlim = range(observed.pg$date))+
    facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
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
  # obs_inc_plot <- ggplot()+
  #   # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
  #   geom_line(data=incidence,aes(x=date_ex,y=inc,color=district,group=district),size=1)+
  #   # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
  #   scale_color_manual(values=colors)+
  #   scale_fill_manual(values=colors)+
  #   scale_x_date(date_labels = "%b %Y")+
  #   # scale_y_continuous(limits=c(0,500))+
  #   # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  #   # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
  #   facet_grid(district~.,labeller = labeller(district = district_labels))+
  #   labs(title='Observed Incidence\nper 10,000 person-months')+
  #   coord_cartesian(ylim=c(0, ifelse(country=='NG',1000,max(incidence$inc)*2)))+
  #   # coord_cartesian(ylim=c(0, max(incidence$inc)*2))+
  #   theme(legend.title = element_blank(),
  #         legend.position = 'none',
  #         axis.title.x = element_blank(),
  #         axis.title.y = element_blank(),
  #         axis.text.x=element_text(angle=45, hjust=1, vjust=1),
  #         axis.ticks.x = element_line(size = 0.5),
  #         axis.ticks.length = unit(3, "pt"))
  # obs_rel_inc_plot <- ggplot()+
  #   geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
  #   # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
  #   geom_line(data=incidence,aes(x=date_ex,y=rel_inc,color=district,group=district),size=1)+
  #   # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
  #   scale_color_manual(values=colors)+
  #   scale_fill_manual(values=colors)+
  #   scale_x_date(date_labels = "%b %Y")+
  #   # scale_y_continuous(limits=c(0,500))+
  #   # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  #   # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
  #   facet_grid(district~.,labeller = labeller(district = district_labels))+
  #   labs(title='Rel. Obs. Cases\nand Rainfall')+
  #   # coord_cartesian(ylim=c(0, ifelse(country=='NG',1000,max(incidence$inc)*2)))+
  #   # coord_cartesian(ylim=c(0, max(incidence$inc)*2))+
  #   theme(legend.title = element_blank(),
  #         legend.position = 'none',
  #         axis.title.x = element_blank(),
  #         axis.title.y = element_blank(),
  #         axis.text.y = element_blank(),
  #         axis.text.x=element_text(angle=45, hjust=1, vjust=1),
  #         axis.ticks.x = element_line(size = 0.5),
  #         axis.ticks.length = unit(3, "pt"))
  
  # print('obs_inc_plot')
  # obs_pos_plot <- ggplot()+
  #   # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
  #   geom_line(data=incidence,aes(x=date_ex,y=prop_pos,color=district,group=district),size=1)+
  #   # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
  #   scale_color_manual(values=colors)+
  #   scale_fill_manual(values=colors)+
  #   scale_x_date(date_labels = "%b %Y")+
  #   scale_y_continuous(limits=c(0,1))+
  #   # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  #   # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
  #   # labs(x='Date',y='EIR')+
  #   labs(title = 'Test positivity\nproportion')+
  #   facet_grid(district~.,labeller = labeller(district = district_labels))+
  #   theme(legend.title = element_blank(),
  #         legend.position = 'none',
  #         axis.title.x = element_blank(),
  #         axis.title.y = element_blank(),
  #         axis.text.x=element_text(angle=45, hjust=1, vjust=1),
  #         axis.ticks.x = element_line(size = 0.5),
  #         axis.ticks.length = unit(3, "pt"))
  sample_size_plot <- ggplot()+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=site,group=site),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=site,group=site),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_line(data=observed.all,aes(x=date,y=tested,color=district,linetype=factor(grav,levels=c('pg','mg'))))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
    labs(title = 'Number Tested')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  volatility.density <- ggplot(mcmc.df)+
    geom_violin(aes(x=sites,y=volatility),fill='#6D6A67')+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(title='volatility density by Site',x='Site',y='volatility')
  
  # log_init_EIR.density <- ggplot(mcmc.df)+
  #   geom_violin(aes(x=sites,y=log_init_EIR),fill='#6D6A67')+
  #   theme(axis.text.x = element_text(angle = 45, hjust=1))+
  #   labs(title='Log init_EIR density by Site',x='Site',y='Log init_EIR')
  
  # init_EIR.density <- ggplot(mcmc.df)+
  #   geom_violin(aes(x=sites,y=exp(log_init_EIR)),fill='#6D6A67')+
  #   scale_y_log10()+
  #   theme(axis.text.x = element_text(angle = 45, hjust=1))+
  #   labs(title='init_EIR density by Site',x='Site',y='init_EIR')
  
  rel_dash <- sample_size_plot+est_prev_plot_pg+est_prev_plot_mg+est_eir_plot+rel_est_inc_plot+plot_layout(guides = "collect",ncol=5)
  std_dash <- sample_size_plot+est_prev_plot_pg+est_prev_plot_mg+est_eir_plot+std_est_inc_plot+plot_layout(guides = "collect",ncol=5)
  trans_dash <- est_eir_plot+rel_est_inc_plot+plot_layout(guides = "collect",ncol=2)
  
  overlay_plot <- NULL
  
  if(param == 'betaa'){
    df_betaa <- bind_rows(lapply(1:4,function(i){
      betaa_history <- data.frame(t(results[[i]]$history['betaa', 101:1000, -1]))
      betaa_history$t <- results[[i]]$times[-1]
      betaa_history <- left_join(dates_list[[i]],betaa_history,by=join_by(t==t))
      
      long_betaa_sum <- betaa_history%>%
        reshape2::melt(id=c('t','date'))%>%
        dplyr::rename(time=t)%>%
        group_by(time,date)%>%
        dplyr::summarise(median=median(value),
                         mean=mean(value),
                         upper=quantile(value,probs=0.975),
                         lower=quantile(value,probs=0.025))%>%
        dplyr::mutate(district = districts[[i]],
                      month = as.yearmon(date))
      long_betaa_sum$rel_betaa = long_betaa_sum$median/max(long_betaa_sum$median)
      long_betaa_sum$std_betaa = (long_betaa_sum$median-mean(long_betaa_sum$median))/sd(long_betaa_sum$median)
      
      daily_dates <- data.frame(date = seq(floor_date(min(long_betaa_sum$date), unit = "month"),ceiling_date(max(long_betaa_sum$date), unit = "month")-1,by='day'))
      daily_dates$month <- as.yearmon(daily_dates$date)
      long_betaa_expand <- left_join(daily_dates,long_betaa_sum,by='month')
      
      return(long_betaa_expand)
    }))
    df_betaa_sample <- bind_rows(lapply(1:4,function(i){
      betaa_history <- data.frame(t(results[[i]]$history['betaa', 101:1000, -1]))
      betaa_history <- betaa_history[, sample_index]
      betaa_history$t <- results[[i]]$times[-1]
      betaa_history <- left_join(dates_list[[i]],betaa_history,by=join_by(t==t))
      
      betaa_sample <- betaa_history %>%
        reshape2::melt(id=c('t','date'))%>%
        dplyr::rename(time=t)%>%
        dplyr::mutate(value = value,
                      district = districts[[i]],
                      month = as.yearmon(date))
      daily_dates <- data.frame(date = seq(floor_date(min(betaa_sample$date), unit = "month"),ceiling_date(max(betaa_sample$date), unit = "month")-1,by='day'))
      daily_dates$month <- as.yearmon(daily_dates$date)
      betaa_sample_expand <- left_join(daily_dates,betaa_sample,by='month')
      
      return(betaa_sample_expand)
    }))
    est_betaa_plot <- ggplot()+
      # geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
      geom_line(data=df_betaa_sample,aes(x=as.Date(date.x),y=value,color=district,group=variable),alpha=0.1)+
      # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
      geom_line(data=df_betaa,aes(x=as.Date(date.x),y=median,color=district,group=district),size=1)+
      # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
      scale_color_manual(values=colors)+
      scale_fill_manual(values=colors)+
      scale_x_date(date_labels = "%b %Y")+
      coord_cartesian(ylim=c(0, 80))+
      # scale_y_continuous(limits=c(0,500))+
      # scale_y_log10(breaks=c(.01,.1,1,10,100,1000),labels=c(.01,.1,1,10,100,1000))+
      # coord_cartesian(ylim=c(.01, 800))+
      # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
      # labs(x='Date',y='EIR')+
      facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
      labs(title = 'Mosquito Emergence')+
      theme(legend.title = element_blank(),
            legend.position = 'none',
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1),
            axis.ticks.x = element_line(size = 0.5),
            axis.ticks.length = unit(3, "pt")
      )
    std_betaa_plot <- ggplot()+
      geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=std_rainfall),color="black",size=1)+
      # geom_line(data=df_eir_sample,aes(x=as.Date(month),y=std_eir,color=site,group=variable),alpha=0.1)+
      # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
      geom_line(data=df_betaa,aes(x=as.Date(date.x),y=std_betaa,color=district,group=district),size=1)+
      # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
      scale_color_manual(values=colors)+
      scale_fill_manual(values=colors)+
      scale_x_date(date_labels = "%b %Y")+
      # coord_cartesian(ylim=c(0, 1000))+
      # scale_y_continuous(limits=c(0,500))+
      # scale_y_log10(breaks=c(.1,1,10,100,1000),labels=c(.1,1,10,100,1000))+
      # coord_cartesian(ylim=c(0, 1.2))+
      # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
      # labs(x='Date',y='EIR')+
      facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
      labs(title = 'Std. Mosquito Emergence and\nRainfall')+
      theme(legend.title = element_blank(),
            legend.position = 'none',
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1),
            axis.ticks.x = element_line(size = 0.5),
            axis.ticks.length = unit(3, "pt")
      )
    
    indicator_palette <- c(`Rainfall` = 'black',
                           `PG Prevalence` = "#999999",
                           `MG Prevalence` = "#999999",
                           `Incidence` = "#1B9E77",
                           `EIR` = "#D95F02",
                           `Mosquito Emergence` = "#377EB8")
    
    overlay_plot <- ggplot()+
      geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=std_rainfall,color='Rainfall'),size=1)+
      geom_line(data=df_prev,aes(x=date,y=std_prev_pg,group=district,color='PG Prevalence'),linetype='solid',linewidth=1)+
      geom_line(data=df_prev,aes(x=date,y=std_prev_mg,group=district,color='MG Prevalence'),linetype='dashed',linewidth=1)+
      geom_line(data=df_inc,aes(x=date,y=std_inc,group=district,color='Incidence'),size=1)+
      geom_line(data=df_eir,aes(x=date,y=std_eir,group=district,color='EIR'),size=1)+
      geom_line(data=df_betaa,aes(x=as.Date(date.x),y=std_betaa,group=district,color='Mosquito Emergence'),size=1)+
      scale_color_manual(values=indicator_palette)+
      scale_x_date(date_labels = "%b %Y")+
      # coord_cartesian(ylim=c(0, 1000))+
      # scale_y_continuous(limits=c(0,500))+
      # scale_y_log10(breaks=c(.1,1,10,100,1000),labels=c(.1,1,10,100,1000))+
      # coord_cartesian(ylim=c(0, 1))+
      # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
      # labs(x='Date',y='EIR')+
      facet_grid(factor(district,levels=c('Ghana', 'Burkina Faso', 'Mali', 'Gambia'))~.)+
      labs(title = 'Indicator Synchrony')+
      theme(legend.title = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1),
            axis.ticks.x = element_line(size = 0.5),
            axis.ticks.length = unit(3, "pt")
      )
    rel_dash <- est_prev_plot_pg+est_prev_plot_mg+est_eir_plot+est_betaa_plot+rel_est_inc_plot+plot_layout(guides = "collect",ncol=5)
    trans_dash <- est_eir_plot+est_betaa_plot+rel_est_inc_plot+plot_layout(guides = "collect",ncol=3)
    std_dash <- est_prev_plot_pg+est_prev_plot_mg+est_eir_plot+est_betaa_plot+std_est_inc_plot+plot_layout(guides = "collect",ncol=5)
  }
  # print('obs_pos_plot')
  # # obs_test_plot <- ggplot()+
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
  full_dash <- est_prev_plot_pg+est_prev_plot_mg+est_inc_plot+ plot_layout(guides = "collect",ncol=3)
  obs_data_dash <- obs_prev_plot_pg+obs_prev_plot_mg+est_inc_plot + plot_layout(guides = "collect",ncol=3)
  pres_dash <- est_prev_plot_pg+est_prev_plot_mg+est_inc_plot+ plot_layout(guides = "collect",ncol=3)
  # mcmc_dash <- volatility.density+init_EIR.density+plot_layout(ncol=1)
  
  return(list(full_dash = full_dash,
              obs_data_dash = obs_data_dash,
              pres_dash = pres_dash,
              rel_dash = rel_dash,
              std_dash = std_dash,
              # mcmc_dash = mcmc_dash,
              overlay_plot = overlay_plot,
              trans_dash = trans_dash))
  # obs_prev_plot_mg
}
