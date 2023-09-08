create_dashboard_plots_trial_ea <- function(results,
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
  districts <- c("Kenya","Malawi" )
  district_labels <- c("Kenya","Malawi")
  colors <- c(Kenya = "#D95F02", Malawi = "#377EB8")
  
  
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
  
  
  mcmc.df <- bind_rows(lapply(1:2, 
                              function(x){
                                df <- results[[x]]$mcmc[101:1000,]
                                df$sites <- names(results[x])
                                df$step <- 1:nrow(df)
                                return(df)
                              }))
  
  df_prev <- bind_rows(lapply(1:2,function(i){
    prev_history <- data.frame(t(results[[i]]$history['prev_05', 51:1000, -1]))
    prev_history$t <- c((start_times[,i]$t+1):(start_times[,i]$t+nrow(prev_history)))
    prev_history <- left_join(dates_list[[i]],prev_history,by=join_by(t==t))
    long_prev_sum <- prev_history%>%
      melt(id=c('t','date'))%>%
      dplyr::rename(time=t)%>%
      mutate(logodds_child = log(get_odds_from_prev(value)),
             prev_pg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
             }))),
             prev_mg = rowMedians(as.matrix(sapply(1:nrow(coefs_mg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_mg_df$gradient[x]*(logodds_child-coefs_mg_df$av_lo_child[x])+coefs_mg_df$intercept[x])
             }))))%>%
      group_by(time,date)%>%
      dplyr::summarise(median_pg=median(prev_pg),
                       mean_pg=mean(prev_pg),
                       upper_pg=quantile(prev_pg,probs=0.975),
                       lower_pg=quantile(prev_pg,probs=0.025),
                       median_mg=median(prev_mg),
                       mean_mg=mean(prev_mg),
                       upper_mg=quantile(prev_mg,probs=0.975),
                       lower_mg=quantile(prev_mg,probs=0.025))%>%
      mutate(district = districts[[i]])
  }))
  df_prev_sample <- bind_rows(lapply(1:2,function(i){
    prev_history <- data.frame(t(results[[i]]$history['prev_05', 51:1000, -1]))
    prev_history <- prev_history[, sample(ncol(prev_history), 100)]
    prev_history$t <- c((start_times[,i]$t+1):(start_times[,i]$t+nrow(prev_history)))
    prev_history <- left_join(dates_list[[i]],prev_history,by=join_by(t==t))
    prev_sample <- prev_history%>%
      melt(id=c('t','date'))%>%
      dplyr::rename(time=t)%>%
      mutate(logodds_child = log(get_odds_from_prev(value)),
             prev_pg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
             }))),
             prev_mg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_mg_df$gradient[x]*(logodds_child-coefs_mg_df$av_lo_child[x])+coefs_mg_df$intercept[x])
             }))))%>%
      mutate(district = districts[[i]])
  }))
  
  observed.pg <- bind_rows(lapply(1:2,
                                  function(x){
                                    df <- prev_pg[[x]]
                                    df$district <- names(prev_pg[x])
                                    return(df)
                                  }))
  observed.pg <- addCIs(observed.pg,observed.pg$positive,observed.pg$tested)
  observed.pg$date <- as.Date(observed.pg$month,frac = 0.5)

  observed.mg <- bind_rows(lapply(1:2,
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
  
  df_inc <- bind_rows(lapply(1:2,function(i){
    inc_history <- data.frame(t(results[[i]]$history['clininc_all', 51:1000, -1]))
    
    long_inc_sum <- inc_history%>%
      dplyr::mutate(t=c((start_times[,i]$t+1):(start_times[,i]$t+nrow(inc_history))))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(inc.median=median(value),
                       inc.mean=mean(value),
                       inc.upper=quantile(value,probs=0.975),
                       inc.lower=quantile(value,probs=0.025))%>%
      mutate(date = as.Date(c((start_times[,i]$date+1):(start_times[,i]$date + nrow(inc_history)))))%>%
      mutate(month = zoo::as.Date(zoo::as.yearmon(date),frac=0.5))%>%
      ungroup()%>%
      group_by(month)%>%
      dplyr::summarise(inc.median=sum(inc.median),
                inc.mean=sum(inc.mean),
                inc.upper=sum(inc.upper),
                inc.lower=sum(inc.lower))%>%
      mutate(district = districts[[i]])
    long_inc_sum <- left_join(dates_list[[i]],long_inc_sum,by=join_by(date==month))%>%
      mutate(month=as.yearmon(date),
             rel_inc = inc.median/max(inc.median))
  }))
  df_inc_sample <- bind_rows(lapply(1:2,function(i){
    inc_history <- data.frame(t(results[[i]]$history['clininc_all', 51:1000, -1]))
    
    inc_sample <- inc_history[, sample(ncol(inc_history), 100)] %>%
      dplyr::mutate(t=c((start_times[,i]$t+1):(start_times[,i]$t+nrow(inc_history))))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      mutate(date = rep(as.Date(c((start_times[,i]$date+1):(start_times[,i]$date + nrow(inc_history)))),100))%>%
      mutate(month = zoo::as.Date(zoo::as.yearmon(date),frac=0.5))%>%
      mutate(district = districts[[i]]) %>%
      group_by(district,month,variable)%>%
      dplyr::summarise(value=sum(value))
    
    inc_sample <- left_join(dates_list[[i]],inc_sample,by=join_by(date==month))%>%
      mutate(month=as.yearmon(date))
  }))
  
  df_eir <- bind_rows(lapply(1:2,function(i){
    eir_history <- data.frame(t(results[[i]]$history['EIR', 51:1000, -1]))

    long_eir_sum <- eir_history%>%
      mutate(t=c(1:nrow(eir_history)))%>%
      melt(id=c('t'))%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value),
                       mean=mean(value),
                       upper=quantile(value,probs=0.975),
                       lower=quantile(value,probs=0.025))%>%
      dplyr::mutate(district = districts[[i]],
                    date = start_times[['date',i]]+time-1,
                    rel_eir = median/max(median))
  }))
  df_eir_sample <- bind_rows(lapply(1:2,function(i){
    eir_history <- data.frame(t(results[[i]]$history['EIR', 51:1000, -1]))
    eir_sample <- eir_history[, sample(ncol(eir_history), 100)] %>%
      mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value,
                    district = districts[[i]],
                    date = start_times[['date',i]]+time-1)
  }))
  
  rainfall <- read_csv('./trial/Data/EA_ISTp_rainfall.csv') %>%
    group_by(Country) %>%
    filter(between(as.yearmon(Month), min(observed.pg$month), max(observed.pg$month)))%>%
    mutate(rel_rainfall = Rainfall/max(Rainfall))%>%
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
    facet_grid(factor(district,levels=c("Kenya","Malawi" ))~.,scales = 'free_y')+
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
    facet_grid(factor(district,levels=c("Kenya","Malawi" ))~.)+
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
  
  est_inc_plot <- ggplot()+
    geom_line(data=df_inc_sample,aes(x=date,y=value*10000,color=district,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_inc,aes(x=date,y=inc.median*10000,color=district,group=district),size=1)+
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
    facet_grid(factor(district,levels=c("Kenya","Malawi" ))~.)+
    labs(title = 'Estimated Incidence\nper 10,000 person-months')+
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
    facet_grid(factor(district,levels=c("Kenya","Malawi" ))~.)+
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
    facet_grid(factor(district,levels=c("Kenya","Malawi" ))~.)+
    scale_y_continuous(limits = c(0,1))+
    coord_cartesian(ylim=c(0, 1),xlim = range(observed.mg$date))+
    labs(title = 'ANC Prevalence\nMultigrav')+
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
    facet_grid(factor(district,levels=c("Kenya","Malawi" ))~.)+
    scale_y_continuous(limits = c(0,1))+
    coord_cartesian(ylim=c(0, 1),xlim = range(observed.mg$date))+
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
    facet_grid(factor(district,levels=c("Kenya","Malawi" ))~.)+
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
    facet_grid(factor(district,levels=c("Kenya","Malawi" ))~.)+
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
    facet_grid(factor(district,levels=c("Kenya","Malawi" ))~.)+
    labs(title = 'Number Tested')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
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
  
  rel_dash <- sample_size_plot+est_prev_plot_pg+est_prev_plot_mg+est_eir_plot+rel_est_inc_plot+plot_layout(guides = "collect",ncol=5)
  if(param == 'betaa'){
    df_betaa <- bind_rows(lapply(1:2,function(i){
      betaa_history <- data.frame(t(results[[i]]$history['betaa', 51:1000, -1]))
      
      long_betaa_sum <- betaa_history%>%
        mutate(t=c(1:nrow(betaa_history)))%>%
        melt(id=c('t'))%>%
        dplyr::rename(time=t)%>%
        group_by(time)%>%
        dplyr::summarise(median=median(value),
                         mean=mean(value),
                         upper=quantile(value,probs=0.975),
                         lower=quantile(value,probs=0.025))%>%
        dplyr::mutate(district = districts[[i]],
                      date = start_times[['date',i]]+time-1,
                      rel_betaa = median/max(median))
    }))
    df_betaa_sample <- bind_rows(lapply(1:2,function(i){
      betaa_history <- data.frame(t(results[[i]]$history['betaa', 51:1000, -1]))
      betaa_sample <- betaa_history[, sample(ncol(betaa_history), 100)] %>%
        mutate(t=c(1:nrow(betaa_history)))%>%
        melt(id='t')%>%
        dplyr::rename(time=t)%>%
        dplyr::mutate(value = value,
                      district = districts[[i]],
                      date = start_times[['date',i]]+time-1)
    }))
    est_betaa_plot <- ggplot()+
      # geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
      geom_line(data=df_betaa_sample,aes(x=as.Date(date),y=value,color=district,group=variable),alpha=0.1)+
      # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
      geom_line(data=df_betaa,aes(x=as.Date(date),y=median,color=district,group=district),size=1)+
      # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
      scale_color_manual(values=colors)+
      scale_fill_manual(values=colors)+
      scale_x_date(date_labels = "%b %Y")+
      # coord_cartesian(ylim=c(0, 1000))+
      # scale_y_continuous(limits=c(0,500))+
      # scale_y_log10(breaks=c(.01,.1,1,10,100,1000),labels=c(.01,.1,1,10,100,1000))+
      # coord_cartesian(ylim=c(.01, 800))+
      # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
      # labs(x='Date',y='EIR')+
      facet_grid(factor(district,levels=c("Kenya","Malawi" ))~.)+
      labs(title = 'Mosquito Emergence')+
      theme(legend.title = element_blank(),
            legend.position = 'none',
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1),
            axis.ticks.x = element_line(size = 0.5),
            axis.ticks.length = unit(3, "pt")
      )
    rel_betaa_plot <- ggplot()+
      geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
      # geom_line(data=df_eir_sample,aes(x=as.Date(month),y=rel_eir,color=site,group=variable),alpha=0.1)+
      # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
      geom_line(data=df_betaa,aes(x=as.Date(date),y=rel_betaa,color=district,group=district),size=1)+
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
      facet_grid(factor(district,levels=c("Kenya","Malawi" ))~.)+
      labs(title = 'Relative Mosquito Emergence and\nRainfall')+
      theme(legend.title = element_blank(),
            legend.position = 'none',
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1),
            axis.ticks.x = element_line(size = 0.5),
            axis.ticks.length = unit(3, "pt")
      )
    rel_dash <- est_prev_plot_pg+est_prev_plot_mg+est_eir_plot+est_betaa_plot+rel_est_inc_plot+plot_layout(guides = "collect",ncol=5)
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
  mcmc_dash <- EIR_SD.density+init_EIR.density+plot_layout(ncol=1)
  
  return(list(full_dash = full_dash,
              obs_data_dash = obs_data_dash,
              pres_dash = pres_dash,
              rel_dash = rel_dash,
              mcmc_dash = mcmc_dash))
  # obs_prev_plot_mg
}
