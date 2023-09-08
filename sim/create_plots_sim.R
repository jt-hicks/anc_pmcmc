create_dashboard_plots_sim <- function(results,
                                       prev_all=NULL,
                                       burnin = 50,
                                       chain_length = 1000,
                                       sample_size = 100,
                                       start_pf_time = 30,
                                       param = 'betaa'){

   
  start_obs <- min(zoo::as.Date(zoo::as.yearmon(prev_all$month)))#Month of first observation (in Date format)
  time_origin <- as.Date(paste0(year(start_obs)-1,'-01-01')) #January 1 of year before observation (in Date format)
  date <- zoo::as.Date(zoo::as.yearmon(prev_all$month), frac = 0.5) #Convert dates to middle of month
  t <- as.integer(difftime(date,time_origin,units="days"))
  
  dates_list <- data.frame(date=date,t=t)
  
  start_times <- list(t=min(dates_list$t)-start_pf_time,date=min(dates_list$date)-start_pf_time)

  
  mcmc.df <- results$mcmc[(burnin+1):chain_length,]
  mcmc.df$step <- 1:nrow(mcmc.df)

    prev_history <- data.frame(t(results$history['prev_05', (burnin+1):chain_length, -1]))
    prev_history$t <- c((start_times$t+1):(start_times$t+nrow(prev_history)))
    prev_history <- left_join(dates_list,prev_history,by=join_by(t==t))
    df_prev <- prev_history%>%
      melt(id=c('t','date'))%>%
      dplyr::rename(time=t)%>%
      group_by(time,date)%>%
      dplyr::summarise(median=median(value),
                       mean=mean(value),
                       upper=quantile(value,probs=0.975),
                       lower=quantile(value,probs=0.025))
  
    prev_history <- data.frame(t(results$history['prev_05', (burnin+1):chain_length, -1]))
    prev_history <- prev_history[, sample(ncol(prev_history), sample_size)]
    prev_history$t <- c((start_times$t+1):(start_times$t+nrow(prev_history)))
    prev_history <- left_join(dates_list,prev_history,by=join_by(t==t))
    df_prev_sample <- prev_history%>%
      melt(id=c('t','date'))%>%
      dplyr::rename(time=t)

  observed.prev <- prev_all
  observed.prev <- addCIs(observed.prev,observed.prev$positive,observed.prev$tested)
  observed.prev$date <- as.Date(observed.prev$month,frac = 0.5)


  inc_history <- data.frame(t(results$history['clininc_all', (burnin+1):chain_length, -1]))
  long_inc_sum <- inc_history%>%
    dplyr::mutate(t=c((start_times$t+1):(start_times$t+nrow(inc_history))))%>%
    melt(id='t')%>%
    dplyr::rename(time=t)%>%
    group_by(time)%>%
    dplyr::summarise(inc.median=median(value),
                     inc.mean=mean(value),
                     inc.upper=quantile(value,probs=0.975),
                     inc.lower=quantile(value,probs=0.025))%>%
    mutate(date = as.Date(c((start_times$date+1):(start_times$date + nrow(inc_history)))))%>%
    mutate(month = zoo::as.Date(zoo::as.yearmon(date),frac=0.5))%>%
    ungroup()%>%
    group_by(month)%>%
    dplyr::summarise(inc.median=sum(inc.median),
                     inc.mean=sum(inc.mean),
                     inc.upper=sum(inc.upper),
                     inc.lower=sum(inc.lower))
  df_inc <- left_join(dates_list,long_inc_sum,by=join_by(date==month))%>%
    mutate(month=as.yearmon(date),
           rel_inc = inc.median/max(inc.median))


    inc_sample <- inc_history[, sample(ncol(inc_history), sample_size)] %>%
      dplyr::mutate(t=c((start_times$t+1):(start_times$t+nrow(inc_history))))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      mutate(date = rep(as.Date(c((start_times$date+1):(start_times$date + nrow(inc_history)))),sample_size))%>%
      mutate(month = zoo::as.Date(zoo::as.yearmon(date),frac=0.5))%>%
      group_by(month,variable)%>%
      dplyr::summarise(value=sum(value))
    
    df_inc_sample <- left_join(dates_list,inc_sample,by=join_by(date==month))%>%
      mutate(month=as.yearmon(date))

    eir_history <- data.frame(t(results$history['EIR', (burnin+1):chain_length, -1]))
    df_eir <- eir_history%>%
      mutate(t=c(1:nrow(eir_history)))%>%
      melt(id=c('t'))%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value),
                       mean=mean(value),
                       upper=quantile(value,probs=0.975),
                       lower=quantile(value,probs=0.025))%>%
      dplyr::mutate(date = start_times$date+time-1,
                    rel_eir = median/max(median))
    df_eir_sample <- eir_history[, sample(ncol(eir_history), sample_size)] %>%
      mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value,
                    date = start_times$date+time-1)

  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  est_eir_plot <- ggplot()+
    geom_line(data=df_eir_sample,aes(x=as.Date(date),y=value,group=variable),alpha=0.1)+
    geom_line(data=df_eir,aes(x=as.Date(date),y=median),size=1)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_log10(breaks=c(.01,.1,1,10,100,1000),labels=c(.01,.1,1,10,100,1000))+
    coord_cartesian(ylim=c(.01, 800))+
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
    geom_line(data=df_eir,aes(x=as.Date(date),y=rel_eir),size=1)+
    scale_x_date(date_labels = "%b %Y")+
    coord_cartesian(ylim=c(0, 1.2))+
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
    geom_line(data=df_inc_sample,aes(x=date,y=value*10000,group=variable),alpha=0.1)+
    geom_line(data=df_inc,aes(x=date,y=inc.median*10000),size=1)+
    scale_x_date(date_labels = "%b %Y")+
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
    geom_line(data=df_inc,aes(x=date,y=rel_inc),size=1)+
    scale_x_date(date_labels = "%b %Y")+
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
  obs_prev_plot <- ggplot()+
    geom_point(data=observed.prev,aes(x=date,y=mean),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=observed.prev,aes(x=date,ymin=lower,ymax=upper),width = 0,position=position_dodge(width=10))+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits = c(0,1))+
    coord_cartesian(ylim=c(0, 1),xlim = range(observed.prev$date))+
    labs(title = 'ANC Prevalence\nMultigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  est_prev_plot <- ggplot()+
    geom_line(data=df_prev_sample,aes(x=date,y=value,group=variable),alpha=0.2)+
    geom_line(data=df_prev,aes(x=date,y=median),linewidth=1)+
    geom_point(data=observed.prev,aes(x=date,y=mean),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=observed.prev,aes(x=date,ymin=lower,ymax=upper),width = 0,position=position_dodge(width=10))+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits = c(0,1))+
    coord_cartesian(ylim=c(0, 1),xlim = range(observed.prev$date))+
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
  # #   # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district),alpha=0.2)+
  # #   # geom_line(aes(x=month,y=median),size=1)+
  # #   # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # #   geom_point(data=df_data_all,aes(x=month,y=mean),pch = 19,position=position_dodge(width=10))+
  # #   geom_errorbar(data=df_data_all,aes(x=month,ymin=lower,ymax=upper),width = 0,position=position_dodge(width=10))+
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
    geom_line(data=df_prev_sample,aes(x=date,y=prev_pg,group=variable),alpha=0.2)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district),alpha=0.2)+
    geom_line(data=df_prev,aes(x=date,y=median_pg),linewidth=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=observed.pg,aes(x=date,y=mean),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=observed.pg,aes(x=date,ymin=lower,ymax=upper),width = 0,position=position_dodge(width=10))+
    geom_point(data=cs_data,aes(x=month,y=mean),size=1.5,pch=2,color="black")+
    geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="black")+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits=c(0,1))+
    coord_cartesian(xlim = range(observed.pg$date))+
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
  obs_prev_plot_pg <- ggplot()+
    geom_point(data=observed.pg,aes(x=date,y=mean),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=observed.pg,aes(x=date,ymin=lower,ymax=upper),width = 0,position=position_dodge(width=10))+
    geom_point(data=cs_data,aes(x=month,y=mean),size=1.5,pch=2,color="black")+
    geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="black")+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    scale_y_continuous(limits=c(0,1))+
    coord_cartesian(xlim = range(observed.pg$date))+
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
  obs_inc_plot <- ggplot()+
    # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
    geom_line(data=incidence,aes(x=date_ex,y=inc),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median),size=1,linetype='dashed')+
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
  obs_rel_inc_plot <- ggplot()+
    geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
    # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
    geom_line(data=incidence,aes(x=date_ex,y=rel_inc),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median),size=1,linetype='dashed')+
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
  obs_pos_plot <- ggplot()+
    # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
    geom_line(data=incidence,aes(x=date_ex,y=prop_pos),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median),size=1,linetype='dashed')+
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
  sample_size_plot <- ggplot()+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=site,group=site),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=site,group=site),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_line(data=observed.all,aes(x=date,y=tested,linetype=factor(grav,levels=c('pg','mg'))))+
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
  volatility.density <- ggplot(mcmc.df)+
    geom_violin(aes(x=sites,y=volatility),fill='#6D6A67')+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(title='volatility density by Site',x='Site',y='volatility')
  
  log_init_EIR.density <- ggplot(mcmc.df)+
    geom_violin(aes(x=sites,y=log_init_EIR),fill='#6D6A67')+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(title='Log init_EIR density by Site',x='Site',y='Log init_EIR')
  
  init_EIR.density <- ggplot(mcmc.df)+
    geom_violin(aes(x=sites,y=exp(log_init_EIR)),fill='#6D6A67')+
    scale_y_log10()+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(title='init_EIR density by Site',x='Site',y='init_EIR')
  
  rel_dash <- sample_size_plot+est_prev_plot_pg+est_prev_plot_mg+est_eir_plot+rel_est_inc_plot+obs_rel_inc_plot+plot_layout(guides = "collect",ncol=6)
  full_dash <- est_prev_plot_pg+est_prev_plot_mg+est_eir_plot+est_inc_plot+obs_inc_plot+obs_pos_plot+ plot_layout(guides = "collect",ncol=6)
  overlay_plot <- NULL
  
  if(param == 'betaa'){
    df_betaa <- bind_rows(lapply(start:(start+number-1),function(i){
      betaa_history <- data.frame(t(results$history['betaa', (burnin+1):chain_length, -1]))
      
      long_betaa_sum <- betaa_history%>%
        mutate(t=c(1:nrow(betaa_history)))%>%
        melt(id=c('t'))%>%
        dplyr::rename(time=t)%>%
        group_by(time)%>%
        dplyr::summarise(median=median(value),
                         mean=mean(value),
                         upper=quantile(value,probs=0.975),
                         lower=quantile(value,probs=0.025))%>%
        dplyr::mutate(district = districts[[i-start+1]],
                      date = start_times[['date',i]]+time-1,
                      rel_betaa = median/max(median))
    }))
    df_betaa_sample <- bind_rows(lapply(start:(start+number-1),function(i){
      betaa_history <- data.frame(t(results$history['betaa', (burnin+1):chain_length, -1]))
      betaa_sample <- betaa_history[, sample(ncol(betaa_history), sample_size)] %>%
        mutate(t=c(1:nrow(betaa_history)))%>%
        melt(id='t')%>%
        dplyr::rename(time=t)%>%
        dplyr::mutate(value = value,
                      district = districts[[i-start+1]],
                      date = start_times[['date',i]]+time-1)
    }))
    est_betaa_plot <- ggplot()+
      # geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall),color="black",size=1)+
      geom_line(data=df_betaa_sample,aes(x=as.Date(date),y=value,group=variable),alpha=0.1)+
      # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district),alpha=0.2)+
      geom_line(data=df_betaa,aes(x=as.Date(date),y=median),size=1)+
      # geom_line(data=df_eir,aes(x=month,y=median),size=1,linetype='dashed')+
      scale_color_manual(values=colors)+
      scale_fill_manual(values=colors)+
      scale_x_date(date_labels = "%b %Y")+
      # coord_cartesian(ylim=c(0, 1000))+
      # scale_y_continuous(limits=c(0,500))+
      # scale_y_log10(breaks=c(.01,.1,1,10,100,1000),labels=c(.01,.1,1,10,100,1000))+
      coord_cartesian(ylim=c(0, ifelse(country=='MZ',20,130)))+
      # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
      # labs(x='Date',y='EIR')+
      facet_grid(district~.,labeller = labeller(district = district_labels))+
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
      # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district),alpha=0.2)+
      geom_line(data=df_betaa,aes(x=as.Date(date),y=rel_betaa),size=1)+
      # geom_line(data=df_eir,aes(x=month,y=median),size=1,linetype='dashed')+
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
    indicator_palette <- c(`Relative Rainfall` = 'black',
                           `PG Prevalence` = "#999999",
                           `MG Prevalence` = "#999999",
                           `Relative Incidence` = "#1B9E77",
                           `Relative EIR` = "#D95F02",
                           `Relative Mosquito Emergence` = "#377EB8")

    
    overlay_plot <- ggplot()+
      geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall,color='Relative Rainfall'),size=1)+
      geom_line(data=df_prev,aes(x=date,y=median_pg,color='PG Prevalence'),linetype='solid',linewidth=1)+
      geom_line(data=df_prev,aes(x=date,y=median_mg,color='MG Prevalence'),linetype='dashed',linewidth=1)+
      geom_line(data=df_inc,aes(x=date,y=rel_inc,color='Relative Incidence'),size=1)+
      geom_line(data=df_eir,aes(x=date,y=rel_eir,color='Relative EIR'),size=1)+
      geom_line(data=df_betaa,aes(x=as.Date(date),y=rel_betaa,color='Relative Mosquito Emergence'),size=1)+
      scale_color_manual(values=indicator_palette)+
      scale_x_date(date_labels = "%b %Y")+
      # coord_cartesian(ylim=c(0, 1000))+
      # scale_y_continuous(limits=c(0,500))+
      # scale_y_log10(breaks=c(.1,1,10,100,1000),labels=c(.1,1,10,100,1000))+
      coord_cartesian(ylim=c(0, 1))+
      # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
      # labs(x='Date',y='EIR')+
      facet_grid(.~district,labeller = labeller(district = district_labels))+
      labs(title = 'Indicator Synchrony')+
      theme(legend.title = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1),
            axis.ticks.x = element_line(size = 0.5),
            axis.ticks.length = unit(3, "pt")
      )
    rel_dash <- est_prev_plot_pg+est_prev_plot_mg+est_eir_plot+est_betaa_plot+rel_est_inc_plot+obs_rel_inc_plot+plot_layout(guides = "collect",ncol=6)
    full_dash <- est_prev_plot_pg+est_prev_plot_mg+est_betaa_plot+est_inc_plot+obs_inc_plot+obs_pos_plot+ plot_layout(guides = "collect",ncol=6)
  }
  # print('obs_pos_plot')
  # # obs_test_plot <- ggplot()+
  # #   # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
  # #   geom_line(aes(x=date_ex,y=inc),size=1)+
  # #   # geom_line(data=df_eir,aes(x=month,y=median),size=1,linetype='dashed')+
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
  obs_data_dash <- obs_prev_plot_pg+obs_prev_plot_mg+est_inc_plot+obs_inc_plot+ plot_layout(guides = "collect",ncol=4)
  pres_dash <- est_prev_plot_pg+est_prev_plot_mg+est_inc_plot+obs_inc_plot+ plot_layout(guides = "collect",ncol=4)
  mcmc_dash <- volatility.density+init_EIR.density+plot_layout(ncol=1)
  
  return(list(full_dash = full_dash,
              obs_data_dash = obs_data_dash,
              pres_dash = pres_dash,
              rel_dash = rel_dash,
              mcmc_dash = mcmc_dash,
              overlay_plot = overlay_plot))
  # obs_prev_plot_mg
}
