create_dashboard_plots_nnp <- function(results,
                                       prev_pg,
                                       prev_mg,
                                       prev_all=NULL,
                                       coefs_pg_df,
                                       coefs_mg_df,
                                       incidence,
                                       cs_data_list,
                                       burnin = 50,
                                       chain_length = 1000,
                                       sample_size = 100,
                                       country=c('BF','MZ','NG'),
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
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  district_labels <- c('Banfora - IG2','Gaoua - Standard','Orodara - PBO',
                       'Changara - PBO','Chemba - Standard','Guro - IG2',
                       'Asa - IG2','Ejigbo - Standard','Ife North - PBO','Moro - RG')
  admin_key <- data.frame(district=c('Banfora','Gaoua','Orodara',
                                     'Changara','Chemba','Guro',
                                     'Asa','Ejigbo','Ife North','Moro'),
                          admin=c('Cascades','Sud-Ouest','Haut-Bassins',
                                  'Tete','Sofala','Manica',
                                  'Kwara','Osun','Osun','Kwara'))
  names(district_labels) <- unlist(district_list)
  
  
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
  
  
  mcmc.df <- bind_rows(lapply(start:(start+number-1), 
                              function(x){
                                df <- results[[x]]$mcmc[(burnin+1):chain_length,]
                                df$sites <- names(results[x])
                                df$step <- 1:nrow(df)
                                return(df)
                              }))
  
  df_prev <- bind_rows(lapply(start:(start+number-1),function(i){
    prev_history <- data.frame(t(results[[i]]$history['prev_05', (burnin+1):chain_length, -1]))
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
      mutate(district = districts[[i-start+1]])
  }))
  df_prev_sample <- bind_rows(lapply(start:(start+number-1),function(i){
    prev_history <- data.frame(t(results[[i]]$history['prev_05', (burnin+1):chain_length, -1]))
    prev_history <- prev_history[, sample(ncol(prev_history), sample_size)]
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
      mutate(district = districts[[i-start+1]])
  }))
  
  observed.pg <- bind_rows(lapply(start:(start+number-1),
                                  function(x){
                                    df <- prev_pg[[x]]
                                    df$district <- names(prev_pg[x])
                                    return(df)
                                  }))
  observed.pg <- addCIs(observed.pg,observed.pg$positive,observed.pg$tested)
  observed.pg$date <- as.Date(observed.pg$month,frac = 0.5)

  observed.mg <- bind_rows(lapply(start:(start+number-1),
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
  
  df_inc <- bind_rows(lapply(start:(start+number-1),function(i){
    inc_history <- data.frame(t(results[[i]]$history['clininc_all', (burnin+1):chain_length, -1]))
    
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
      mutate(district = districts[[i-start+1]])
    long_inc_sum <- left_join(dates_list[[i]],long_inc_sum,by=join_by(date==month))%>%
      mutate(month=as.yearmon(date),
             rel_inc = inc.median/max(inc.median))
  }))
  df_inc_sample <- bind_rows(lapply(start:(start+number-1),function(i){
    inc_history <- data.frame(t(results[[i]]$history['clininc_all', (burnin+1):chain_length, -1]))
    
    inc_sample <- inc_history[, sample(ncol(inc_history), sample_size)] %>%
      dplyr::mutate(t=c((start_times[,i]$t+1):(start_times[,i]$t+nrow(inc_history))))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      mutate(date = rep(as.Date(c((start_times[,i]$date+1):(start_times[,i]$date + nrow(inc_history)))),sample_size))%>%
      mutate(month = zoo::as.Date(zoo::as.yearmon(date),frac=0.5))%>%
      mutate(district = districts[[i-start+1]]) %>%
      group_by(district,month,variable)%>%
      dplyr::summarise(value=sum(value))
    
    inc_sample <- left_join(dates_list[[i]],inc_sample,by=join_by(date==month))%>%
      mutate(month=as.yearmon(date))
  }))
  
  df_eir <- bind_rows(lapply(start:(start+number-1),function(i){
    eir_history <- data.frame(t(results[[i]]$history['EIR', (burnin+1):chain_length, -1]))

    long_eir_sum <- eir_history%>%
      mutate(t=c(1:nrow(eir_history)))%>%
      melt(id=c('t'))%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value),
                       mean=mean(value),
                       upper=quantile(value,probs=0.975),
                       lower=quantile(value,probs=0.025))%>%
      dplyr::mutate(district = districts[[i-start+1]],
                    date = start_times[['date',i]]+time-1,
                    rel_eir = median/max(median))
  }))
  df_eir_sample <- bind_rows(lapply(start:(start+number-1),function(i){
    eir_history <- data.frame(t(results[[i]]$history['EIR', (burnin+1):chain_length, -1]))
    eir_sample <- eir_history[, sample(ncol(eir_history), sample_size)] %>%
      mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value,
                    district = districts[[i-start+1]],
                    date = start_times[['date',i]]+time-1)
  }))
  
  rainfall <- read_csv('./nnp/data/NNP_rainfall.csv') %>%
    left_join(admin_key,by=join_by(NAME_1==admin))%>%
    filter(district %in% district_list[[country]])%>%
    group_by(district) %>%
    filter(between(as.yearmon(Month), 
                   as.yearmon(min(bind_rows(dates_list[start:(start+number-1)])$date)), 
                   as.yearmon(max(bind_rows(dates_list[start:(start+number-1)])$date))))%>%
    mutate(rel_rainfall = Rainfall/max(Rainfall))
  
  incidence <- incidence%>%
    filter(district %in% district_list[[country]])%>%
    group_by(district)%>%
    mutate(rel_inc = (inc/max(inc)))
  
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
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
    geom_line(data=df_inc_sample,aes(x=date,y=value*10000,color=district,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_inc,aes(x=date,y=inc.median*10000,color=district,group=district),size=1)+
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
    geom_point(data=observed.mg,aes(x=date,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=observed.mg,aes(x=date,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),pch= 2, size=1.5,color="black")+
    geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color='black')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
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
  est_prev_plot_mg <- annotations[[country]]+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_line(data=df_prev_sample,aes(x=date,y=prev_mg,color=district,group=variable),alpha=0.2)+
    geom_line(data=df_prev,aes(x=date,y=median_mg,color=district,group=district),linewidth=1)+
    geom_point(data=observed.mg,aes(x=date,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=observed.mg,aes(x=date,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),pch= 2, size=1.5,color="black")+
    geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color='black')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    facet_grid(district~.,labeller = labeller(district = district_labels))+
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
    geom_line(data=df_prev_sample,aes(x=date,y=prev_pg,color=district,group=variable),alpha=0.2)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_prev,aes(x=date,y=median_pg,color=district,group=district),linewidth=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=observed.pg,aes(x=date,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=observed.pg,aes(x=date,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),size=1.5,pch=2,color="black")+
    geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color="black")+
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
  obs_prev_plot_pg <- annotations[[country]]+
    geom_point(data=observed.pg,aes(x=date,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=observed.pg,aes(x=date,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),size=1.5,pch=2,color="black")+
    geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color="black")+
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
    geom_line(data=observed.all,aes(x=date,y=tested,color=district,linetype=factor(grav,levels=c('pg','mg'))))+
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
      betaa_history <- data.frame(t(results[[i]]$history['betaa', (burnin+1):chain_length, -1]))
      
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
      betaa_history <- data.frame(t(results[[i]]$history['betaa', (burnin+1):chain_length, -1]))
      betaa_sample <- betaa_history[, sample(ncol(betaa_history), sample_size)] %>%
        mutate(t=c(1:nrow(betaa_history)))%>%
        melt(id='t')%>%
        dplyr::rename(time=t)%>%
        dplyr::mutate(value = value,
                      district = districts[[i-start+1]],
                      date = start_times[['date',i]]+time-1)
    }))
    est_betaa_plot <- annotations[[country]]+
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
    rel_betaa_plot <- annotations[[country]]+
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

    
    overlay_plot <- annotations[[country]]+
      geom_line(data=rainfall, aes(x=as.Date(as.yearmon(Month)),y=rel_rainfall,color='Relative Rainfall'),size=1)+
      geom_line(data=df_prev,aes(x=date,y=median_pg,group=district,color='PG Prevalence'),linetype='solid',linewidth=1)+
      geom_line(data=df_prev,aes(x=date,y=median_mg,group=district,color='MG Prevalence'),linetype='dashed',linewidth=1)+
      geom_line(data=df_inc,aes(x=date,y=rel_inc,group=district,color='Relative Incidence'),size=1)+
      geom_line(data=df_eir,aes(x=date,y=rel_eir,group=district,color='Relative EIR'),size=1)+
      geom_line(data=df_betaa,aes(x=as.Date(date),y=rel_betaa,group=district,color='Relative Mosquito Emergence'),size=1)+
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
