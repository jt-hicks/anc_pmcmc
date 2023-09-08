create_summary_plots <- function(results,data_list,rainfall,itn_access,
                                 level=c('Region','Council'),
                                 start_pf_time,
                                 date_limits = as.Date(c(NA,NA))){
  date_limits <- as.Date(date_limits)
  itn_access <- itn_access %>%
    rename(itn_access_prop = itn_access)%>%
    mutate(date = as.Date(date_decimal(date)),
           itn_access_prop = ifelse(itn_access_prop > 100,100,itn_access_prop),
           net_access = 'Yes')
  # itn_no_access <- itn_access %>%
  #   mutate(itn_access_prop = 100-itn_access_prop,
  #          net_access = 'No')
  # itn_access <- rbind(itn_access,itn_no_access) %>%
  #   mutate(net_access = factor(net_access,levels = c('No','Yes')))
  dates_list <- lapply(1:length(data_list),
                        function(x,start_pf_time){
                          start_obs <- min(zoo::as.Date(zoo::as.yearmon((data_list[[x]]$month))))#Month of first observation (in Date format)
                          time_origin <- as.Date(paste0(year(start_obs)-1,'-01-01')) #January 1 of year before observation (in Date format)
                          date <- zoo::as.Date(zoo::as.yearmon(data_list[[x]]$month), frac = 0.5) #Convert dates to middle of month
                          t <- as.integer(difftime(date,time_origin,units="days"))

                          return(data.frame(date=date,t=t))
                        })
  
  start_times <- sapply(dates_list, function(x,start_pf_time){
    return(list(t=min(x$t)-start_pf_time,date=min(x$date)-start_pf_time))
  },start_pf_time=start_pf_time)
  
  observed.df <- bind_rows(lapply(1:length(data_list),
                                  function(x){
                                    df <- data_list[[x]]
                                    df$sites <- names(data_list)[x]
                                    return(df)
                                  }))
  observed.df <- addCIs(observed.df,observed.df$positive,observed.df$tested)
  observed.df$date <- as.Date(observed.df$month,frac = 0.5)
  
  mcmc.df <- bind_rows(lapply(1:length(results), 
                              function(x){
                                df <- results[[x]]$mcmc[101:1000,]
                                df$sites <- names(results[x])
                                df$step <- 1:nrow(df)
                                return(df)
                              }))
  
  # rainfall$sites <- sapply(1:nrow(rainfall), function(x) gsub(' Region','',rainfall[x,level]))
  rainfall$month <- as.Date(as.yearmon(rainfall$Month))
  rainfall$yearmon <- as.yearmon(rainfall$Month)

  inc.rainfall.df <- bind_rows(lapply(1:length(results), 
                                      function(x,rainfall,date_key,start_times){
                                        history.df <- as.data.frame(t(results[[x]]$history['clininc_all', 101:1000, -1]))
                                        inc_history <- history.df%>%
                                          dplyr::mutate(t=c((start_times[,x]$t+1):(start_times[,x]$t+nrow(history.df))))%>%
                                          melt(id='t')%>%
                                          dplyr::rename(time=t)%>%
                                          group_by(time)%>%
                                          dplyr::summarise(inc.median=median(value),
                                                           inc.mean=mean(value),
                                                           inc.upper=quantile(value,probs=0.975),
                                                           inc.lower=quantile(value,probs=0.025))%>%
                                          mutate(date = as.Date(c((start_times[,x]$date+1):(start_times[,x]$date + nrow(history.df)))))%>%
                                          mutate(month = zoo::as.Date(zoo::as.yearmon(date),frac=0.5))%>%
                                          ungroup()%>%
                                          group_by(month)%>%
                                          summarise(inc.median=sum(inc.median),
                                                    inc.mean=sum(inc.mean),
                                                    inc.upper=sum(inc.upper),
                                                    inc.lower=sum(inc.lower))%>%
                                          mutate(sites = names(results[x]))

                                        inc_history <- left_join(dates_list[[x]],inc_history,by=join_by(date==month))%>%
                                          mutate(month=as.yearmon(date))
                                        inc_plus_rainfall <- left_join(inc_history,rainfall,by=join_by(month==yearmon,sites==NAME_2))%>%
                                          mutate(rainfall_norm = Rainfall/max(Rainfall,na.rm = TRUE),
                                                 rainfall_maxinc = Rainfall * (max(inc.median,na.rm = TRUE)/max(Rainfall,na.rm = TRUE)))
                                        history.df.prev <- as.data.frame(t(results[[x]]$history['prev_05', 101:1000, -1]))
                                        prev_history <- history.df.prev%>%
                                          dplyr::mutate(t=c((start_times[,x]$t+1):(start_times[,x]$t+nrow(history.df))))%>%
                                          melt(id='t')%>%
                                          dplyr::rename(time=t)%>%
                                          group_by(time)%>%
                                          dplyr::summarise(prev.median=median(value),
                                                           prev.mean=mean(value),
                                                           prev.upper=quantile(value,probs=0.975),
                                                           prev.lower=quantile(value,probs=0.025))%>%
                                          dplyr::mutate(sites = names(results[x]))
                                        prev_history <- left_join(dates_list[[x]],prev_history,by=join_by(t==time))%>%
                                          mutate(month=as.yearmon(date))
                                        all <- left_join(inc_plus_rainfall,prev_history,by=c('month','t','date','sites'))%>%
                                          mutate(prev_maxinc = prev.median * (max(inc.median)),
                                                 upper_maxinc = prev.upper *max(inc.median),
                                                 lower_maxinc = prev.lower *max(inc.median),
                                                 rainfall_maxprev = Rainfall * (max(prev.median)/max(Rainfall)),
                                                 inc_maxinc = inc.median/ max(inc.median),
                                                 upper_maxinc = inc.upper /max(inc.median),
                                                 lower_maxinc = inc.lower /max(inc.median),
                                                 rainfall_maxprev = Rainfall * (max(prev.median)/max(Rainfall)))
                                        return(all)
                                      },rainfall=rainfall,date_key=dates_list,start_times=start_times))
  # inc.rainfall.df <- addCIs(inc.rainfall.df,inc.rainfall.df$positive,inc.rainfall.df$tested)
  
  eir.df <- bind_rows(lapply(1:length(results), 
                             function(x,start_times){
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
                                        date = start_times[['date',x]]+time-1)

                               return(eir_history)
                             },start_time=start_times))
  
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

  province_grid<-as.data.frame(read_csv("./tanz/processed_inputs/tanga_grid.csv"))
  province_grid$name <- province_grid$code <- province_grid$sites

  ylim.prev <- c(0, 0.5)   # ANC prev limits
  ylim.sec <- c(0, 1)    # net access limits
  
  b.prev <- diff(ylim.prev)/diff(ylim.sec)
  a.prev <- ylim.prev[1] - b.prev*ylim.sec[1]
  
  obs_prev_plot <- ggplot(inc.rainfall.df)+
    geom_area(data=itn_access, aes(x=as.Date(date),y=(itn_access_prop/100)*b.prev+a.prev),
              fill="#9ECAE1",alpha = .7)+
    geom_line(aes(x=date,y=rainfall_norm*0.5),col="#1582AD",size=0.8)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper),alpha=0.2)+
    # geom_line(aes(x=month,y=median),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=observed.df,aes(x=date,y=mean),pch = 19,color='#6D6A67')+
    geom_errorbar(data=observed.df,aes(x=date,ymin=lower,ymax=upper),width = 0,color='#6D6A67')+
    # scale_color_manual(values=colors)+
    # scale_fill_manual(values=colors)+
    facet_geo(~ sites, grid = province_grid%>%
                select(row,col,name,code))+
    scale_x_date(date_labels = "'%y")+
    scale_y_continuous("ANC Prevalence", sec.axis = sec_axis(~ (. - a.prev)/b.prev, name = "Proportion of population\nwith ITN access"),limits = c(0,NA)) +
    coord_cartesian(xlim=date_limits, ylim=c(0,0.5))+
    labs(x='Date')
  est_prev_plot <- ggplot(inc.rainfall.df)+
    geom_area(data=itn_access, aes(x=as.Date(date),y=(itn_access_prop/100)*b.prev+a.prev),
              fill="#9ECAE1",alpha = .7)+
    geom_line(aes(x=date,y=rainfall_norm*0.5),col="#1582AD",size=0.8)+
    geom_point(data=observed.df,aes(x=date,y=mean),pch = 19,color='#6D6A67')+
    geom_errorbar(data=observed.df,aes(x=date,ymin=lower,ymax=upper),width = 0,color='#6D6A67')+
    geom_ribbon(aes(x=date,ymin=prev.lower,ymax=prev.upper),alpha=0.2,fill="#EFBB12")+
    geom_line(aes(x=date,y=prev.median),size=1,color="#EFBB12")+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    # scale_color_manual(values=colors)+
    # scale_fill_manual(values=colors)+
    # facet_wrap(~ sites)+
    facet_geo(~ sites, grid = province_grid%>%
                select(row,col,name,code))+
    scale_x_date(date_labels = "'%y")+
    scale_y_continuous("ANC Prevalence", sec.axis = sec_axis(~ (. - a.prev)/b.prev, name = "Proportion of population\nwith ITN access"),limits = c(0,NA)) +
    coord_cartesian(xlim=date_limits,ylim=c(0,0.5))+
    labs(x='Date')

  corr_plot <- ggplot(inc.rainfall.df)+
    geom_point(aes(x=Rainfall,y=inc.median))+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    # geom_point(data=inc.rainfall.df,aes(x=month,y=mean),pch = 19)+
    # geom_errorbar(data=inc.rainfall.df,aes(x=month,ymin=lower,ymax=upper),width = 0)+
    # scale_color_manual(values=colors)+
    # scale_fill_manual(values=colors)+
    facet_wrap(~ sites)+
    # facet_geo(~ sites, grid = province_grid%>%
    #             select(row,col,code,name))+
    labs(x='rainfall',y='Incidence')+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  inc_rf_corr <- inc.rainfall.df%>%
    group_by(sites)%>%
    dplyr::summarise(corr=cor(inc.median,Rainfall))
  
  # corr_map <- tm_shape(tz) +
  #   tm_borders()
  ylim.inc <- c(0, 0.25)   # incidence limits

  b.inc <- diff(ylim.inc)/diff(ylim.sec)
  a.inc <- ylim.inc[1] - b.inc*ylim.sec[1]
  
  inc.rainfall <- ggplot(inc.rainfall.df)+
    geom_area(data=itn_access, aes(x=as.Date(date),y=(itn_access_prop/100)*b.inc+a.inc),
              fill="#9ECAE1",alpha = .7)+
    geom_line(aes(x=date,y=rainfall_norm*0.25),col="#1582AD",size=0.8)+
    geom_ribbon(aes(x=date,ymin=inc.lower,ymax=inc.upper),alpha=0.4,fill="#CE5126")+
    geom_line(aes(x=date,y=inc.median),size=0.8,col="#CE5126")+
    scale_x_date(date_labels = "'%y")+
    scale_y_continuous("Estimated incidence/normalised rainfall", sec.axis = sec_axis(~ (. - a.inc)/b.inc, name = "Proportion of population\nwith ITN access"),limits = c(0,NA)) +
    # facet_wrap(~ sites)+
    facet_geo(~ sites, grid = province_grid%>%
                select(row,col,code,name))+
    xlab("Year")+
    coord_cartesian(xlim=date_limits,ylim=c(0,0.25))
  
  inc.rainfall.3 <- ggplot(inc.rainfall.df)+
    geom_area(data=itn_access, aes(x=as.Date(date),y=(itn_access_prop/100)),
              fill="#9ECAE1",alpha = .7)+
    geom_line(aes(x=date,y=rainfall_norm),col="#1582AD",size=0.8)+
    geom_ribbon(aes(x=date,ymin=lower_maxinc,ymax=upper_maxinc),alpha=0.4,fill="#CE5126")+
    geom_line(aes(x=date,y=inc_maxinc),size=0.8,col="#CE5126")+
    # facet_wrap(~ sites)+
    facet_geo(~ sites, grid = province_grid%>%
                select(row,col,code,name))+
    # geom_ribbon(aes(x=month,ymin=lower_maxinc,ymax=upper_maxinc),alpha=0.2,fill="#414487FF")+
    # geom_line(aes(x=month,y=prev_maxinc),size=1,color="#414487FF")+
    scale_fill_brewer(palette = "Blues") +
    scale_x_date(date_labels = "'%y")+
    scale_y_continuous("Normalised incidence/normalised rainfall", sec.axis = sec_axis(~ .*1, name = "Proportion of population\nwith ITN access"),limits = c(0,NA)) +
    xlab("Year")+
    coord_cartesian(xlim=date_limits,ylim=c(0,1))+
    theme(axis.text.y = element_blank())
  
  ylim.eir <- c(0, max(eir.df$median))   # school prev limits

  b.eir <- diff(ylim.eir)/diff(ylim.sec)
  a.eir <- ylim.eir[1] - b.eir*ylim.sec[1]
  
  eir.rainfall <- ggplot(inc.rainfall.df)+
    geom_area(data=itn_access, aes(x=as.Date(date),y=(itn_access_prop/100)*b.eir+a.eir),
              fill="#9ECAE1",alpha = .7)+
    geom_line(aes(x=date,y=rainfall_norm*max(eir.df$median)),col="#1582AD",size=0.8)+
    geom_ribbon(data=eir.df ,aes(x=date,ymin=lower,ymax=upper),alpha=0.4,fill="#712F79")+
    geom_line(data=eir.df, aes(x=date,y=median),size=0.8,col="#712F79")+
    # facet_wrap(~ sites)+
    facet_geo(~ sites, grid = province_grid%>%
                select(row,col,code,name))+
    # geom_ribbon(aes(x=month,ymin=lower_maxinc,ymax=upper_maxinc),alpha=0.2,fill="#414487FF")+
    # geom_line(aes(x=month,y=prev_maxinc),size=1,color="#414487FF")+
    scale_fill_brewer(palette = "Blues") +
    scale_x_date(date_labels = "'%y")+
    scale_y_continuous("EIR/normalised rainfall", sec.axis = sec_axis(~ (. - a.eir)/b.eir, name = "Proportion of population\nwith ITN access"),limits = c(0,NA)) +
    xlab("Year")+
    coord_cartesian(xlim=date_limits,ylim=c(0,max(eir.df$median)))

    
  return(list(EIR_SD.density = EIR_SD.density,
              log_init_EIR.density =log_init_EIR.density,
              init_EIR.density = init_EIR.density,
              obs_prev_plot = obs_prev_plot,
              est_prev_plot = est_prev_plot,
              inc.rainfall = inc.rainfall,
              inc.rainfall.3 = inc.rainfall.3,
              eir.rainfall = eir.rainfall,
              inc_rf_corr = inc_rf_corr
  ))
}
