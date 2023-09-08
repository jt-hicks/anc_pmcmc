get_prev_from_log_odds<-function(log_odds){
  return(exp(log_odds)/(1+exp(log_odds)))
}

## fn to return odds from prevalence
get_odds_from_prev<-function(prev){
  return(prev/(1-prev))
}

date_convert <- function(x,start_pf_time){
  start_obs <- min(zoo::as.Date(zoo::as.yearmon((prev_pg[[x]]$month))))#Month of first observation (in Date format)
  time_origin <- as.Date(paste0(year(start_obs)-1,'-01-01')) #January 1 of year before observation (in Date format)
  date <- zoo::as.Date(zoo::as.yearmon(prev_pg[[x]]$month), frac = 0.5) #Convert dates to middle of month
  t <- as.integer(difftime(date,time_origin,units="days"))
  
  return(data.frame(date=date,t=t))
}

mcmc_history_process <- function(x,par_name,burnin,chain_lgth,comparison){
  ## Get the trajectories from the pmcmc output
  par_history <- data.frame(t(results[[i]]$history[par_name, (burnin+1):chain_lgth, -1]))
  
  ##Only need an instantaneous measure of prevalence so just take modelled prevalence
  ##  at time of observation (to match how model is fit in comparison function)
  if(str_detect(par_name,'prev')){
    par_history$t <- c((start_times[,i]$t+1):(start_times[,i]$t+nrow(par_history)))
    par_history <- left_join(dates_list[[i]],par_history,by=join_by(t==t))
  }

  long_par_history <- par_history%>%
    melt(id=c('t','date'))%>%
    dplyr::rename(time=t)
  
  if(str_detect(par_name,'prev_05')){
    par_history <- convert_prev(long_par_history,comparison=comparison,coefs_pg_df=coefs_pg_df,coefs_mg_df=coefs_mg_df)
  }
  
  frame <- unique(long_par_history$variable)
  selection <- sample(length(frame), 100)
  sample_history <- long_par_history %>%
    filter(variable %in% frame[selection])
}

convert_prev <- function(df,comparison,coefs_pg_df,coefs_mg_df=NULL){
  logodds_child <- log(get_odds_from_prev(df$value))
  df$prev_pg <- rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
    prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
  })))
  if(comparison %in% c('pgmg','pgsg')){
    df$prev_mg <- rowMedians(as.matrix(sapply(1:nrow(coefs_mg_df), function(x){
      prev_preg <- get_prev_from_log_odds(logodds_child+coefs_mg_df$gradient[x]*(logodds_child-coefs_mg_df$av_lo_child[x])+coefs_mg_df$intercept[x])
    })))
  }
  return(df)
}