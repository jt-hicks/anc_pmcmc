gen_seasonal_sim <- function(init_EIR=100,
                             max_param=125,
                             model_file= "shared/odin_model_stripped_seasonal.R",
                             country = 'Burkina Faso',
                             admin_unit = 'Cascades'){
  #Provide age categories, proportion treated, and number of heterogeneity brackets
  #Create model parameter list. Also loads seasonality profile data file to match to desired admin_unit and country
  cat('Country = ',country,' Admin = ',admin_unit,' init_EIR = ',init_EIR,'\n')
  init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
  
  prop_treated <- 0.4
  het_brackets <- 5
  
  ##set up the simulation for the simualted data 
  time<- 10*365
  out_step=1
  
  mpl <- sifter::model_param_list_create(init_EIR = init_EIR,
                                         init_ft = prop_treated,
                                         betaa_times=betaa_times,
                                         betaa_vals=betaa_vals,
                                         country=country,
                                         admin_unit = admin_unit,
                                         lag_rates = 10,
                                         volatility = volatility,
                                         max_param = max_param,
                                         state_check = 0
  )
  
  pars <- sifter::equilibrium_init_create_stripped(age_vector = init_age,
                                                   init_EIR = init_EIR,
                                                   ft = prop_treated,
                                                   model_param_list = mpl,
                                                   het_brackets = het_brackets)
  
  ##The malaria model 
  generator <- odin(model_file)
  state_use <- pars[names(pars) %in% coef(generator)$name]
  
  # create model with initial values
  mod <- generator(user = state_use, use_dde = TRUE)
  tt <- seq(0, time, out_step)
  
  # run the simulation to base the data
  mod_run <- mod$run(tt, step_max_n = 1e7,
                     atol = 1e-5,
                     rtol = 1e-5)
  
  # shape output
  out <- mod$transform_variables(mod_run)
  
  # plot data and generate data
  # plot(out$t,out$betaa_out,type='l',col="red",ylim=c(0,125))
  # lines(out$t,out$prev*100,col="blue",lwd=4)
  out_df <- data.frame(t=out$t,
                       prev=out$prev,
                       date=as.character(seq.Date(from = as.Date('2015-01-01'),by = 'day',length.out = length(tt))))
  months <- unique(as.yearmon(out_df$date))
  midmonth_dates <- data.frame(date=as.character(as.Date(months,frac=0.5)))
  monthly_data <- left_join(midmonth_dates,out_df,by='date')
  ##Use only last 2 years for simulated data
  monthly_data <- monthly_data[(nrow(monthly_data)-23):nrow(monthly_data),]
  
  monthly_data$tested<-round(rnorm(nrow(monthly_data),500,50))
  monthly_data$positive<-rbinom(nrow(monthly_data),monthly_data$tested,monthly_data$prev)
  monthly_data$month <- zoo::as.yearmon(monthly_data$date)
  if(any(monthly_data$positive>monthly_data$tested)) {stop('Number of positive is greater than number tested')}
  
  true_data <- data.frame(t=out$t,
                          date=seq.Date(from = as.Date('2015-01-01'),by = 'day',length.out = length(tt)),
                          prev05_true=out$prev,
                          EIR_true=out$EIR_out,
                          betaa_true=out$betaa_out,
                          inc05_true=out$inc05,
                          inc_all_true = out$inc,
                          prev_all_true = out$prev_all)
  
  return(list(volatility=volatility,
              init_EIR = init_EIR,
              true_data=true_data,
              data_raw=monthly_data))
}
