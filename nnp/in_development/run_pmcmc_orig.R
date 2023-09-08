run_pmcmc_orig <- function(data_raw,
                      n_particles=200,
                      proposal_matrix,
                      max_EIR=1000,
                      # EIR_vol,
                      # proposal_dist,
                      # init_EIR = 100,
                      max_steps = 1e7,
                      atol = 1e-3,
                      rtol = 1e-6,
                      n_steps = 500,
                      n_threads = 4,
                      lag_rates = 10,
                      state_check = 0,## Run equilibrium checks
                      # If state_check = 1, returns expected deriv values which should equal 0 and sets stochastic model to have EIR constant at init_EIR
                      # If state_check = 1 and seasonality_on = 1, then the deterministic seasonal model is still run, but theta2 is forced to 1, forcing a constant seasonality profile
                      # If state_check = 0, no values are printed
                      country = NULL,
                      admin_unit = NULL,
                      preyears = 2, #Length of time in years the deterministic seasonal model should run before Jan 1 of the year observations began
                      seasonality_on = 1,  ## state_check = 1 runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                      ## If seasonality_on = 0, runs the stochastic model based on the standard equilibrium solution
                      seasonality_check = 0,##If 1, saves values of seasonality equilibrium
                      seed = 1L){
  ## Modify dates from data
  start_obs <- min(as.Date(data_raw$month))#Month of first observation (in Date format)
  time_origin <- as.Date(ifelse(month(start_obs)!=1,paste0(year(start_obs),'-01-01'),paste0(year(start_obs)-1,'-01-01'))) #January 1 of the first year of observation (in Date format)
  start_stoch <- as.Date(as.yearmon(start_obs)) #Start of stochastic schedule (First day of the month of first observation)
  data_raw_time <- data_raw %>%
    mutate(date = as.Date(as.yearmon(month), frac = 0.5))%>% #Convert dates to middle of month
    mutate(t = as.integer(difftime(date,time_origin,units="days"))) #Calculate date as number of days since January 1 of first year of observation
  initial_time <- min(data_raw_time$t) - 30 #Start particle_filter_data one month before first ime in data
  data_raw_time <- plyr::rbind.fill(data.frame(t=initial_time,date=zoo::as.Date(start_obs - start_pf_time)),data_raw_time)
  data <- mcstate::particle_filter_data(data_raw_time, time = "t", rate = NULL, initial_time = initial_time) #Declares data to be used for particle filter fitting
  # print('Data processed')
  
  # Compare function to calculate likelihood
  compare <- function(state, observed, pars = NULL) {
    # print('in compare function')
    if(is.na(observed$positive)) {return(numeric(length(state[1,])))}
    dbinom(x = observed$positive,
           size = observed$tested,
           prob = state[1,],
           log = TRUE)
  }
  
  ##Output from particle filter
  ##    run: output used for likelihood calculation
  ##    state: output used for visualization
  index <- function(info) {
    list(run = c(prev = info$index$prev),
         state = c(prev = info$index$prev,
                   EIR = info$index$EIR_out,
                   inc = info$index$inc,
                   Sout = info$index$Sout,
                   Tout = info$index$Tout,
                   Dout = info$index$Dout,
                   Aout = info$index$Aout,
                   Uout = info$index$Uout,
                   Pout = info$index$Pout,
                   p_det_out = info$index$p_det_out,
                   phi_out = info$index$phi_out,
                   b_out = info$index$b_out))
  }
  
  ## Provide schedule for changes in stochastic process (in this case EIR)
  ## Converts a sequence of dates (from start_stoch to 1 month after last observation point) to days since January 1 of the first year of observation
  stochastic_schedule <- as.integer(difftime(seq.Date(start_stoch,max(as.Date(data_raw_time$date+30)),by='month'),time_origin,units="days"))#[-1]
  # print('stochastic_schedule assigned')
  
  #Provide age categories, proportion treated, and number of heterogeneity brackets
  init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
  prop_treated <- 0.4
  het_brackets <- 5
  
  #Create model parameter list. Also loads seasonality profile data file to match to desired admin_unit and country
  mpl_pf <- model_param_list_create(init_age = init_age,
                                    pro_treated = prop_treated,
                                    het_brackets = het_brackets,
                                    max_EIR = max_EIR,
                                    state_check = state_check,
                                    lag_rates = lag_rates,
                                    country = country,
                                    admin_unit = admin_unit,
                                    start_stoch = start_stoch,
                                    time_origin = time_origin,
                                    seasonality_on = seasonality_on)
  # print('model parameter list created')
  # print(mpl_pf$state_check)
  # print(mpl_pf$ssa0)
  
  ## If a deterministic seasonal model is needed prior to the stochastic model, this loads the deterministic odin model
  if(seasonality_on == 1){
    season_model <- odin::odin("shared/odin_model_stripped_seasonal.R")
  }
  
  ## Transformation function that calculates initial values for stohastic model
  transform <- function(mpl,season_model){ ## Wraps transformation function in a 'closure' environment so you can pass other parameters that you aren't fitting with the pMCMC
    function(theta) {
      ## theta: particle filter parameters that are being fitted (and so are changing at each MCMC step)
      # print('in transform function')
      init_EIR <- exp(theta[["log_init_EIR"]]) ## Exponentiate EIR since MCMC samples on the log scale for EIR
      EIR_vol <- theta[["EIR_SD"]]
      mpl <- append(mpl_pf,list(EIR_SD = EIR_vol)) ## Add MCMC parameters to model parameter list
      
      ## Run equilibrium function
      state <- equilibrium_init_create_stripped(age_vector = mpl$init_age,
                                                init_EIR = init_EIR,
                                                ft = prop_treated,
                                                model_param_list = mpl,
                                                het_brackets = het_brackets,
                                                state_check = mpl$state_check)
      # print('equilibrium state calculated')
      # print(state)
      ##run seasonality model first if seasonality_on == 1
      if(seasonality_on==1){
        # print('creating seasonality equilirium')
        state_use <- state[names(state) %in% coef(season_model)$name]
        
        # create model with initial values
        mod <- season_model$new(user = state_use, use_dde = TRUE)
        
        # tt <- c(0, preyears*365+as.integer(difftime(mpl$start_stoch,mpl$time_origin,units="days")))
        tt <- seq(0, preyears*365+as.integer(difftime(mpl$start_stoch,mpl$time_origin,units="days")),length.out=100)
        
        # run seasonality model
        mod_run <- mod$run(tt, verbose=FALSE,step_size_max=9)
        
        # shape output
        out <- mod$transform_variables(mod_run)
        # windows(10,8)
        # plot(out$t,out$prev,type='l')
        # View(out)
        
        # Transform seasonality model output to match expected input of the stochastic model
        init4pmcmc <- transform_init(out)
        # print(init4pmcmc)
        # cat('prev equilibrium: ',state$prev,'\n')
        # cat('prev seasonal: ',init4pmcmc$prev,'\n')
        
        #Print some equilibrium checks if state_check==1
        if(state_check==1){
          # print('running equilibrium checks')
          H <- sum(init4pmcmc$init_S) + sum(init4pmcmc$init_T) + sum(init4pmcmc$init_D) + sum(init4pmcmc$init_A) + sum(init4pmcmc$init_U) + sum(init4pmcmc$init_P)
          
          deriv_S11 <- -init4pmcmc$FOI_eq[1,1]*init4pmcmc$init_S[1,1] + init4pmcmc$rP*init4pmcmc$init_P[1,1] + init4pmcmc$rU*init4pmcmc$init_U[1,1] +
            init4pmcmc$eta*H*init4pmcmc$het_wt[1] - (init4pmcmc$eta+init4pmcmc$age_rate[1])*init4pmcmc$init_S[1,1]
          cat('deriv S check: ',deriv_S11,'\n')
          b <- init4pmcmc$b0 * ((1-init4pmcmc$b1)/(1+(init4pmcmc$init_IB[1,1]/init4pmcmc$IB0)^init4pmcmc$kB)+init4pmcmc$b1)
          print(b)
          EIR_eq11 <- init4pmcmc$init_EIR/365 * init4pmcmc$rel_foi[1] * init4pmcmc$foi_age[1]
          print(EIR_eq11)
          FOI_lag <- EIR_eq11 * (if(init4pmcmc$init_IB[1,1]==0) init4pmcmc$b0 else b)
          print(FOI_lag)
          deriv_FOI111 <- (init4pmcmc$lag_rates/init4pmcmc$dE)*FOI_lag - (init4pmcmc$lag_rates/init4pmcmc$dE)*init4pmcmc$FOI_eq[1,1]
          cat('deriv FOI check: ',deriv_FOI111,'\n')
          
          # cat('S check: ',init4pmcmc$init_S-state_use$init_S,'\n')
          # cat('T check: ',init4pmcmc$init_T-state_use$init_T,'\n')
          # cat('D check: ',init4pmcmc$init_D-state_use$init_D,'\n')
          # cat('A check: ',init4pmcmc$init_A-state_use$init_A,'\n')
          # cat('U check: ',init4pmcmc$init_U-state_use$init_U,'\n')
          # cat('P check: ',init4pmcmc$init_P-state_use$init_P,'\n')
          cat('init_EIR: ',state_use$init_EIR,'\n')
          cat('init_EIR seasonal: ',init4pmcmc$init_EIR,'\n')
          cat('prev seasonal: ',state$prev,'\n')
          cat('prev seasonal: ',init4pmcmc$prev,'\n')
          
          # mpl['init_EIR'] <- NULL
          # View(init4pmcmc)
          # View(state_use)
        }
        return(init4pmcmc) #Append all parameters from model parameter list for stochastic model
      }
      else{
        return(state)
      }
    }
  }
  
  ## Load stochastic model in odin.dust  
  # print('about to load stochastic model')
  model <- odin.dust::odin_dust("shared/odinmodelmatchedstoch.R")
  # print('loaded stochastic model')
  
  set.seed(seed) #To reproduce pMCMC results
  
  ### Set particle filter
  # print('about to set up particle filter')
  pf <- mcstate::particle_filter$new(data, model, n_particles, compare,
                                     index = index, seed = seed,
                                     stochastic_schedule = stochastic_schedule,
                                     ode_control = dust::dust_ode_control(max_steps = max_steps, atol = atol, rtol = rtol),
                                     n_threads = n_threads)
  # print('set up particle filter')
  
  # print('about to set up pmcmc control')
  ### Set pmcmc control
  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE,
    n_chains = 1,
    n_workers = 1,
    n_threads_total = n_threads,
    rerun_every = 50,
    rerun_random = TRUE)
  # print('set up pmcmc control')
  
  ### Set pmcmc parameters
  EIR_SD <- mcstate::pmcmc_parameter("EIR_SD", 0.3, min = 0,max=2.5,
                                     prior = function(p) dlnorm(p, meanlog = -.2, sdlog = 0.5, log = TRUE))
  log_init_EIR <- mcstate::pmcmc_parameter("log_init_EIR", 1.5, min = -8.5, max = 8.5,
                                           prior = function(p) dnorm(p, mean = 0, sd = 10, log = TRUE) + p) #Add p to adjust for sampling on log scale
  
  pars = list(EIR_SD = EIR_SD, log_init_EIR = log_init_EIR) ## Put pmcmc parameters into a list
  
  mcmc_pars <- mcstate::pmcmc_parameters$new(pars,
                                             proposal_matrix,
                                             transform = transform(mpl_pf,season_model)) ## Calls transformation function based on pmcmc parameters
  # print('parameters set')
  ### Run pMCMC
  # print('starting pmcmc run')
  start.time <- Sys.time()
  pmcmc_run <- mcstate::pmcmc(mcmc_pars, pf, control = control)
  run_time <- difftime(Sys.time(),start.time,units = 'secs')
  print(run_time)
  pars <- pmcmc_run$pars
  probs <- pmcmc_run$probabilities
  mcmc <- coda::as.mcmc(cbind(probs, pars))
  
  ##Save seasonality equilibrium trajectories if checking equilibrium
  seas_pretime <- NULL
  if(seasonality_check==1){
    check_seasonality <- function(theta,mpl_pf,season_model){
      init_EIR <- exp(theta[["log_init_EIR"]]) ## Exponentiate EIR since MCMC samples on the log scale for EIR
      EIR_vol <- theta[["EIR_SD"]]
      mpl <- append(mpl_pf,list(EIR_SD = EIR_vol)) ## Add MCMC parameters to model parameter list
      
      ## Run equilibrium function
      state <- equilibrium_init_create_stripped(age_vector = mpl$init_age,
                                                init_EIR = init_EIR,
                                                ft = prop_treated,
                                                model_param_list = mpl,
                                                het_brackets = het_brackets,
                                                state_check = mpl$state_check)
      # print(state)
      ##run seasonality model first if seasonality_on == 1
      state_use <- state[names(state) %in% coef(season_model)$name]
      
      # create model with initial values
      mod <- season_model$new(user = state_use, use_dde = TRUE)
      
      # tt <- c(0, preyears*365+as.integer(difftime(mpl$start_stoch,mpl$time_origin,units="days")))
      tt <- seq(0, preyears*365+as.integer(difftime(mpl$start_stoch,mpl$time_origin,units="days")),length.out=100)
      
      # run seasonality model
      mod_run <- mod$run(tt, verbose=FALSE,step_size_max=9)
      
      # shape output
      out <- mod$transform_variables(mod_run)
      out.df <- data.frame(t=out$t,
                           prev = out$prev,
                           prev_all = out$prev_all,
                           inc05 = out$inc05,
                           inc = out$inc)
      return(out.df)
    }
    
    seas_pretime <- lapply(1:nrow(pars), function(x) check_seasonality(theta=pars[x,],mpl_pf=mpl_pf,season_model=season_model))
  }
  to_return <- list(threads = n_threads,
                    particles = n_particles,
                    run_time = run_time,
                    mcmc = as.data.frame(mcmc),
                    pars = as.data.frame(pars),
                    probs = as.data.frame(probs),
                    times = pmcmc_run$trajectories$time,
                    history = pmcmc_run$trajectories$state,
                    seas_history = seas_pretime)
  
  return(to_return)
}