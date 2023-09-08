#------------------------------------------------
#' run_model_example
#'
#' \code{run_model_example} runs model using declared age, EIR, ft, country, admin
#'
#' @param age Vector of age brackets.
#'   Default=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60)
#' @param EIR Numeric of annual EIR. Default=10
#' @param ft Numeric of proportion of symptomatic cases recieving treatment.
#'   Default=0.4
#' @param time Numeric of length of time that the model will run for in days.
#'   Default=365
#' @param admin2 Character of admin unit. Default = "Tororo"
#'
#' @importFrom ggplot2 aes scale_colour_manual scale_x_continuous
#'   xlab ylab geom_line ggplot
#' @importFrom reshape2 melt
#' @importFrom odin odin
#' @importFrom stats coef
#' @useDynLib ICDMM
#'
#' @export


run_model_example <- function(age = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60),
                              EIR = 10,
                              ft = 0.4,
                              admin2 = "Tororo",
                              time = 365){
  
  mpl <- model_param_list_create(num_int = 1)
  
  # generate initial state variables from equilibrium solution
  state <- equilibrium_init_create(age_vector = age,
                                   EIR = EIR,
                                   ft = ft,
                                   model_param_list = mpl,
                                   het_brackets = 5,
                                   admin_unit = admin2)
  
  # create odin generator
  generator <- odin_model
  
  # There are many parameters used that should not be passed through
  # to the model.
  state_use <- state[names(state) %in% coef(generator)$name]
  
  # create model with initial values
  mod <- generator(user = state_use, use_dde = TRUE)
  tt <- seq(0, time, 1)
  
  # run model
  mod_run <- mod$run(tt)
  
  # shape output
  out <- mod$transform_variables(mod_run)
  hum <- list("S"=out$S,"T"=out$T,"D"=out$D,"A"=out$A,"U"=out$U,"P"=out$P)
  humsum <- lapply(hum,function(x){apply(x,sum,MARGIN=1)})
  df <- as.data.frame(matrix(unlist(humsum),nrow=length(tt)))
  colnames(df) <- c("S","T","D","A","U","P")
  df$t <- tt
  brks <- seq(from=0,to=time,by=365)
  ret <- ggplot(melt(df,id.vars="t"),
                aes(x=t,y=.data$value,col=.data$variable)) +
    geom_line() +
    scale_colour_manual(
      values=c("#000000","#CC0000","#339900","#3333FF","#CC9933","#CC3399"),
      name="Compartment",
      labels=c("S","T","D","A","U","P")
    ) +
    scale_x_continuous(breaks=brks,labels=0:(length(brks)-1)) +
    xlab("Years") + ylab("Proportion of population")
  return(list("plot"=ret,"dat"=out))
  
}

#------------------------------------------------
#' load_file
#'
#' \code{load_file} loads package file
#'
#' @description Load a file from within the inst/extdata folder of the
#'   ICMDMM package. File extension must be one of .csv, .txt, or .rds.
#'
#' @param name the name of a file within the inst/extdata folder.
#'
#' @importFrom utils read.csv
#'
#' @export

load_file <- function(name) {
  
  # check that valid file extension
  ext <- strsplit(name, "\\.")[[1]]
  ext <- ext[length(ext)]
  if(is.element(ext, c("csv", "rds")) == FALSE){
    stop("file extension not valid")
  }
  
  # get full file path
  name_full <- paste0(getwd(),"/shared/", name)
  
  # read in file
  if (ext == "rds") {
    ret <- readRDS(name_full)
  } else {
    ret <-  read.csv(file=name_full, header=TRUE, sep=",")
  }
  
  return(ret)
}

#------------------------------------------------
#' match clean
#'
#' @param a First string to compare. Default = NULL
#' @param b Second string to compare. Default = NULL
#'
#' @importFrom RecordLinkage levenshteinSim
#'
#' @export

match_clean <- function(a, b){
  
  a <- gsub("[[:punct:][:space:]]", "", tolower(stringi::stri_trans_general(a, "latin-ascii")))
  b <- gsub("[[:punct:][:space:]]", "", tolower(stringi::stri_trans_general(b, "latin-ascii")))
  
  ret <- which(b %in% a)
  
  if(length(ret) == 0){
    distance <- levenshteinSim(a, b)
    ret <- which.max(distance)
  }
  return(ret)
}

#------------------------------------------------
#' match admin region
#'
#' \code{admin_match} Matches the user input admin unit and country with data
#'
#' @param country Character for country within which admin unit is in.
#'   Default = NULL
#' @param admin_unit Character for admin region. Some fuzzy logic will be used to
#'   match. If not provided then no seasonality is introduced. Default = NULL
#' @param admin_units_seasonal Dataframe of seasonality data for country and admin unit
#'
#' @export

admin_match <- function(admin_unit = NULL, country = NULL,
                        admin_units_seasonal){
  
  # intialise admin match as no match
  admin_matches <- 0
  
  if (!is.null(admin_unit)) {
    
    # if there is no country given then search for the admin unit
    if (is.null(country)) {
      
      # find exact match
      admin_matches <- which(tolower(admin_units_seasonal$admin1) %in% tolower(admin_unit))
      
      # if exact does not match try closest match
      if (length(admin_matches) < 1) {
        admin_matches <- match_clean(admin_unit, admin_units_seasonal$admin1)
      } else if (length(admin_matches) > 1){
        stop('Please specify the country of admin unit.  There are multiple with same name.')
      }
      
      # if we do have a country though find that match first and then find admin
    } else {
      
      # first find an exact match
      country_matches <- which(tolower(admin_units_seasonal$country) %in% tolower(country))
      
      if (length(country_matches) < 1) {
        country_name <- admin_units_seasonal$country[match_clean(country, admin_units_seasonal$country)]
        country_matches <- which(tolower(admin_units_seasonal$country) %in% tolower(country_name))
      }
      
      sub_admin_units_seasonal <- admin_units_seasonal[country_matches,]
      
      # find exact match
      admin_sub_matches <- which(tolower(sub_admin_units_seasonal$admin1) %in% tolower(admin_unit))
      
      # if exact does not match try closest match
      if (length(admin_sub_matches) != 1) {
        admin_sub_matches <- match_clean(admin_unit,
                                         sub_admin_units_seasonal$admin1)
      } else if (length(admin_sub_matches) > 1){
        stop('There are multiple admins with same name - check the data file!')
      }
      
      admin_matches <- country_matches[admin_sub_matches]
    }
    
    message("Requested: ", admin_unit,
            "\nReturned: ", admin_units_seasonal$admin1[admin_matches], ", ",
            admin_units_seasonal$country[admin_matches])
  }
  
  return(admin_matches)
}

#------------------------------------------------
#' transform final state of a seasonal model into the initial state of a stochastic model
#'
#' \code{transform_init} Transform final state of a seasonal model into the initial state of a stochastic model
#'
#' @param final_state Final time point of a seasonal model
#'   Default = NULL
#'
#' @export

transform_init <- function(final_state = NULL){
  f <- final_state
  last <- length(f$FOI_init[,1,1])
  list(FOI_eq = f$FOI_init[last,,],
       init_EIR = f$EIR_init[last,1,1],
       init_S = f$S_init[last,,],
       init_T = f$T_init[last,,],
       init_D = f$D_init[last,,],
       init_A = f$A_init[last,,],
       init_U = f$U_init[last,,],
       init_P = f$P_init[last,,],
       init_IB = f$IB_init[last,,],
       init_ID = f$ID_init[last,,],
       init_ICA = f$ICA_init[last,,],
       ICM_age = f$ICM_age_init[last,],
       age_rate = f$age_rate_init[last,],
       het_wt = f$het_wt_init[last,],
       foi_age = f$foi_age_init[last,],
       rel_foi = f$rel_foi_init[last,],
       na = f$na_init[last],
       nh = f$nh_init[last],
       x_I = f$x_I_init[last,],
       omega = f$omega_init[last],
       den = f$den_init[last,],
       age59 = f$age59_init[last],
       age05 = f$age05_init[last],
       age = f$age_init[last,],
       ft = f$ft_init[last],
       age20l = f$age20l_init[last],
       age20u = f$age20u_init[last],
       age_20_factor = f$age_20_factor_init[last],
       eta = f$eta_init[last],
       rA = f$rA_init[last],
       rT = f$rT_init[last],
       rD = f$rD_init[last],
       rU = f$rU_init[last],
       rP = f$rP_init[last],
       uCA = f$uCA_init[last],
       dCA = f$dCA_init[last],
       dB = f$dB_init[last],
       uB = f$uB_init[last],
       dID = f$dID_init[last],
       uD = f$uD_init[last],
       PM = f$PM_init[last],
       phi0 = f$phi0_init[last],
       phi1 = f$phi1_init[last],
       IC0 = f$IC0_init[last],
       kC = f$kC_init[last],
       b0 = f$b0_init[last],
       b1 = f$b1_init[last],
       kB = f$kB_init[last],
       IB0 = f$IB0_init[last],
       aD = f$aD_init[last],
       fD0 = f$fD0_init[last],
       gammaD = f$gammaD_init[last],
       d1 = f$d1_init[last],
       ID0 = f$ID0_init[last],
       kD = f$kD_init[last],
       dE = f$dE_init[last],
       DY = f$DY_init[last],
       volatility = f$volatility_init[last],
       lag_rates = f$lag_rates_init[last],
       max_param = f$max_param_init[last],
       Q0 = f$Q0_init[last],
       state_check = f$state_check_init[last],
       tau1 = f$tau1_init[last],
       tau2 = f$tau2_init[last],
       prev = f$prev[last]
       
  )
}