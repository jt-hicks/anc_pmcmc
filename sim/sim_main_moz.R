##Load libraries
library('tidyverse')
library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(didehpc)
library(pkgdepends)
library("coda")
library(binom)
library(ggplot2)
library(bayesplot)
library(reshape2)
library(ggpubr)
library(zoo)
library(patchwork)
library(RColorBrewer)
theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))

##Load sources
#Required functions
source('sim/data_gen_moz.R')

##Generate data sets of random walk mosquito emergence
test <- data_gen_moz(volatility = 1)
View(test$data_raw)

pars <- expand.grid(EIR_levels = c(1,10,100,1000), vol_levels = c(0.25,0.75,1.0,2.0))
pars_list <- bind_rows(pars,pars,pars,pars,pars,pars)

rand_walk_data <- lapply(1:nrow(pars_list),function(i){
    data_gen_moz(volatility = pars_list$vol_levels[i], init_EIR = pars_list$EIR_levels[i])
})
saveRDS(rand_walk_data,'sim/rand_walk_data_220823.rds')
rand_walk_data <- readRDS('sim/rand_walk_data_220823.rds')
View(rand_walk_data)
# Check data set
# Turn into a large dataframe
rand_walk_true_df <- bind_rows(lapply(1:length(rand_walk_data),function(i){
  data <- rand_walk_data[[i]]$true_data
  data$init_EIR <- rand_walk_data[[i]]$init_EIR
  data$volatility <- rand_walk_data[[i]]$volatility
  data$trial <- ceiling(i/16)
  return(data)
}))
rand_walk_raw_df <- bind_rows(lapply(1:length(rand_walk_data),function(i){
  data <- rand_walk_data[[i]]$data_raw
  data$init_EIR <- rand_walk_data[[i]]$init_EIR
  data$volatility <- rand_walk_data[[i]]$volatility
  data$trial <- ceiling(i/16)
  return(data)
}))
windows(15,10)
ggplot(rand_walk_true_df)+
  geom_line(aes(x=t,y=prev05_true,color=factor(trial)))+
  facet_grid(init_EIR~volatility)+
  scale_color_brewer(palette = 'Paired')+
  labs(color='Trial',x='Days',y='Prevalence u5')
ggplot(rand_walk_true_df)+
  geom_line(aes(x=t,y=betaa_true,color=factor(trial)))+
  facet_grid(init_EIR~volatility)+
  scale_color_brewer(palette = 'Paired')+
  labs(color='Trial',x='Days',y='Mosquito Emergence')
ggplot(rand_walk_true_df)+
  geom_line(aes(x=t,y=EIR_true,color=factor(trial)))+
  facet_grid(init_EIR~volatility)+
  scale_color_brewer(palette = 'Paired')+
  labs(color='Trial',x='Days',y='EIR')

##Generate data sets of seasonally fluctuating mosquito emergence
source('sim/data_gen_seas.R')
test_ky <- gen_seasonal_sim(init_EIR=100,country='Kenya',admin_unit='Siaya')
plot(test_ky$true_data$t,test_ky$true_data$EIR_true)
test_mw <- gen_seasonal_sim(init_EIR=100,country='Malawi',admin_unit='Blantyre')
plot(test_mw$true_data$t,test_mw$true_data$EIR_true)
test_bf <- gen_seasonal_sim(init_EIR=100,country='Burkina Faso',admin_unit='Plateau-Central')
plot(test_bf$true_data$t,test_bf$true_data$EIR_true)
test_gm <- gen_seasonal_sim(init_EIR=100,country='Gambia',admin_unit='Upper River')
plot(test_gm$true_data$t,test_gm$true_data$EIR_true)
test_gh <- gen_seasonal_sim(init_EIR=100,country='Ghana',admin_unit='Upper East')
plot(test_gh$true_data$t,test_gh$true_data$EIR_true)
test_ml <- gen_seasonal_sim(init_EIR=100,country='Mali',admin_unit='Koulikoro')
plot(test_ml$true_data$t,test_ml$true_data$EIR_true)

admins_seasonal <- data.frame(country = c('Kenya', #20%
                                          'Malawi', #80% #Two sites are in Blantyre, one is in neighboring Chikwawa
                                          'Burkina Faso', #80%
                                          'Gambia', #90%
                                          'Ghana', #70%
                                          'Mali'),
                              admin = c('Siaya', #20%
                                        'Blantyre', #80% #Two sites are in Blantyre, one is in neighboring Chikwawa
                                        'Plateau-Central', #80%
                                        'Upper River', #90%
                                        'Upper East', #70%
                                        'Koulikoro')) #80%##There were 3 Mali sites in the trial. Chose the admin area in the middle
eir_vals <- c(1,10,100,1000)
pars_seas <- bind_rows(lapply(1:4, function(x){
  admins_seasonal$EIR <- eir_vals[x]
  return(admins_seasonal)
}))

sim_seas_data <- lapply(1:nrow(pars_seas),function(i){
  gen_seasonal_sim(init_EIR = pars_seas$EIR[i],
                   country = pars_seas$country[i],
                   admin_unit = pars_seas$admin[i])
})
saveRDS(sim_seas_data,'sim/sim_seas_data_220823.rds')

##Set up cluster
##Run long runs
##Configure cluster settings
ctx_sim <- context::context_save("T:/jth/contexts.sim",
                                  packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde','RecordLinkage','plyr'),
                                  package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate','jt-hicks/sifter@issue-10')))
config_sim <- didehpc::didehpc_config(template = "24Core",cores =6, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
# config_dide <- didehpc::didehpc_config(template = "8Core",cores =1, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj_sim <- didehpc::queue_didehpc(ctx_sim,config = config_sim)
obj_sim$cluster_load(TRUE)
obj_sim$login()

##Submit datasets to cluster
sim_rw_bulk_betaa_1108 <- obj_sim$enqueue_bulk(1:length(rand_walk_data), function(i,sim_data){
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    proposal_matrix = diag(0.1, 2),
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 2,
                    lag_rates = 10,
                    seasonality_on = 0,
                    state_check = 0,
                    start_pf_time = 30,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 1008)
},sim_data=rand_walk_data)

test_sifter <- obj_sim$enqueue(sifter::run_pmcmc(data_raw = sifter::data_sim, #I've added data_sim to the package for an easy test
                  n_particles = 100,
                  proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                  max_param=1000,
                  max_steps = 1e7,
                  atol = 1e-5,
                  rtol = 1e-6,
                  n_steps = 100,
                  n_threads = 2,
                  lag_rates = 10,
                  country = 'Burkina Faso',
                  admin_unit = 'Cascades',
                  seasonality_on = 0,
                  state_check = 0,
                  seasonality_check = 0,
                  comparison = 'u5',
                  stoch_param = 'betaa'))
test_sifter$status()
test_sifter$log()
obj_sim
sim_rw_bulk_betaa_1108$status() #'nonmathematic_dromaeosaur' submitted 5:40pm on 11 Aug
sim_rw_bulk_betaa_1108$tasks[[1]]$log()
sim_rw_bulk_betaa_1108$times()
# obj_sim$unsubmit(sim_rw_bulk_betaa_1108$ids)
obj_sim$login()
sim_rw_bulk_test <- obj_sim$enqueue_bulk(1:length(rand_walk_data), function(i,sim_data){
  print(sim_data[[i]]$data_raw$prev)},sim_data=rand_walk_data)
sim_rw_bulk_test$status()
sim_rw_bulk_test$tasks[[5]]$result()

sim_rw_bulk_betaa_1408_test <- obj_sim$enqueue_bulk(1:3, function(i,sim_data){
  print(sim_data[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    proposal_matrix = diag(0.1, 2),
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 2,
                    lag_rates = 10,
                    seasonality_on = 0,
                    state_check = 0,
                    start_pf_time = 30,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 1008)
},sim_data=rand_walk_data)
sim_rw_bulk_betaa_1408$status() #'sound_neonbluehermitcrab' submitted at 2:56pm on 14 Aug
sim_rw_bulk_betaa_1408_test$status()
obj_sim$unsubmit(sim_rw_bulk_betaa_1408$ids)
sim_rw_bulk_betaa_1408$tasks[['c258188b621c5cc434fb93f4d0e072b7']]$log()

sim_rw_bulk_betaa_1708_test <- obj_sim$enqueue_bulk(1:5, function(i,sim_data){
  print(sim_data[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    proposal_matrix = diag(0.1, 2),
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 2,
                    lag_rates = 10,
                    seasonality_on = 0,
                    state_check = 0,
                    start_pf_time = 30,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 1008)
},sim_data=rand_walk_data)
sim_rw_bulk_betaa_1708_test$status() #'trisyllabical_darklingbeetle' submitted 17 Aug 12:38pm
sim_rw_bulk_betaa_1708_test$times()
sim_rw_bulk_betaa_1708_test$tasks[[4]]$log()
# obj_sim$unsubmit(sim_rw_bulk_betaa_1708_test$ids)
cov(sim_rw_bulk_betaa_1708_test$tasks[[1]]$result()$pars)
proposal_mat <- cov(sim_rw_bulk_betaa_1708_test$tasks[[2]]$result()$pars) #0.089089447 0.003156369 0.003156369 0.047585393
cov(sim_rw_bulk_betaa_1708_test$tasks[[3]]$result()$pars)
cov(sim_rw_bulk_betaa_1708_test$tasks[[5]]$result()$pars)

sim_rw_bulk_betaa_2108_test <- obj_sim$enqueue_bulk(1:16, function(i,sim_data,proposal_mat){
  print(sim_data[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    init_EIR = sim_data[[i]]$init_EIR,
                    proposal_matrix = proposal_mat,
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 0,
                    state_check = 0,
                    start_pf_time = 90,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 1008)
},sim_data=rand_walk_data,proposal_mat=proposal_mat)
sim_rw_bulk_betaa_2108_test$status() #palaeoecological_altiplanochinchillamouse submitted at 1:39pm on 21 Aug
sim_rw_bulk_betaa_2108_test$tasks[[12]]$log()
sim_rw_bulk_betaa_2108_test$times()
obj_sim$unsubmit(sim_rw_bulk_betaa_2108_test$ids)

test_5 <- sim_rw_bulk_betaa_2108_test$tasks[[5]]$result()
matplot(test_5$times[-1],t(test_5$history['prev_05',,-1]),col='black',type='l',ylim = c(0,1))
lines(rand_walk_data[[5]]$data_raw$t+350,rand_walk_data[[5]]$data_raw$prev,type='l',col='red')

sim_rw_bulk_betaa_2208_test <- obj_sim$enqueue_bulk(1:16, function(i,sim_data,proposal_mat){
  print(sim_data[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    init_EIR = sim_data[[i]]$init_EIR,
                    proposal_matrix = proposal_mat,
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 0,
                    state_check = 0,
                    start_pf_time = 90,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 1008)
},sim_data=rand_walk_data,proposal_mat=proposal_mat)
sim_rw_bulk_betaa_2208_test$status() #skyborne_slug submitted at 9:47am on 22 Aug
sim_rw_bulk_betaa_2208_test$tasks[['77a7526d9f09691480155d7b64888851']]$log()
sim_rw_bulk_betaa_2208_test$tasks[[8]] <- sim_rw_bulk_betaa_2208_test4.8$tasks[[2]]
prop_mat_rw_list_2208 <- lapply(c(1:16),function(i) cov(sim_seas_bulk_betaa_2208$tasks[[i]]$result()$pars))
all_shapes_rw <- lapply(c(1:16),function(i) cov_var_time(sim_rw_bulk_betaa_2208_test$tasks[[i]]$result()$pars))
all_shapes_rw_mat <- do.call(rbind, all_shapes_rw)
all_shapes_rw_mat[,999]
matplot(t(all_shapes_rw_mat),type='l')

sim_rw_bulk_betaa_2208_test4.8 <- obj_sim$enqueue_bulk(c(4,8), function(i,sim_data,proposal_mat){
  print(sim_data[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    init_EIR = sim_data[[i]]$init_EIR,
                    proposal_matrix = proposal_mat,
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 0,
                    state_check = 0,
                    start_pf_time = 90,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 2308)
},sim_data=rand_walk_data,proposal_mat=proposal_mat)
# obj_sim$unsubmit(sim_rw_bulk_betaa_2208_test4.8$ids)
sim_rw_bulk_betaa_2208_test4.8$status()

sim_rw_bulk_betaa_2508 <- obj_sim$enqueue_bulk(1:16, function(i,sim_data,proposal_mat){
  print(sim_data[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    init_EIR = sim_data[[i]]$init_EIR,
                    proposal_matrix = proposal_mat[[i]],
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 0,
                    state_check = 0,
                    start_pf_time = 90,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 2508)
},sim_data=rand_walk_data,proposal_mat=prop_mat_rw_list_2208)
sim_rw_bulk_betaa_2508$status() #mid_goitered submitted 10:01am on 25 Aug

sim_rw_bulk_betaa_3008_test <- obj_sim$enqueue_bulk(1:16, function(i,sim_data,proposal_mat){
  print(sim_data[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    init_EIR = sim_data[[i]]$init_EIR,
                    proposal_matrix = proposal_mat,
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 0,
                    state_check = 0,
                    start_pf_time = 90,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 1008)
},sim_data=rand_walk_data,proposal_mat=proposal_mat)
sim_rw_bulk_betaa_3008_test$status() #dermic_anchovy submitted at 1:32pm on 30 Aug


sim_seas_bulk_betaa_2208 <- obj_sim$enqueue_bulk(1:24, function(i,sim_data,proposal_mat){
  print(sim_data[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    init_EIR = sim_data[[i]]$init_EIR,
                    proposal_matrix = proposal_mat,
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 1,
                    preyears = 8,
                    state_check = 0,
                    start_pf_time = 90,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 1008)
},sim_data=sim_seas_data,proposal_mat=proposal_mat)
sim_seas_bulk_betaa_2208$status() #deductive_funnelweaverspider submitted at 3:33am on 22 Aug
sim_seas_bulk_betaa_2208$times()
# obj_sim$unsubmit(sim_rw_bulk_betaa_2108_test$ids)
prop_mat_list_2208 <- lapply(c(1:24),function(i) cov(sim_seas_bulk_betaa_2208$tasks[[i]]$result()$pars))
source('./test_scripts/cov_shape.R')
test_cov <- cov_var_time(sim_seas_bulk_betaa_2208$tasks[[1]]$result()$pars)
str(sim_seas_bulk_betaa_2208$tasks[[1]]$result()$pars)
all_shapes <- lapply(c(1:24),function(i) cov_var_time(sim_seas_bulk_betaa_2208$tasks[[i]]$result()$pars))
all_shapes_mat <- do.call(rbind, all_shapes)
matplot(t(all_shapes_mat),type='l')
obj_sim$login()

sim_seas_bulk_betaa_2408 <- obj_sim$enqueue_bulk(1:24, function(i,sim_data,proposal_mat){
  print(sim_data[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    init_EIR = sim_data[[i]]$init_EIR,
                    proposal_matrix = proposal_mat[[i]],
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 1,
                    preyears = 8,
                    state_check = 0,
                    start_pf_time = 90,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 2408)
},sim_data=sim_seas_data,proposal_mat=prop_mat_list_2208)
sim_seas_bulk_betaa_2408$status() # tenuous_curassow submitted at 11:40am on 24 Aug
# obj_sim$unsubmit(sim_seas_bulk_betaa_2408$ids)

sim_seas_bulk_betaa_2408_flat <- obj_sim$enqueue_bulk(1:24, function(i,sim_data,proposal_mat){
  print(sim_data[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    init_EIR = sim_data[[i]]$init_EIR,
                    proposal_matrix = proposal_mat[[i]],
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 0,
                    preyears = 8,
                    state_check = 0,
                    start_pf_time = 90,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 2408)
},sim_data=sim_seas_data,proposal_mat=prop_mat_list_2208)
sim_seas_bulk_betaa_2408_flat$status()

sim_seas_bulk_betaa_3008 <- obj_sim$enqueue_bulk(1:24, function(i,sim_data,proposal_mat){
  print(sim_data[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    init_EIR = sim_data[[i]]$init_EIR,
                    proposal_matrix = proposal_mat,
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 1,
                    preyears = 8,
                    state_check = 0,
                    start_pf_time = 90,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 1008)
},sim_data=sim_seas_data,proposal_mat=proposal_mat)
sim_seas_bulk_betaa_3008$status() #piny_cockatiel submitted at 3:33am on 22 Aug

prop_mat_list_3008 <- lapply(c(1:24),function(i) cov(sim_seas_bulk_betaa_3008$tasks[[i]]$result()$pars))
obj_sim$login()
obj_sim$config
sim_seas_bulk_betaa_0109 <- obj_sim$enqueue_bulk(1:24, function(i,sim_data,proposal_mat){
  print(sim_data[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    init_EIR = sim_data[[i]]$init_EIR,
                    proposal_matrix = proposal_mat[[i]],
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 1,
                    preyears = 8,
                    state_check = 0,
                    start_pf_time = 30,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 1008)
},sim_data=sim_seas_data,proposal_mat=prop_mat_list_3008)
sim_seas_bulk_betaa_0109$status() #transcendent_irishredandwhitesetter submitted at 4:34pm on 1 Sep

sim_seas_bulk_betaa_3008_flat <- obj_sim$enqueue_bulk(1:24, function(i,sim_data,proposal_mat){
  print(sim_data[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = sim_data[[i]]$data_raw,
                    n_particles = 200,
                    init_EIR = sim_data[[i]]$init_EIR,
                    proposal_matrix = proposal_mat,
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 0,
                    preyears = 8,
                    state_check = 0,
                    start_pf_time = 90,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 1008)
},sim_data=sim_seas_data,proposal_mat=proposal_mat)
sim_seas_bulk_betaa_3008_flat$status() #innocuous_indri submitted at 3:33am on 22 Aug

##Check runs
sim_status <- sim_rw_bulk_betaa_1108$status()
complete_index <- which(sim_status=='COMPLETE')
sim_status.df <- enframe(sim_status)
prop.table(table(sim_status.df$value))
sim_complete <- sim_status.df %>%
  filter(value=='COMPLETE')
sim_rw_bulk_betaa_1108_results <- lapply(sim_complete$name, function(i){
  sim_rw_bulk_betaa_1108$tasks[[i]]$result()
})
names(sim_rw_bulk_betaa_1108_results) <- complete_index

sim_seas_bulk_betaa_2208_results <- lapply(1:24, function(i){
  sim_seas_bulk_betaa_2208$tasks[[i]]$result()
})

sim_rw_bulk_betaa_2208_test_results <- lapply(1:16, function(i){
  sim_rw_bulk_betaa_2208_test$tasks[[i]]$result()
})
sim_rw_bulk_betaa_2508_results <- lapply(1:16, function(i){
  sim_rw_bulk_betaa_2508$tasks[[i]]$result()
})
sim_seas_bulk_betaa_2408_results <- lapply(1:24, function(i){
  sim_seas_bulk_betaa_2408$tasks[[i]]$result()
})
sim_rw_bulk_betaa_3008_results <- lapply(1:16, function(i){
  sim_rw_bulk_betaa_3008_test$tasks[[i]]$result()
})
sim_rw_bulk_betaa_3008_results <- lapply(1:16, function(i){
  sim_rw_bulk_betaa_3008_test$tasks[[i]]$result()
})
sim_seas_bulk_betaa_3008_results <- lapply(1:24, function(i){
  sim_seas_bulk_betaa_3008$tasks[[i]]$result()
})
sim_seas_bulk_betaa_3008_flat_results <- lapply(1:24, function(i){
  sim_seas_bulk_betaa_3008_flat$tasks[[i]]$result()
})

##Summarize runs
create_diag_sum <- function(result,country){
  mcmc_chain <- result$mcmc[501:1000]
  ar <- 1 - coda::rejectionRate(as.mcmc(mcmc_chain))
  ess <- coda::effectiveSize(as.mcmc(mcmc_chain))
  return(data.frame(ar = ar,
                    ess = ess,
                    param = names(ar),
                    country = country))
}
create_diag_sum <- function(result,country){
  mcmc_chain <- result$mcmc[501:1000,]
  ar <- 1 - coda::rejectionRate(as.mcmc(mcmc_chain))
  ess <- coda::effectiveSize(as.mcmc(mcmc_chain))
  return(data.frame(ar = ar,
                    ess = ess,
                    param = names(ar),
                    country = country))
}
create_diag_figs <- function(result,country){
  print('acceptance rate')
  ar <- 1 - coda::rejectionRate(as.mcmc(result$mcmc))
  print(ar)
  print('effective size')
  ess <- coda::effectiveSize(as.mcmc(result$mcmc))
  print(ess)
  
  title <- paste0('Diagnostic plots for seasonal model - ',country)
  pars_list <- c('log_prior','log_likelihood','log_posterior','volatility','init_betaa')
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
test <- create_diag_figs(sim_rw_bulk_betaa_1708_test$tasks[[5]]$result(),country = 1)
test_5 <- sim_rw_bulk_betaa_1708_test$tasks[[5]]$result()
matplot(test_5$times[-1],t(test_5$history['prev_05',,-1]),col='black',type='l',ylim = c(0,1))
lines(rand_walk_data[[5]]$data_raw$t+350,rand_walk_data[[5]]$data_raw$prev,type='l',col='red')
test_2 <- sim_rw_bulk_betaa_1708_test$tasks[[2]]$result()
matplot(test_2$times[-1],t(test_2$history['prev_05',,-1]),col='black',type='l',ylim = c(0,1))
lines(rand_walk_data[[2]]$data_raw$t+350,rand_walk_data[[2]]$data_raw$prev,type='l',col='red')
matplot(test_2$times[],t(test_2$history['betaa',,]),col='black',type='l')
matplot(test_2$times[],t(test_2$history['Dout',,]),col='black',type='l')
matplot(test_2$times[],t(test_2$history['spz_rate',,]),col='black',type='l')
matplot(test_2$times[],t(test_2$history['clinic_all',,]),col='black',type='l')
matplot(test_2$times[],t(test_2$history['phi_out',,]),col='black',type='l')
matplot(test_2$times[],t(test_2$history['p_det_out',,]),col='black',type='l')
test_1 <- create_diag_figs(sim_rw_bulk_betaa_3008_test$tasks$`d0ef4d2cbdb8bca6eb1c4f5813309187`$result(),country = 1)
test <- create_diag_figs(sim_seas_bulk_betaa_3008$tasks[['f27477038e9b1114200588b43a9c56de']]$result(),country = 1)
test
diag_sim_rw_bulk_betaa_1108_plots <- lapply(1:length(sim_rw_bulk_betaa_1108_results),function(i) create_diag_figs(sim_rw_bulk_betaa_1108_results[[i]],country = complete_index[[i]]))

diag_sim_seas_bulk_betaa_2208 <- lapply(1:length(sim_seas_bulk_betaa_2208_results),function(i) create_diag_figs(sim_seas_bulk_betaa_2208_results[[i]],country = i))
diag_sim_seas_bulk_betaa_2208

for(i in c(1:length(diag_sim_seas_bulk_betaa_2208))){
  ggsave(paste0('Q:/anc_pmcmc/sim/figures/diag/',i,'_seas_bulk_betaa_2208.pdf'),plot=diag_sim_seas_bulk_betaa_2208[[i]], width = 10, height = 7)
}
results <- sim_rw_bulk_betaa_1108_results[[1]]

diag_sim_rw_bulk_betaa_2208 <- lapply(1:length(sim_rw_bulk_betaa_2208_test_results),function(i) create_diag_figs(sim_rw_bulk_betaa_2208_test_results[[i]],country = i))
for(i in c(1:length(diag_sim_rw_bulk_betaa_2208))){
  ggsave(paste0('Q:/anc_pmcmc/sim/figures/diag/',i,'_rw_bulk_betaa_2208.pdf'),plot=diag_sim_rw_bulk_betaa_2208[[i]], width = 10, height = 7)
}

diag_sim_seas_bulk_betaa_2408 <- lapply(1:length(sim_seas_bulk_betaa_2408_results),function(i) create_diag_figs(sim_seas_bulk_betaa_2408_results[[i]],country = i))
for(i in c(1:length(diag_sim_seas_bulk_betaa_2408))){
  ggsave(paste0('Q:/anc_pmcmc/sim/figures/diag/',i,'_seas_bulk_betaa_2408.pdf'),plot=diag_sim_seas_bulk_betaa_2408[[i]], width = 10, height = 7)
}

sim_rw_bulk_betaa_3008_results
diag_sim_rw_bulk_betaa_3008 <- lapply(1:length(sim_rw_bulk_betaa_3008_results),function(i) create_diag_figs(sim_rw_bulk_betaa_3008_results[[i]],country = i))

for(i in c(1:length(diag_sim_rw_bulk_betaa_3008))){
  ggsave(paste0('Q:/anc_pmcmc/sim/figures/diag/',i,'_rw_bulk_betaa_3008.pdf'),plot=diag_sim_rw_bulk_betaa_3008[[i]], width = 10, height = 7)
}

diag_sim_seas_bulk_betaa_3008 <- lapply(1:length(sim_seas_bulk_betaa_3008_results),function(i) create_diag_figs(sim_seas_bulk_betaa_3008_results[[i]],country = i))
for(i in c(1:length(diag_sim_seas_bulk_betaa_3008))){
  ggsave(paste0('Q:/anc_pmcmc/sim/figures/diag/',i,'_seas_bulk_betaa_3008.pdf'),plot=diag_sim_seas_bulk_betaa_3008[[i]], width = 10, height = 7)
}
matplot(sim_seas_bulk_betaa_2408_results[[10]]$times[-1],t(sim_seas_bulk_betaa_2408_results[[10]]$history['prev_05',,-1]),col='black',type='l',ylim = c(0,1))
lines(sim_seas_data[[10]]$data_raw$t-2557,sim_seas_data[[10]]$data_raw$prev,type='l',col='red')
matplot(sim_seas_bulk_betaa_2408_results[[24]]$times[-1],t(sim_seas_bulk_betaa_2408_results[[24]]$history['betaa',,-1]),col='black',type='l')

diag_sim_rw_bulk_betaa_2508_sum <- bind_rows(lapply(1:length(sim_rw_bulk_betaa_2508_results),function(i) create_diag_sum(sim_rw_bulk_betaa_2508_results[[i]],country = i)))
diag_sim_rw_bulk_betaa_2508_sum$model_type <- 'randomwalk'
pars_list$country <- c(1:nrow(pars_list))
diag_sim_rw_bulk_betaa_2508_sum <- left_join(diag_sim_rw_bulk_betaa_2508_sum,pars_list,by='country')
diag_sim_rw_bulk_betaa_2508_sum$prop_mat <- 'informed'

diag_sim_rw_bulk_betaa_2208_sum <- bind_rows(lapply(1:length(sim_rw_bulk_betaa_2208_test_results),function(i) create_diag_sum(sim_rw_bulk_betaa_2208_test_results[[i]],country = i)))
diag_sim_rw_bulk_betaa_2208_sum$model_type <- 'randomwalk'
pars_list$country <- c(1:nrow(pars_list))
diag_sim_rw_bulk_betaa_2208_sum <- left_join(diag_sim_rw_bulk_betaa_2208_sum,pars_list,by='country')
diag_sim_rw_bulk_betaa_2208_sum$prop_mat <- 'naive'

sim_rw_both_diag <- rbind(diag_sim_rw_bulk_betaa_2208_sum,diag_sim_rw_bulk_betaa_2508_sum)

ggplot(sim_rw_both_diag[sim_rw_both_diag$param=='init_betaa',])+
  geom_point(aes(x=as.factor(EIR_levels),y=ess,color=prop_mat))+
  facet_wrap(.~vol_levels)
ggplot(sim_rw_both_diag[sim_rw_both_diag$param=='init_betaa',])+
  geom_point(aes(x=as.factor(EIR_levels),y=ar,color=prop_mat))+
  facet_wrap(.~vol_levels)
diag_sim_rw_bulk_betaa_2508 <- lapply(1:length(sim_rw_bulk_betaa_2508_results),function(i) create_diag_figs(sim_rw_bulk_betaa_2508_results[[i]],country = i))
for(i in c(1:length(diag_sim_rw_bulk_betaa_2508))){
  ggsave(paste0('Q:/anc_pmcmc/sim/figures/diag/',i,'_rw_bulk_betaa_2508.pdf'),plot=diag_sim_rw_bulk_betaa_2508[[i]], width = 10, height = 7)
}
matplot(test_2$times[-1],t(test_2$history['prev_05',,-1]),col='black',type='l',ylim = c(0,1))

## get median q1 and 3-rd quartile q2 of p
fn <- function(par){
  a <- par[1]
  b <- par[2]
  print(par)
  print(pgamma(c(q025, q50, q975), a, b))
  return(sum((pgamma(c(q025, q50, q975), a, b) - c(0.025,0.5,0.975))^2))
}
q50 <- 1
q025 <- 0.1
q975 <- 2.5
optim(c(1, 1), fn)

q50 <- 6.2 #EIR = 100
q025 <- 0.05 #EIR = 0.1
q975 <- 124 #EIR = 2000
optim(c(1, 1), fn)
