library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(didehpc)
library(pkgdepends)
library(tidyverse)
library("coda")
library(binom)
library(ggplot2)
library(bayesplot)
library(zoo)
library(lubridate)

#Required functions
source('shared/run_pmcmc.R')
source('nnp/in_development/run_pmcmc_pg.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')
source('shared/utils.R')

##Get saved processed data
WA_pg_data_list <- readRDS('trial/Data/WA_pg_data_list.rds')
WA_sg_data_list <- readRDS('trial/Data/WA_sg_data_list.rds')
WA_all_data_list <- readRDS('trial/Data/WA_all_data_list.rds')

EA_pg_data_list <- readRDS('trial/Data/EA_pg_data_list.rds')
EA_mg_data_list <- readRDS('trial/Data/EA_mg_data_list.rds')

names(WA_pg_data_list)
admins <- c(`Burkina Faso` = 'Plateau-Central',
            Gambia = 'Upper River',
            Ghana = 'Upper East',
            Mali = 'Koulikoro') ##There were 3 Mali sites in the trial. Chose the 
                                ##admin area in the middle
##test
test_run_pg_wa <- run_pmcmc_pg(data_raw = WA_pg_data_list[[1]],
                               n_particles = 10,
                               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                               max_EIR=1000,
                               max_steps = 1e7,
                               atol = 1e-5,
                               rtol = 1e-6,
                               n_steps = 20,
                               n_threads = 2,
                               lag_rates = 10,
                               country = names(WA_pg_data_list[1]),
                               admin_unit = admins[[names(WA_pg_data_list[1])]],
                               seasonality_on = 1,
                               state_check = 0)
windows(10,7)
plot_particle_filter(test_run_pg_wa$history,true_history=WA_pg_data_list[[1]],times=c(1:nrow(WA_pg_data_list[[1]])))



##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("shared/run_pmcmc.R",
             "nnp/in_development/run_pmcmc_pg.R",
             "nnp/in_development/run_pmcmc_mg.R",
             "shared/model_parameters.R",
             "shared/equilibrium-init-create-stripped.R",
             'shared/utils.R')
packages <- list(loaded = 'plyr',attached = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'))
ctx <- context::context_save("T:/jth/contexts_may23", sources = sources,
                             packages = packages,
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate')))
config_1 <- didehpc::didehpc_config(template = "32Core",cores =4, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
obj <- didehpc::queue_didehpc(ctx,config = config_1)
obj$login()
obj$cluster_load(TRUE)

wa_pg_bulk_seas <- obj$enqueue_bulk(1:4, function(i,data_pg,admins){
  run_pmcmc_pg(data_raw = data_pg[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 4,
               lag_rates = 10,
               seasonality_on = 1,
               country = names(data_pg[i]),
               admin_unit = admins[[names(data_pg[i])]],
               state_check = 0)
},data_pg=WA_pg_data_list,admins=admins)
wa_pg_bulk_seas$status()#'didactic_rainbowfish'
wa_pg_bulk_seas$tasks$`4c0c66f66ddac58dbeb57dc8a754cff5`$log()
obj$unsubmit(wa_pg_bulk_seas$ids)
wa_pg_bulk_seas$tasks$f7ca68b6142eb5a7b89d8effea87e586$log()
wa_pg_bulk_seas$tasks$`821caff3e86e94b8eb66e2c3a1348a41`$log()
wa_pg_seas_result_list <- lapply(1:4, function(id){
  wa_pg_bulk_seas$tasks[[id]]$result()
})

wa_pg_bulk_seas_2 <- obj$enqueue_bulk(1:4, function(i,data_pg,admins){
  run_pmcmc_pg(data_raw = data_pg[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 4,
               lag_rates = 10,
               seasonality_on = 1,
               country = names(data_pg[i]),
               admin_unit = admins[[names(data_pg[i])]],
               state_check = 0,
               seed = 2L)
},data_pg=WA_pg_data_list,admins=admins)
wa_pg_bulk_seas_2 <- obj$task_bundle_get('elder_oryx')
wa_pg_bulk_seas_2$status()#'elder_oryx'
wa_pg_bulk_seas_2$tasks[[3]]$log()

wa_pg_seas_2_result_list <- lapply(1:4, function(id){
  wa_pg_bulk_seas_2$tasks[[id]]$result()
})

admins[[names(WA_pg_data_list[1])]]

wa_pg_bulk_std <- obj$enqueue_bulk(1:4, function(i,data_pg,admins){
  run_pmcmc_pg(data_raw = data_pg[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 4,
               lag_rates = 10,
               seasonality_on = 0,
               country = names(data_pg[i]),
               admin_unit = admins[[names(data_pg[i])]],
               state_check = 0,
               seed = 2L)
},data_pg=WA_pg_data_list,admins=admins)
wa_pg_bulk_std <- obj$task_bundle_get('defunct_wryneck')
wa_pg_bulk_std$status()#''defunct_wryneck''
wa_pg_bulk_std$tasks[[1]]$log()

wa_pg_std_result_list <- lapply(1:4, function(id){
  wa_pg_bulk_std$tasks[[id]]$result()
})


wa_pg_bulk_seas_padded <- obj$enqueue_bulk(1:4, function(i,data_pg,admins){
  run_pmcmc_pg(data_raw = data_pg[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 4,
               lag_rates = 10,
               seasonality_on = 1,
               country = names(data_pg[i]),
               admin_unit = admins[[names(data_pg[i])]],
               state_check = 0,
               seed = 2L,
               start_pf_time = 30*4)
},data_pg=WA_pg_data_list,admins=admins)
wa_pg_bulk_seas_padded<-obj$task_bundle_get('unhumane_hammerheadshark')
wa_pg_bulk_seas_padded$status()#'unhumane_hammerheadshark'
#obj$unsubmit(wa_pg_bulk_seas_padded$ids)
wa_pg_seas_padded_result_list <- lapply(1:4, function(id){
  wa_pg_bulk_seas_padded$tasks[[id]]$result()
})

wa_pg_bulk_seas_padded$tasks[[4]]$log()

data_gamb_all <- rbind(WA_pg_data_list$Gambia,WA_sg_data_list$Gambia)%>%
  group_by(month)%>%
  summarise(tested = sum(tested),
            positive = sum(positive))

wa_pg_gamb_seas <- obj$enqueue(run_pmcmc_pg(data_raw = data_gamb_all,
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 4,
               lag_rates = 10,
               seasonality_on = 1,
               country = 'Gambia',
               admin_unit = 'Upper River',
               state_check = 0,
               seed = 2L,
               start_pf_time = 30))
wa_pg_gamb_seas$id #"4cd1aad5848796687a13e1234594c54f"
wa_pg_gamb_seas <- obj$task_get('4cd1aad5848796687a13e1234594c54f')
wa_pg_gamb_seas$status()
wa_pg_gamb_seas$log()

data_bf_all <- rbind(WA_pg_data_list$`Burkina Faso`,WA_sg_data_list$`Burkina Faso`)%>%
  group_by(month)%>%
  summarise(tested = sum(tested),
            positive = sum(positive))

wa_pg_bfall_seas <- obj$enqueue(run_pmcmc_pg(data_raw = data_bf_all,
                                            n_particles = 200,
                                            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                            max_EIR=1000,
                                            max_steps = 1e7,
                                            atol = 1e-5,
                                            rtol = 1e-6,
                                            n_steps = 250,
                                            n_threads = 4,
                                            lag_rates = 10,
                                            seasonality_on = 1,
                                            seasonality_check = 1,
                                            country = 'Burkina Faso',
                                            admin_unit = 'Plateau-Central',
                                            state_check = 0,
                                            seed = 2L,
                                            start_pf_time = 30*4))
wa_pg_bfall_seas$status()
wa_pg_bfall_seas$log()

wa_pg_bf_seas <- obj$enqueue(run_pmcmc_pg(data_raw = WA_pg_data_list[[1]],
                                             n_particles = 200,
                                             proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                             max_EIR=1000,
                                             max_steps = 1e7,
                                             atol = 1e-5,
                                             rtol = 1e-6,
                                             n_steps = 500,
                                             n_threads = 4,
                                             lag_rates = 10,
                                             seasonality_on = 1,
                                             seasonality_check = 1,
                                             country = 'Burkina Faso',
                                             admin_unit = 'Plateau-Central',
                                             state_check = 0,
                                             seed = 2L,
                                             start_pf_time = 30*4))
wa_pg_bf_seas$status()
wa_pg_bf_seas$times(unit_elapsed = 'mins')
wa_pg_bf_seas$log()

wa_pg_bfall_seas_16 <- obj$enqueue(run_pmcmc_pg(data_raw = data_bf_all,
                                             n_particles = 200,
                                             proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                             max_EIR=1000,
                                             max_steps = 1e7,
                                             atol = 1e-5,
                                             rtol = 1e-6,
                                             n_steps = 250,
                                             n_threads = 16,
                                             lag_rates = 10,
                                             seasonality_on = 1,
                                             seasonality_check = 1,
                                             country = 'Burkina Faso',
                                             admin_unit = 'Plateau-Central',
                                             state_check = 0,
                                             seed = 2L,
                                             start_pf_time = 30*4))
wa_pg_bfall_seas_16$status()
wa_pg_bfall_seas_16$times(unit_elapsed = 'mins')
wa_pg_bfall_seas_16$log()

wa_pg_bf_seas_16 <- obj$enqueue(run_pmcmc_pg(data_raw = WA_pg_data_list[[1]],
                                          n_particles = 200,
                                          proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                          max_EIR=1000,
                                          max_steps = 1e7,
                                          atol = 1e-5,
                                          rtol = 1e-6,
                                          n_steps = 500,
                                          n_threads = 16,
                                          lag_rates = 10,
                                          seasonality_on = 1,
                                          seasonality_check = 1,
                                          country = 'Burkina Faso',
                                          admin_unit = 'Plateau-Central',
                                          state_check = 0,
                                          seed = 2L,
                                          start_pf_time = 30*4))
wa_pg_bf_seas$status()
wa_pg_bf_seas$log()

wa_pg_bfall_seas_pc_2 <- run_pmcmc_pg(data_raw = data_bf_all,
                                             n_particles = 50,
                                             proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                             max_EIR=1000,
                                             max_steps = 1e7,
                                             atol = 1e-5,
                                             rtol = 1e-6,
                                             n_steps = 20,
                                             n_threads = 2,
                                             lag_rates = 10,
                                             seasonality_on = 1,
                                             seasonality_check = 1,
                                             country = 'Burkina Faso',
                                             admin_unit = 'Plateau-Central',
                                             state_check = 0,
                                             seed = 2L,
                                             start_pf_time = 30*4)
wa_pg_bfall_seas_pc_2$history[1,,]
wa_pg_bfall_seas_pc$seas_history
wa_pg_bfall_seas_pc_2$history[1,,]

wa_pg_bf_seas_16_u5prev <- obj$enqueue(run_pmcmc(data_raw = WA_pg_data_list[[1]],
                                             n_particles = 200,
                                             proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                             max_EIR=1000,
                                             max_steps = 1e7,
                                             atol = 1e-5,
                                             rtol = 1e-6,
                                             n_steps = 500,
                                             n_threads = 16,
                                             lag_rates = 10,
                                             seasonality_on = 1,
                                             seasonality_check = 0,
                                             country = 'Burkina Faso',
                                             admin_unit = 'Plateau-Central',
                                             state_check = 0,
                                             seed = 2L))
wa_pg_bf_seas_16_u5prev$status()
wa_pg_bf_seas_16_u5prev$log()

wa_pg_bfall_seas_16_u5prev <- obj$enqueue(run_pmcmc(data_raw = data_bf_all,
                                                 n_particles = 200,
                                                 proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                                 max_EIR=1000,
                                                 max_steps = 1e7,
                                                 atol = 1e-5,
                                                 rtol = 1e-6,
                                                 n_steps = 500,
                                                 n_threads = 16,
                                                 lag_rates = 10,
                                                 seasonality_on = 1,
                                                 seasonality_check = 0,
                                                 country = 'Burkina Faso',
                                                 admin_unit = 'Plateau-Central',
                                                 state_check = 0,
                                                 seed = 2L))
wa_pg_bfall_seas_16_u5prev$status()
wa_pg_bfall_seas_16_u5prev$log()

data_bf_all_standardized <- data_bf_all%>%
  rename(positive_old = positive,
         tested_old = tested)%>%
  mutate(prev=positive_old/tested_old,
         tested = 200,
         positive = round(prev*tested))
data_bf_all_standardized_cis <- addCIs(df=data_bf_all_standardized,Ys=data_bf_all_standardized$positive,data_bf_all_standardized$tested)
wa_pg_bfall_seas_16_sampsz <- obj$enqueue(run_pmcmc_pg(data_raw = data_bf_all_standardized,
                                                       n_particles = 200,
                                                       proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                                       max_EIR=1000,
                                                       max_steps = 1e7,
                                                       atol = 1e-5,
                                                       rtol = 1e-6,
                                                       n_steps = 250,
                                                       n_threads = 16,
                                                       lag_rates = 10,
                                                       seasonality_on = 1,
                                                       seasonality_check = 0,
                                                       country = 'Burkina Faso',
                                                       admin_unit = 'Plateau-Central',
                                                       state_check = 0,
                                                       seed = 2L,
                                                       start_pf_time = 30*4))
wa_pg_bfall_seas_16_sampsz$status()
wa_pg_bfall_seas_16_sampsz$log()
#obj$unsubmit(wa_pg_bfall_seas_16_sampsz$id)
obj$login()
##Check seasonality equilibrium
check_seas_pg_wa <- run_pmcmc_pg(data_raw = WA_pg_data_list[[1]],
                               n_particles = 10,
                               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                               max_EIR=1000,
                               max_steps = 1e7,
                               atol = 1e-5,
                               rtol = 1e-6,
                               n_steps = 20,
                               n_threads = 2,
                               lag_rates = 10,
                               country = names(WA_pg_data_list[1]),
                               admin_unit = admins[[names(WA_pg_data_list[1])]],
                               seasonality_on = 1,
                               seasonality_check = 1)
check_seas_pg_wa$seas_history


####
wa_pg_bulk_seas_090623 <- obj$enqueue_bulk(1:4, function(i,data_pg,admins){
  run_pmcmc_pg(data_raw = data_pg[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 4,
               lag_rates = 10,
               seasonality_on = 1,
               country = names(data_pg[i]),
               admin_unit = admins[[names(data_pg[i])]],
               state_check = 0,
               seed = 2L,
               start_pf_time = 30*4)
},data_pg=WA_pg_data_list,admins=admins)
wa_pg_bulk_seas_090623$status() #''disastrous_andeancondor''
wa_pg_bulk_seas_090623 <- obj$task_bundle_get('disastrous_andeancondor')

wa_pg_bulk_seas_090623$tasks[[1]]$log()

wa_pg_bulk_seas_090623_result_list <- lapply(1:4, function(id){
  wa_pg_bulk_seas_090623$tasks[[id]]$result()
})
wa_pg_bulk_seas_090623_result_list[[1]]$history[1,1,1]
wa_all_bulk_seas_090623 <- obj$enqueue_bulk(1:4, function(i,data_pg,admins){
  run_pmcmc_pg(data_raw = data_pg[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 4,
               lag_rates = 10,
               seasonality_on = 1,
               country = names(data_pg[i]),
               admin_unit = admins[[names(data_pg[i])]],
               state_check = 0,
               seed = 2L,
               start_pf_time = 30*4)
},data_pg=WA_all_data_list,admins=admins)
wa_all_bulk_seas_090623 <- obj$task_bundle_get('citric_snowmonkey')
wa_all_bulk_seas_090623$status() #'citric_snowmonkey'

wa_all_bulk_seas_090623_result_list <- lapply(1:4, function(id){
  wa_all_bulk_seas_090623$tasks[[id]]$result()
})

#### With Sifter #####
##Configure cluster settings
ctx_sifter3 <- context::context_save("T:/jth/contexts.sift3",
                                     packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde','RecordLinkage'),
                                     package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate','jt-hicks/sifter')))
config_1 <- didehpc::didehpc_config(template = "24Core",cores =6, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
# config_dide <- didehpc::didehpc_config(template = "8Core",cores =1, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj_trial <- didehpc::queue_didehpc(ctx_sifter3,config = config_1)
obj_trial$cluster_load(TRUE)
obj_trial$login()

library(foresite)
library(site)
library(malariasimulation)
ghana_file$sites
ghana_file <- get_site('GHA')
ghana_site <- site::single_site(ghana_file,c(13,14))

## get median q1 and 3-rd quartile q2 of p

admins_wa <- c(`Burkina Faso` = 'Plateau-Central',
            Gambia = 'Upper River',
            Ghana = 'Upper East',
            Mali = 'Koulikoro')
country_wa <- names(WA_sg_data_list)
0.583962944
0.127724696
0.681794701
0.297692604
#MAP estimates for 2010
prevs_wa <- c(`Burkina Faso`=0.583962944,Gambia=0.127724696,Ghana=0.681794701,Mali=0.297692604)
init_EIRs <- c(`Burkina Faso`=49.81896,Gambia=1.712536,Ghana=111.0336,Mali=7.419249)

49.81896
1.712536
111.0336
7.419249
admins_wa <- c(`Burkina Faso` = 'Plateau-Central',
               Gambia = 'Upper River',
               Ghana = 'Upper East',
               Mali = 'Koulikoro')
admins_ea <- c(Kenya = 'Siaya',
               Malawi = 'Blantyre') #Two sites are in Blantyre, one is in neighboring Chikwawa
country_ea <- names(EA_pg_data_list)
prevs_ea <- c(0.42,0.391)
wa_pgsg_bulk_sifter_eir <- obj_trial$enqueue_bulk(1:4, function(i,data_pg,data_mg,country,admin){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_EIR=1000,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 6,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 1,
                    stoch_param = 'EIR',
                    comparison = 'pgsg')
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa)
wa_pgsg_bulk_sifter_eir$status() #'theosophic_hogget' - submitted 18 July 5:17pm
wa_pgsg_bulk_sifter_eir$tasks[[1]]$log()
wa_pgsg_bulk_sifter_eir <- obj_trial$task_bundle_get('theosophic_hogget')
wa_pgsg_bulk_sifter_eir$tasks[[4]] <- wa_pgsg_bulk_sifter_eir.4$tasks[[1]]
wa_pgsg_bulk_sifter_eir_results <- lapply(1:4, function(id){
  wa_pgsg_bulk_sifter_eir$tasks[[id]]$result()
})
obj_trial$login()
wa_pgsg_bulk_sifter_eir.4 <- obj_trial$enqueue_bulk(4, function(i,data_pg,data_mg,country,admin){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_param=1000,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 6,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 1,
                    stoch_param = 'EIR',
                    comparison = 'pgsg',
                    seed=2007)
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa)
wa_pgsg_bulk_sifter_eir.4$status() #'predesirous_mantis'
wa_pgsg_bulk_sifter_eir.4$tasks[[1]]$log()
wa_pgsg_bulk_sifter_eir.4 <- obj_trial$task_bundle_get('demoniac_killdeer')
wa_pgsg_bulk_sifter_eir_results[4]<-wa_pgsg_bulk_sifter_eir.4$tasks[[1]]$result()

wa_pgsg_bulk_sifter_eir_max2000 <- obj_trial$enqueue_bulk(1:4, function(i,data_pg,data_mg,country,admin){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_param=2000,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 6,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 1,
                    stoch_param = 'EIR',
                    comparison = 'pgsg')
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa)
wa_pgsg_bulk_sifter_eir_max2000 <- obj_trial$task_bundle_get('enchanting_cob')
wa_pgsg_bulk_sifter_eir_max2000$status() #'enchanting_cob' submitted 21Jul 10:14am
wa_pgsg_bulk_sifter_eir_max2000$tasks[[3]]$log()

obj_trial$login()
wa_pgsg_bulk_sifter_eir_max2000.3 <- obj_trial$enqueue_bulk(3, function(i,data_pg,data_mg,country,admin){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_param=2000,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 6,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 1,
                    stoch_param = 'EIR',
                    comparison = 'pgsg',
                    seed = 2407)
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa)
wa_pgsg_bulk_sifter_eir_max2000.3 <- obj_trial$task_bundle_get('declinate_daddylonglegs')
wa_pgsg_bulk_sifter_eir_max2000.3$status() #'declinate_daddylonglegs' submitted 24Jul 9:54am
# obj_trial$unsubmit(wa_pgsg_bulk_sifter_eir_max2000.3$ids)

##betaa##

wa_pgsg_bulk_sifter_betaa <- obj_trial$enqueue_bulk(1:4, function(i,data_pg,data_mg,country,admin){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_param=62,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 6,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 1,
                    stoch_param = 'betaa',
                    comparison = 'pgsg')
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa)
wa_pgsg_bulk_sifter_betaa$status() #'weatherbeaten_antarcticfurseal' - submitted 18 July 5:19pm
wa_pgsg_bulk_sifter_betaa$status() #'mushy_globefish' - submitted 20 July 3:25pm
wa_pgsg_bulk_sifter_betaa$tasks[[1]]$log()
wa_pgsg_bulk_sifter_betaa <- obj_trial$task_bundle_get('weatherbeaten_antarcticfurseal')
wa_pgsg_bulk_sifter_betaa_results <- lapply(1:4, function(id){
  wa_pgsg_bulk_sifter_betaa$tasks[[id]]$result()
})

obj_trial$login()
wa_pgsg_bulk_sifter_betaa_max125 <- obj_trial$enqueue_bulk(1:4, function(i,data_pg,data_mg,country,admin){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 6,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 1,
                    stoch_param = 'betaa',
                    comparison = 'pgsg')
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa)
wa_pgsg_bulk_sifter_betaa_max125 <- obj_trial$task_bundle_get('preoccupied_basil')
wa_pgsg_bulk_sifter_betaa_max125$status() #preoccupied_basil submitted 21Jul 10:11am
wa_pgsg_bulk_sifter_betaa_max125_results <- lapply(1:4, function(id){
  wa_pgsg_bulk_sifter_betaa_max125$tasks[[id]]$result()
})
names(wa_pgsg_bulk_sifter_betaa_max125_results) <- names(WA_pg_data_list)
proposals <- lapply(c(1:4), function(i) cov(wa_pgsg_bulk_sifter_betaa_max125_results[[i]]$pars))
wa_pgsg_bulk_sifter_betaa_2407 <- obj_trial$enqueue_bulk(1:4, function(i,data_pg,data_mg,country,admin, proposal_matrix){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = proposal_matrix,
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 6,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 1,
                    stoch_param = 'betaa',
                    comparison = 'pgsg',
                    seed = 2407)
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa,proposal_matrix=proposals[[i]])
wa_pgsg_bulk_sifter_betaa_2407$status() #lawabiding_moa submitted 25Jul 3:40pm
wa_pgsg_bulk_sifter_betaa_2407$tasks[[4]]$log()
wa_pgsg_bulk_sifter_betaa_2407 <- obj_trial$task_bundle_get('lawabiding_moa')
wa_pgsg_bulk_sifter_betaa_2407$tasks[[3]] <- wa_pgsg_bulk_sifter_betaa_2407.3$tasks[[1]]
wa_pgsg_bulk_sifter_betaa_2407_results <- lapply(1:4, function(i){
  wa_pgsg_bulk_sifter_betaa_2407$tasks[[i]]$result()
})
names(wa_pgsg_bulk_sifter_betaa_2407_results) <- names(WA_pg_data_list)
wa_proposals_2407 <- lapply(c(1:4), function(i) cov(wa_pgsg_bulk_sifter_betaa_2407_results[[i]]$pars))

wa_pgsg_bulk_sifter_betaa_2407.3 <- obj_trial$enqueue_bulk(3, function(i,data_pg,data_mg,country,admin, proposal_matrix){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = proposal_matrix,
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 6,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 1,
                    stoch_param = 'betaa',
                    comparison = 'pgsg',
                    seed = 2607)
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa,proposal_matrix=proposals[[i]])
wa_pgsg_bulk_sifter_betaa_2407.3$status() #rotatable_monarch submitted 26Jul 10:13am
wa_pgsg_bulk_sifter_betaa_2407.3 <- obj_trial$task_bundle_get('rotatable_monarch')

##W Kenya and Malawi
#EIR
obj_trial$login()
ea_pgmg_bulk_sifter_eir_max2000 <- obj_trial$enqueue_bulk(1:2, function(i,data_pg,data_mg,country,admin){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_param=2000,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 6,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 1,
                    stoch_param = 'EIR',
                    comparison = 'pgmg')
},data_pg=EA_pg_data_list,data_mg=EA_mg_data_list,country=country_ea,admin=admins_ea)
# obj_trial$unsubmit(ea_pgmg_bulk_sifter_eir_max2000$ids)
ea_pgmg_bulk_sifter_eir_max2000 <- obj_trial$task_bundle_get('astraphobic_mutt')
ea_pgmg_bulk_sifter_eir_max2000$status() #'astraphobic_mutt' submitted 21 July 11:21am
ea_pgmg_bulk_sifter_eir_max2000_results <- lapply(1:2, function(id){
  ea_pgmg_bulk_sifter_eir_max2000$tasks[[id]]$result()
})
names(ea_pgmg_bulk_sifter_eir_max2000_results) <- names(EA_pg_data_list)

#betaa
ea_pgmg_bulk_sifter_betaa_max125 <- obj_trial$enqueue_bulk(1:2, function(i,data_pg,data_mg,country,admin){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 6,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 1,
                    stoch_param = 'betaa',
                    comparison = 'pgmg')
},data_pg=EA_pg_data_list,data_mg=EA_mg_data_list,country=country_ea,admin=admins_ea)
ea_pgmg_bulk_sifter_betaa_max125 <- obj_trial$task_bundle_get('uncongestive_duckling')
ea_pgmg_bulk_sifter_betaa_max125$status() #'uncongestive_duckling' Submitted 21 July 11:23am
ea_pgmg_bulk_sifter_betaa_max125_results <- lapply(1:2, function(id){
  ea_pgmg_bulk_sifter_betaa_max125$tasks[[id]]$result()
})
names(ea_pgmg_bulk_sifter_betaa_max125_results) <- names(EA_pg_data_list)

ea_pgmg_bulk_sifter_betaa_max125_ky <- obj_trial$enqueue_bulk(1, function(i,data_pg,data_mg,country,admin){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 6,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 1,
                    stoch_param = 'betaa',
                    comparison = 'pgmg')
},data_pg=EA_pg_data_list,data_mg=EA_mg_data_list,country=country_ea,admin=admins_ea)
ea_pgmg_bulk_sifter_betaa_max125_ky$status() #'unmystified_olingo' submitted 26 Jul 10:55am

proposals_ea <- lapply(c(1:2), function(i) cov(ea_pgmg_bulk_sifter_betaa_max125_results[[i]]$pars))
proposals_ky <- cov(ea_pgmg_bulk_sifter_betaa_max125_ky$tasks[[1]]$result()$pars)
ea_pgmg_bulk_sifter_betaa_2407 <- obj_trial$enqueue_bulk(1:2, function(i,data_pg,data_mg,country,admin, proposal_matrix){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = proposal_matrix[[i]],
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 6,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 1,
                    stoch_param = 'betaa',
                    comparison = 'pgmg',
                    seed = 2407)
},data_pg=EA_pg_data_list,data_mg=EA_mg_data_list,country=country_ea,admin=admins_ea,proposal_matrix=proposals_ea)
ea_pgmg_bulk_sifter_betaa_2407 <- obj_trial$task_bundle_get('ferrety_bird')
ea_pgmg_bulk_sifter_betaa_2407$status() #'ferrety_bird' submitted 24 July 3:48pm
ea_pgmg_bulk_sifter_betaa_2407$tasks[[1]]$log()
ea_pgmg_bulk_sifter_betaa_2407$tasks[[1]] <- ea_pgmg_bulk_sifter_betaa_2407_ky$tasks[[1]]
ea_pgmg_bulk_sifter_betaa_2407_results <- lapply(1:2, function(id){
  ea_pgmg_bulk_sifter_betaa_2407$tasks[[id]]$result()
})
names(ea_pgmg_bulk_sifter_betaa_2407_results) <- names(EA_pg_data_list)
ea_proposals_2407 <- lapply(c(1:2), function(i) cov(ea_pgmg_bulk_sifter_betaa_2407_results[[i]]$pars))

ea_pgmg_bulk_sifter_betaa_2407_ky <- obj_trial$enqueue_bulk(1, function(i,data_pg,data_mg,country,admin, proposal_matrix){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = proposal_matrix,
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 6,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 1,
                    stoch_param = 'betaa',
                    comparison = 'pgmg',
                    seed = 2407)
},data_pg=EA_pg_data_list,data_mg=EA_mg_data_list,country=country_ea,admin=admins_ea,proposal_matrix=proposals_ky)
ea_pgmg_bulk_sifter_betaa_2407_ky$status() #apiarian_corydorascatfish 27 July 8:56am
ea_pgmg_bulk_sifter_betaa_2407_ky <- obj_trial$task_bundle_get('apiarian_corydorascatfish')

##Initial Runs##
##using obj_init - see nnp_mg_submit.R for initialization
wa_proposals_2407
obj_init$login()
test_wa_pgsg_bulk_sifter_betaa <- obj_init$enqueue_bulk(2, function(i,data_pg,data_mg,country,admin, proposal_matrix){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = proposal_matrix,
                    target_prev = 0.73,
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 100,
                    n_threads = 8,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'pgsg',
                    seed = 1008)
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa,proposal_matrix=proposal_mat)
test_wa_pgsg_bulk_sifter_betaa$status() #augmented_halcyon
test_wa_pgsg_bulk_sifter_betaa$tasks[[1]]$log()
obj_init$unsubmit(test_wa_pgsg_bulk_sifter_betaa$ids)
test_result <- test_wa_pgsg_bulk_sifter_betaa$tasks[[1]]$result()
matplot(test_result$times[-1],t(test_result$history['betaa',,-1]),col='black',type='l')
ghana_data <- WA_pg_data_list[[3]]
ghana_data$date <- as.Date(ghana_data$month)
start_obs <- min(zoo::as.Date(zoo::as.yearmon((ghana_data$month))))#Month of first observation (in Date format)
time_origin <- as.Date(paste0(year(start_obs)-1,'-01-01')) #January 1 of year before observation (in Date format)
ghana_data$date <- zoo::as.Date(zoo::as.yearmon(ghana_data$month), frac = 0.5) #Convert dates to middle of month
ghana_data$t <- as.integer(difftime(ghana_data$date,time_origin,units="days"))

matplot(test_result$times[-1],t(test_result$history['prev_pg',,-1]),col='black',type='l',ylim = c(0,1))

lines(ghana_data$t,WA_pg_data_list[[3]]$positive/WA_pg_data_list[[3]]$tested,type='p',col='red')
smc_times <- ghana_data[month(ghana_data$date)%in%c(7,8,9,10,11),]$t
ghana_betaa <- bind_cols(c(data.frame(t=test_result$times[-1],
                                      date=ghana_data$date),
                           data.frame(t(test_result$history['betaa',,-1]))))
ghana_eir <- bind_cols(c(data.frame(t=test_result$times[-1],
                                      date=ghana_data$date),
                           data.frame(t(test_result$history['EIR',,-1]))))
saveRDS(smc_times,'./trial/ghana_smc_times.rds')
saveRDS(ghana_betaa,'./trial/ghana_betaa_test.rds')
saveRDS(ghana_eir,'./trial/ghana_eir_test.rds')

wa_pgsg_bulk_sifter_betaa_0109 <- obj_init$enqueue_bulk(1:4, function(i,data_pg,data_mg,country,admin, proposal_matrix,baseline_prevs){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = proposal_matrix,
                    target_prev = baseline_prevs[i],
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*1,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'pgsg',
                    seed = 1008)
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa,proposal_matrix=proposal_mat,baseline_prevs=prevs_wa)
wa_pgsg_bulk_sifter_betaa_0109$status() #unsecular_andeancondor submitted 3:39pm on 1 Sep
wa_pgsg_bulk_sifter_betaa_0109$tasks[[4]]$log()
c(`Burkina Faso`=278.5669,Gambia=1.024789,Ghana=8.370505,Mali=37.63816)
wa_pgsg_bulk_sifter_betaa_0109_results <- lapply(1:4, function(i){
  wa_pgsg_bulk_sifter_betaa_0109$tasks[[i]]$result()
})
names(wa_pgsg_bulk_sifter_betaa_0109_results) <- names(WA_pg_data_list)
wa_proposals_0109 <- lapply(c(1:4), function(i) cov(wa_pgsg_bulk_sifter_betaa_0109_results[[i]]$pars))
diag_wa_pgsg_bulk_sifter_betaa_0109 <- lapply(1:length(wa_pgsg_bulk_sifter_betaa_0109_results),function(i) create_diag_figs(wa_pgsg_bulk_sifter_betaa_0109_results[[i]],country = names(wa_pgsg_bulk_sifter_betaa_0109_results)[[i]]))
for(i in c(1:length(diag_wa_pgsg_bulk_sifter_betaa_0109))){
  ggsave(paste0('Q:/anc_pmcmc/trial/figures/Diagnostic/',names(wa_pgsg_bulk_sifter_betaa_0109_results)[[i]],'_wa_pgsg_bulk_sifter_betaa_0109.pdf'),plot=diag_wa_pgsg_bulk_sifter_betaa_0109[[i]], width = 10, height = 7)
}
obj_init$login()
wa_pgsg_bulk_sifter_betaa_0309 <- obj_init$enqueue_bulk(1:4, function(i,data_pg,data_mg,country,admin, proposal_matrix,baseline_prevs){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = proposal_matrix[[i]],
                    target_prev = baseline_prevs[i],
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*1,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'pgsg',
                    seed = 1008)
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa,proposal_matrix=wa_proposals_0109,baseline_prevs=prevs_wa)
wa_pgsg_bulk_sifter_betaa_0309$status() #blotchy_americancrow submitted 4:36pm on 3 Sep
wa_pgsg_bulk_sifter_betaa_0309_results <- lapply(1:4, function(i){
  wa_pgsg_bulk_sifter_betaa_0309$tasks[[i]]$result()
})
names(wa_pgsg_bulk_sifter_betaa_0309_results) <- names(WA_pg_data_list)
diag_wa_pgsg_bulk_sifter_betaa_0309 <- lapply(1:length(wa_pgsg_bulk_sifter_betaa_0309_results),function(i) create_diag_figs(wa_pgsg_bulk_sifter_betaa_0309_results[[i]],country = names(wa_pgsg_bulk_sifter_betaa_0309_results)[[i]]))
for(i in c(1:length(diag_wa_pgsg_bulk_sifter_betaa_0309))){
  ggsave(paste0('Q:/anc_pmcmc/trial/figures/Diagnostic/',names(wa_pgsg_bulk_sifter_betaa_0309_results)[[i]],'_wa_pgsg_bulk_sifter_betaa_0309.pdf'),plot=diag_wa_pgsg_bulk_sifter_betaa_0309[[i]], width = 10, height = 7)
}
create_dashboard_plots_trial(results=wa_pgsg_bulk_sifter_betaa_0309_results)
matplot(WA_pg_data_list[[3]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_results[[3]]$history['prev_05',,-1]),col='black',type='l',ylim = c(0,1))
lines(WA_pg_data_list[[3]]$month,WA_pg_data_list[[3]]$positive/WA_pg_data_list[[3]]$tested,type='l',col='red')
matplot(WA_sg_data_list[[3]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_results[[3]]$history['prev_sg',,-1]),col='black',type='l',ylim = c(0,1))
lines(WA_sg_data_list[[3]]$month,WA_sg_data_list[[3]]$positive/WA_sg_data_list[[3]]$tested,type='l',col='red')
matplot(WA_sg_data_list[[3]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_results[[3]]$history['clininc_all',,-1]),col='black',type='l')
matplot(WA_sg_data_list[[3]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_results[[3]]$history['betaa',,-1]),col='black',type='l')
ghana_betaa <- bind_cols(c(data.frame(t=test_result$times[-1],
                                      date=ghana_data$date),
                           data.frame(t(test_result$history['betaa',,-1]))))
ghana_eir <- bind_cols(c(data.frame(t=test_result$times[-1],
                                    date=ghana_data$date),
                         data.frame(t(test_result$history['EIR',,-1]))))

wa_pgsg_bulk_sifter_betaa_0309_dash <- create_dashboard_plots_trial(results = wa_pgsg_bulk_sifter_betaa_0309_results,
                                                                      prev_pg=WA_pg_data_list,
                                                                      prev_mg=WA_sg_data_list,
                                                                      coefs_pg_df = as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS')),
                                                                      coefs_mg_df = as.data.frame(readRDS('./nnp/Corr/sg_corr_sample.RDS')),
                                                                      start_pf_time = 30*4,
                                                                      param = 'betaa')
windows(10,6)
wa_pgsg_bulk_sifter_betaa_0309_dash$full_dash
wa_pgsg_bulk_sifter_betaa_0309_flat <- obj_init$enqueue_bulk(1:4, function(i,data_pg,data_mg,country,admin, proposal_matrix,baseline_prevs){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = proposal_matrix[[i]],
                    target_prev = baseline_prevs[i],
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 8,
                    lag_rates = 10,
                    seasonality_on = 0,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'pgsg',
                    seed = 1008)
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa,proposal_matrix=wa_proposals_0109,baseline_prevs=prevs_wa)
wa_pgsg_bulk_sifter_betaa_0309_flat$status() #semisolemn_bluefintuna submitted 10:11am on 4 Sep
wa_pgsg_bulk_sifter_betaa_0309_flat$times()
wa_pgsg_bulk_sifter_betaa_0309_flat_results <- lapply(c(1,2,4), function(i){
  wa_pgsg_bulk_sifter_betaa_0309_flat$tasks[[i]]$result()
})
names(wa_pgsg_bulk_sifter_betaa_0309_flat_results) <- names(WA_pg_data_list)[c(1,2,4)]
diag_wa_pgsg_bulk_sifter_betaa_0309_flat <- lapply(1:length(wa_pgsg_bulk_sifter_betaa_0309_flat_results),function(i) create_diag_figs(wa_pgsg_bulk_sifter_betaa_0309_flat_results[[i]],country = names(wa_pgsg_bulk_sifter_betaa_0309_flat_results)[[i]]))
for(i in c(1:length(diag_wa_pgsg_bulk_sifter_betaa_0309_flat))){
  ggsave(paste0('Q:/anc_pmcmc/trial/figures/Diagnostic/',names(wa_pgsg_bulk_sifter_betaa_0309_flat_results)[[i]],'_wa_pgsg_bulk_sifter_betaa_0309_flat.pdf'),plot=diag_wa_pgsg_bulk_sifter_betaa_0309_flat[[i]], width = 10, height = 7)
}
create_dashboard_plots_trial(results=wa_pgsg_bulk_sifter_betaa_0309_results)
matplot(WA_pg_data_list[[1]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_flat_results[[1]]$history['prev_05',,-1]),col='black',type='l',ylim = c(0,1))
lines(WA_pg_data_list[[1]]$month,WA_pg_data_list[[1]]$positive/WA_pg_data_list[[1]]$tested,type='l',col='red')
matplot(WA_sg_data_list[[1]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_flat_results[[1]]$history['prev_sg',,-1]),col='black',type='l',ylim = c(0,1))
lines(WA_sg_data_list[[1]]$month,WA_sg_data_list[[1]]$positive/WA_sg_data_list[[1]]$tested,type='l',col='red')

wa_pgsg_bulk_sifter_betaa_0309_flat2 <- obj_init$enqueue_bulk(1:4, function(i,data_pg,data_mg,country,admin, proposal_matrix,baseline_prevs){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = proposal_matrix[[i]],
                    target_prev = baseline_prevs[i],
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 8,
                    lag_rates = 10,
                    seasonality_on = 0,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'pgsg',
                    seed = 1008)
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa,proposal_matrix=wa_proposals_0109,baseline_prevs=prevs_wa)
wa_pgsg_bulk_sifter_betaa_0309_flat2$status() #squirarchal_walkingstick submitted at 2:41pm on 4 Sept
wa_pgsg_bulk_sifter_betaa_0309_flat2$tasks[[4]]$log()

wa_pgsg_bulk_sifter_betaa_0309_flat3 <- obj_init$enqueue_bulk(1:4, function(i,data_pg,data_mg,country,admin, proposal_matrix,baseline_prevs){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 96,
                    proposal_matrix = proposal_matrix[[i]],
                    target_prev = baseline_prevs[i],
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 8,
                    lag_rates = 10,
                    seasonality_on = 0,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'pgsg',
                    seed = 1008)
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa,proposal_matrix=wa_proposals_0109,baseline_prevs=prevs_wa)
wa_pgsg_bulk_sifter_betaa_0309_flat3 <- obj_init$task_bundle_get('displayed_frillneckedlizard')
wa_pgsg_bulk_sifter_betaa_0309_flat3$status() #displayed_frillneckedlizard
wa_pgsg_bulk_sifter_betaa_0309_flat3$times()
wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[2]]
wa_pgsg_bulk_sifter_betaa_0309_flat3_results <- lapply(1:4, function(i){
  wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[i]]$result()
})
names(wa_pgsg_bulk_sifter_betaa_0309_flat3_results) <- names(WA_pg_data_list)
wa_pgsg_bulk_sifter_betaa_0309_flat3_dash <- create_dashboard_plots_trial(results=wa_pgsg_bulk_sifter_betaa_0309_flat3_results,
                                                               prev_pg=WA_pg_data_list,
                                                               prev_mg=WA_sg_data_list,
                                                               coefs_pg_df = as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS')),
                                                               coefs_mg_df = as.data.frame(readRDS('./nnp/Corr/sg_corr_sample.RDS')),
                                                               start_pf_time = 30*4,
                                                               param = 'betaa')
windows(10,6)
wa_pgsg_bulk_sifter_betaa_0309_flat3_dash$full_dash
wa_pgsg_bulk_sifter_betaa_0309_flat3_dash$obs_data_dash
wa_pgsg_bulk_sifter_betaa_0309_flat3_dash$pres_dash
wa_pgsg_bulk_sifter_betaa_0309_flat3_dash$std_dash
wa_pgsg_bulk_sifter_betaa_0309_flat3_dash$trans_dash
wa_pgsg_bulk_sifter_betaa_0309_flat3_dash$overlay_plot

gambia_flat3 <- create_diag_figs(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[2]]$result(),country = 'Gambia')
mali_flat3 <- create_diag_figs(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[4]]$result(),country = 'Mali')
bf_flat3 <- create_diag_figs(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[1]]$result(),country = 'Mali')
ghana_flat3 <- create_diag_figs(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[3]]$result(),country = 'Ghana')

matplot(WA_pg_data_list[[1]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[1]]$result()$history['prev_05',,-1]),col='black',type='l',ylim = c(0,1))
lines(WA_pg_data_list[[1]]$month,WA_pg_data_list[[1]]$positive/WA_pg_data_list[[1]]$tested,type='l',col='red')
matplot(WA_sg_data_list[[1]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[1]]$result()$history['prev_sg',,-1]),col='black',type='l',ylim = c(0,1))
lines(WA_sg_data_list[[1]]$month,WA_sg_data_list[[1]]$positive/WA_sg_data_list[[1]]$tested,type='l',col='red')

matplot(WA_pg_data_list[[2]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[2]]$result()$history['prev_05',,-1]),col='black',type='l',ylim = c(0,1))
lines(WA_pg_data_list[[2]]$month,WA_pg_data_list[[2]]$positive/WA_pg_data_list[[2]]$tested,type='l',col='red')
matplot(WA_sg_data_list[[2]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[2]]$result()$history['prev_sg',,-1]),col='black',type='l',ylim = c(0,1))
lines(WA_sg_data_list[[2]]$month,WA_sg_data_list[[2]]$positive/WA_sg_data_list[[2]]$tested,type='l',col='red')

matplot(WA_pg_data_list[[3]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[3]]$result()$history['prev_05',,-1]),col='black',type='l',ylim = c(0,1))
lines(WA_pg_data_list[[3]]$month,WA_pg_data_list[[3]]$positive/WA_pg_data_list[[3]]$tested,type='l',col='red')
matplot(WA_sg_data_list[[3]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[3]]$result()$history['prev_sg',,-1]),col='black',type='l',ylim = c(0,1))
lines(WA_sg_data_list[[3]]$month,WA_sg_data_list[[3]]$positive/WA_sg_data_list[[3]]$tested,type='l',col='red')

matplot(WA_pg_data_list[[4]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[4]]$result()$history['prev_05',,-1]),col='black',type='l',ylim = c(0,1))
lines(WA_pg_data_list[[4]]$month,WA_pg_data_list[[4]]$positive/WA_pg_data_list[[4]]$tested,type='l',col='red')
matplot(WA_sg_data_list[[4]]$month,t(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[4]]$result()$history['prev_sg',,-1]),col='black',type='l',ylim = c(0,1))
lines(WA_sg_data_list[[4]]$month,WA_sg_data_list[[4]]$positive/WA_sg_data_list[[4]]$tested,type='l',col='red')
wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[4]]$result()$times

mcmc_sample <- sample(c(201:1000), 100)
start_stoch <- zoo::as.Date(min(WA_pg_data_list[[2]]$month) - start_pf_time) #Start of stochastic schedule; needs to start when particle filter starts


ghana_data <- WA_pg_data_list[[3]]
ghana_data$date <- as.Date(ghana_data$month)
start_obs <- min(zoo::as.Date(zoo::as.yearmon((ghana_data$month))))#Month of first observation (in Date format)
time_origin <- as.Date(paste0(year(start_obs)-1,'-01-01')) #January 1 of year before observation (in Date format)
ghana_data$date <- zoo::as.Date(zoo::as.yearmon(ghana_data$month), frac = 0.5) #Convert dates to middle of month
ghana_data$t <- as.integer(difftime(ghana_data$date,time_origin,units="days"))
ghana_betaa <- bind_cols(c(data.frame(t=wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[3]]$result()$times[-1],
                                       date=ghana_data$date),
                            data.frame(t(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[3]]$result()$history['betaa',mcmc_sample,-1]))))
ghana_smc_times <- ghana_data[month(ghana_data$date)%in%c(7,8,9,10,11),]$t
saveRDS(ghana_smc_times,'./trial/ghana_smc_times.rds')
saveRDS(ghana_betaa,'./trial/ghana_betaa_flat3.rds')
WA_sg_data_list[[2]]
start_obs <- min(zoo::as.Date(zoo::as.yearmon(WA_pg_data_list[[2]]$month)))#Month of first observation (in Date format)
# cat('start_obs: ',as.character(start_obs),'\n')
time_origin <- zoo::as.Date(paste0(lubridate::year(start_obs)-1,'-01-01')) #January 1 of year before observation (in Date format)
# cat('time_origin: ',as.character(time_origin),'\n')
data_raw_time <- WA_pg_data_list[[2]]
data_raw_time$date <- zoo::as.Date(zoo::as.yearmon(WA_pg_data_list[[2]]$month), frac = 0.5) #Convert dates to middle of month
data_raw_time$t <- as.integer(difftime(data_raw_time$date,time_origin,units="days")) #Calculate date as number of days since January 1 of year before observation
# cat('First observed time: ',min(data_raw_time$t),'\n')
initial_time <- min(data_raw_time$t) - start_pf_time #Start particle filter a given time (default = 30d) before first observation
start_stoch <- zoo::as.Date(min(data_raw_time$date) - start_pf_time) #Start of stochastic schedule; needs to start when particle filter starts
data <- mcstate::particle_filter_data(data_raw_time, time = "t", rate = NULL, initial_time = initial_time) #Declares data to be used for particle filter fitting
stoch_update_dates <- seq.Date(start_stoch,max(as.Date(data_raw_time$date+30),na.rm = TRUE),by='month')
stochastic_schedule <- as.integer(difftime(stoch_update_dates,time_origin,units="days"))

gambia_data <- WA_pg_data_list[[2]]
gambia_data$date <- as.Date(gambia_data$month)
start_obs <- min(zoo::as.Date(zoo::as.yearmon((gambia_data$month))))#Month of first observation (in Date format)
time_origin <- as.Date(paste0(year(start_obs)-1,'-01-01')) #January 1 of year before observation (in Date format)
gambia_data$date <- zoo::as.Date(zoo::as.yearmon(gambia_data$month), frac = 0.5) #Convert dates to middle of month
gambia_data$t <- as.integer(difftime(gambia_data$date,time_origin,units="days"))
gambia_betaa <- bind_cols(c(data.frame(t=wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[2]]$result()$times[-1],
                                      date=gambia_data$date),
                           data.frame(t(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[2]]$result()$history['betaa',mcmc_sample,-1]))))

gambia_smc_times <- gambia_data[month(gambia_data$date)%in%c(9,10,11),]$t
saveRDS(gambia_smc_times,'./trial/gambia_smc_times.rds')
saveRDS(gambia_betaa,'./trial/gambia_betaa_flat3.rds')

mali_data <- WA_pg_data_list[[4]]
mali_data$date <- as.Date(mali_data$month)
start_obs <- min(zoo::as.Date(zoo::as.yearmon((mali_data$month))))#Month of first observation (in Date format)
time_origin <- as.Date(paste0(year(start_obs)-1,'-01-01')) #January 1 of year before observation (in Date format)
mali_data$date <- zoo::as.Date(zoo::as.yearmon(mali_data$month), frac = 0.5) #Convert dates to middle of month
mali_data$t <- as.integer(difftime(mali_data$date,time_origin,units="days"))
mali_betaa <- bind_cols(c(data.frame(t=wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[4]]$result()$times[-1],
                                       date=mali_data$date),
                            data.frame(t(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[4]]$result()$history['betaa',mcmc_sample,-1]))))
mali_smc_times <- mali_data[month(mali_data$date)%in%c(8,9,10),]$t
saveRDS(mali_smc_times,'./trial/mali_smc_times.rds')
saveRDS(mali_betaa,'./trial/mali_betaa_flat3.rds')

bf_data <- WA_pg_data_list[[1]]
bf_data$date <- as.Date(bf_data$month)
start_obs <- min(zoo::as.Date(zoo::as.yearmon((bf_data$month))))#Month of first observation (in Date format)
time_origin <- as.Date(paste0(year(start_obs)-1,'-01-01')) #January 1 of year before observation (in Date format)
bf_data$date <- zoo::as.Date(zoo::as.yearmon(bf_data$month), frac = 0.5) #Convert dates to middle of month
bf_data$t <- as.integer(difftime(bf_data$date,time_origin,units="days"))
bf_betaa <- bind_cols(c(data.frame(t=wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[1]]$result()$times[-1],
                                     date=bf_data$date),
                          data.frame(t(wa_pgsg_bulk_sifter_betaa_0309_flat3$tasks[[1]]$result()$history['betaa',mcmc_sample,-1]))))
bf_smc_times <- bf_data[month(bf_data$date)%in%c(8,9,10),]$t
saveRDS(bf_smc_times,'./trial/bf_smc_times.rds')
saveRDS(bf_betaa,'./trial/bf_betaa_flat3.rds')

names(prevs_wa)
# obj_init$unsubmit(wa_pgsg_bulk_sifter_betaa_0109$ids)
WA_sg_data_list[[2]]
ea_pgsg_bulk_sifter_betaa_0109 <- obj_init$enqueue_bulk(1:2, function(i,data_pg,data_mg,country,admin, proposal_matrix,baseline_prevs){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = proposal_matrix,
                    target_prev = baseline_prevs[i],
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*1,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'pgmg',
                    seed = 1008)
},data_pg=EA_pg_data_list,data_mg=EA_mg_data_list,country=country_ea,admin=admins_ea,proposal_matrix=proposal_mat,baseline_prevs=prevs_ea)
ea_pgsg_bulk_sifter_betaa_0109$status() #repealable_amethystinepython submitted 1 Sep 3:39pm
ea_pgsg_bulk_sifter_betaa_0109_results <- lapply(1:2, function(i){
  ea_pgsg_bulk_sifter_betaa_0109$tasks[[i]]$result()
})
names(ea_pgsg_bulk_sifter_betaa_0109_results) <- names(EA_pg_data_list)
ea_proposals_0109 <- lapply(c(1:2), function(i) cov(ea_pgsg_bulk_sifter_betaa_0109_results[[i]]$pars))
ea_pgsg_bulk_sifter_betaa_0309 <- obj_init$enqueue_bulk(1:2, function(i,data_pg,data_mg,country,admin, proposal_matrix,baseline_prevs){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = proposal_matrix[[i]],
                    target_prev = baseline_prevs[i],
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*1,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'pgmg',
                    seed = 1008)
},data_pg=EA_pg_data_list,data_mg=EA_mg_data_list,country=country_ea,admin=admins_ea,proposal_matrix=proposal_mat,baseline_prevs=ea_proposals_0109)
ea_pgsg_bulk_sifter_betaa_0309$status() #sweating_easternnewt submitted 3 Sep 4:39pm


##Long Runs##
##Configure cluster settings
ctx_long <- context::context_save("T:/jth/contexts.long",
                                  packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde','RecordLinkage'),
                                  package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate','jt-hicks/sifter')))
config_long <- didehpc::didehpc_config(template = "16Core",cores =8, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
# config_dide <- didehpc::didehpc_config(template = "8Core",cores =1, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj_long <- didehpc::queue_didehpc(ctx_long,config = config_long)
obj_long$cluster_load(TRUE)

wa_pgsg_bulk_sifter_betaa_1008 <- obj_long$enqueue_bulk(1:4, function(i,data_pg,data_mg,country,admin, proposal_matrix){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = proposal_matrix[[i]],
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 10000,
                    n_threads = 8,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'pgsg',
                    seed = 1008)
},data_pg=WA_pg_data_list,data_mg=WA_sg_data_list,country=country_wa,admin=admins_wa,proposal_matrix=wa_proposals_2407)
wa_pgsg_bulk_sifter_betaa_1008$status() #scribblenautable_pterosaurs submitted 10 Aug 2:45pm

ea_pgsg_bulk_sifter_betaa_1008 <- obj_long$enqueue_bulk(1:2, function(i,data_pg,data_mg,country,admin, proposal_matrix){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = proposal_matrix[[i]],
                    max_param=125,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 10000,
                    n_threads = 8,
                    lag_rates = 10,
                    seasonality_on = 1,
                    country = country[i],
                    admin_unit = admin[i],
                    state_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'pgmg',
                    seed = 1008)
},data_pg=EA_pg_data_list,data_mg=EA_mg_data_list,country=country_ea,admin=admins_ea,proposal_matrix=ea_proposals_2407)
ea_pgsg_bulk_sifter_betaa_1008$status() #scribblenautable_pterosaurs submitted 10 Aug 2:45pm
