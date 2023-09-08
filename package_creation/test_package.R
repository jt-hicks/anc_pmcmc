library('devtools')
##Install these from github, if needed
remotes::install_github("mrc-ide/mcstate", upgrade = TRUE)
remotes::install_github("mrc-ide/dust", upgrade = TRUE)
remotes::install_github("mrc-ide/mode", upgrade = TRUE)
remotes::install_github("mrc-ide/odin.dust", upgrade = TRUE)
remotes::install_github("mrc-ide/didehpc", upgrade = TRUE)
#devtools::install_github('klutometis/roxygen')
#library(roxygen2)

## Load package from github
devtools::install_github('jt-hicks/sifter@issue-5',force = TRUE)
library(sifter)

test_run_sifter <- sifter::run_pmcmc(data_raw = sifter::data_sim, #I've added data_sim to the package for an easy test
                                     n_particles = 200,
                                     proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                     max_EIR=1000,
                                     max_steps = 1e7,
                                     atol = 1e-6,
                                     rtol = 1e-6,
                                     n_steps = 100,
                                     n_threads = 2,
                                     lag_rates = 10,
                                     country = 'Burkina Faso',
                                     admin_unit = 'Cascades',
                                     seasonality_on = 1,
                                     state_check = 0,
                                     seasonality_check = 0,
                                     stoch_param = 'betaa')

###Troubleshoot step too small error

##Configure cluster settings
ctx_trbsht <- context::context_save("T:/jth/contexts.trbsht",
                                  packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde','RecordLinkage'),
                                  package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate','jt-hicks/sifter@issue-9')))
config_trbsht <- didehpc::didehpc_config(template = "16Core",cores =8, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
# config_dide <- didehpc::didehpc_config(template = "8Core",cores =1, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj_trbsht <- didehpc::queue_didehpc(ctx_trbsht,config = config_trbsht)
devtools::install_github('jt-hicks/sifter@25172f2')
obj_trbsht$login()
#Run a version of the sifter package that does not produce daily trajectory output
sim_rw_bulk_betaa_1108_trbsht_test <- obj_trbsht$enqueue_bulk(c(10:12), function(i,data_raw){
  print(data_raw[[i]]$data_raw)
  sifter::run_pmcmc(data_raw = data_raw[[i]]$data_raw,
                    n_particles = 10,
                    proposal_matrix = diag(0.1, 2),
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 10,
                    n_threads = 2,
                    lag_rates = 10,
                    seasonality_on = 0,
                    state_check = 0,
                    start_pf_time = 30,
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'u5',
                    seed = 1008)
},data_raw=rand_walk_data)
sim_rw_bulk_betaa_1108_trbsht$status() #'protozoological_barbet' submitted at 10:02 am on 14 Aug
sim_rw_bulk_betaa_1108_trbsht$times()
sim_rw_bulk_betaa_1108_trbsht$tasks[[2]]$log()
data_raw[[10]]$data_raw
sim_rw_bulk_betaa_1108_trbsht_test$status()
