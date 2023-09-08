remotes::install_github("mrc-ide/mcstate", upgrade = TRUE)
remotes::install_github("mrc-ide/dust", upgrade = TRUE)
remotes::install_github("mrc-ide/mode", upgrade = TRUE)
remotes::install_github("mrc-ide/odin.dust", upgrade = TRUE)
remotes::install_github("mrc-ide/odin", upgrade = TRUE)
remotes::install_github("mrc-ide/didehpc", upgrade = TRUE)
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
source('nnp/in_development/run_pmcmc_mg.R')
source('nnp/in_development/run_pmcmc_orig.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')
source('shared/utils.R')

##Import NNP data (incomplete data sets)
##Primigrav
data_raw_ng_pg_asa <- readRDS('nnp/data/data_raw_ng_pg_asa.RDS')
data_raw_ng_pg_ifenorth <- readRDS('nnp/data/data_raw_ng_pg_ifenorth.RDS')
data_raw_ng_pg_ejigbo <- readRDS('nnp/data/data_raw_ng_pg_ejigbo.RDS')
data_raw_ng_pg_moro <- readRDS('nnp/data/data_raw_ng_pg_moro.RDS')

data_raw_bf_pg_banfora <- readRDS('nnp/data/data_raw_bf_pg_banfora.RDS')
data_raw_bf_pg_orodara <- readRDS('nnp/data/data_raw_bf_pg_orodara.RDS')
data_raw_bf_pg_gaoua <- readRDS('nnp/data/data_raw_bf_pg_gaoua.RDS')

data_raw_mz_pg_guro <- readRDS('nnp/data/data_raw_mz_pg_guro.RDS')
data_raw_mz_pg_chemba <- readRDS('nnp/data/data_raw_mz_pg_chemba.RDS')%>%
  filter(as.yearmon(month) >= as.yearmon('July 2021'))
data_raw_mz_pg_changara <- readRDS('nnp/data/data_raw_mz_pg_changara.RDS')

nnp_pg_list <- list(data_raw_bf_pg_banfora,data_raw_bf_pg_gaoua,data_raw_bf_pg_orodara,
                    data_raw_mz_pg_changara,data_raw_mz_pg_chemba,data_raw_mz_pg_guro,
                    data_raw_ng_pg_asa,data_raw_ng_pg_ejigbo,data_raw_ng_pg_ifenorth,data_raw_ng_pg_moro)
names(nnp_pg_list) <- c('Banfora','Gaoua','Orodara','Changara','Chemba','Guro','Asa','Ejigbo','Ife North','Moro')

##Multigrav
data_raw_ng_mg_asa <- readRDS('nnp/data/data_raw_ng_mg_asa.RDS')
data_raw_ng_mg_ifenorth <- readRDS('nnp/data/data_raw_ng_mg_ifenorth.RDS')
data_raw_ng_mg_ejigbo <- readRDS('nnp/data/data_raw_ng_mg_ejigbo.RDS')
data_raw_ng_mg_moro <- readRDS('nnp/data/data_raw_ng_mg_moro.RDS')

data_raw_bf_mg_banfora <- readRDS('nnp/data/data_raw_bf_mg_banfora.RDS')
data_raw_bf_mg_orodara <- readRDS('nnp/data/data_raw_bf_mg_orodara.RDS')
data_raw_bf_mg_gaoua <- readRDS('nnp/data/data_raw_bf_mg_gaoua.RDS')

data_raw_mz_mg_guro <- readRDS('nnp/data/data_raw_mz_mg_guro.RDS')
data_raw_mz_mg_chemba <- readRDS('nnp/data/data_raw_mz_mg_chemba.RDS')%>%
  filter(as.yearmon(month) >= as.yearmon('July 2021'))
data_raw_mz_mg_changara <- readRDS('nnp/data/data_raw_mz_mg_changara.RDS')

nnp_mg_list <- list(data_raw_bf_mg_banfora,data_raw_bf_mg_gaoua,data_raw_bf_mg_orodara,
                    data_raw_mz_mg_changara,data_raw_mz_mg_chemba,data_raw_mz_mg_guro,
                    data_raw_ng_mg_asa,data_raw_ng_mg_ejigbo,data_raw_ng_mg_ifenorth,data_raw_ng_mg_moro)
names(nnp_mg_list) <- c('Banfora','Gaoua','Orodara','Changara','Chemba','Guro','Asa','Ejigbo','Ife North','Moro')

country <- c('Burkina Faso','Burkina Faso','Burkina Faso',
             'Mozambique','Mozambique','Mozambique',
             'Nigeria','Nigeria','Nigeria','Nigeria')
admin <- c('Cascades','Sud-Ouest','Haut-Bassins',
           'Tete','Sofala','Manica',
           'Kwara','Osun','Osun','Kwara')
##Test run_pmcmc function##
source('nnp/in_development/run_pmcmc_mg.R')

test_run_mg_corr <- run_pmcmc_mg(data_raw_pg = data_raw_bf_pg_banfora,
                                 data_raw_mg = data_raw_bf_mg_banfora,
                                 n_particles = 10,
                                 proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                 max_EIR=1000,
                                 max_steps = 1e7,
                                 atol = 1e-5,
                                 rtol = 1e-6,
                                 n_steps = 3,
                                 n_threads = 2,
                                 lag_rates = 10,
                                 country = 'Burkina Faso',
                                 admin_unit = 'Cascades',
                                 seasonality_on = 0,
                                 state_check = 0)
windows(10,7)
plot_particle_filter(test_run_mg_corr$history,true_history=data_raw_bf_pg_banfora,times=data_raw_bf_pg_banfora$t)
dim(t(test_run_mg_corr$history[1,,-1]))

##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("shared/run_pmcmc.R",
             "nnp/in_development/run_pmcmc_pg.R",
             "nnp/in_development/run_pmcmc_mg.R",
             "shared/model_parameters.R",
             "shared/equilibrium-init-create-stripped.R",
             'shared/utils.R')
ctx <- context::context_save("T:/jth/contexts.temp", sources = sources,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust@8aef08d",'mrc-ide/mcstate')))
config_1 <- didehpc::didehpc_config(template = "32Core",cores =4, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
config_single <- didehpc::didehpc_config(template = "GeneralNodes",cores =1, parallel = FALSE,wholenode = FALSE, cluster = 'fi--didemrchnb')
config_dide <- didehpc::didehpc_config(template = "8Core",cores =4, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj <- didehpc::queue_didehpc(ctx,config = config_1)
obj <- didehpc::queue_didehpc(ctx,config = config_single)
obj$login()
obj$cluster_load(TRUE)

##Submit bulk non-seasonality runs
nnp_mgcorr_bulk_std_test <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg){
  run_pmcmc_mg(data_raw_pg = data_pg[[i]],
               data_raw_mg = data_mg[[i]],
               n_particles = 10,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 10,
               n_threads = 8,
               lag_rates = 10,
               seasonality_on = 0,
               state_check = 0)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list)
nnp_mgcorr_bulk_std_test$status() #'banal_eagle'
nnp_mgcorr_bulk_std_test$tasks[[1]]$log()
obj$unsubmit(nnp_pgcorr_bulk_std_test$ids)

nnp_mgcorr_bulk_std <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg){
  run_pmcmc_mg(data_raw_pg = data_pg[[i]],
               data_raw_mg = data_mg[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 8,
               lag_rates = 10,
               seasonality_on = 0,
               state_check = 0)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list)
nnp_mgcorr_bulk_std$status() #'trifling_bongo'
nnp_mgcorr_bulk_std <- obj$task_bundle_get('trifling_bongo')
nnp_mgcorr_bulk_std_results <- lapply(1:10, function(id){
  nnp_mgcorr_bulk_std$tasks[[id]]$result()
})

nnp_ng_mgcorr_bulk_std <- obj$enqueue_bulk(7:10, function(i,data_pg,data_mg){
  run_pmcmc_mg(data_raw_pg = data_pg[[i]],
               data_raw_mg = data_mg[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 1,
               lag_rates = 10,
               seasonality_on = 0,
               state_check = 0)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list)
nnp_ng_mgcorr_bulk_std$status() #'flannel_easternglasslizard' submitted 2 Mar 10:48am

nnp_mgcorr_bulk_seas <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg,country,admin){
  run_pmcmc_mg(data_raw_pg = data_pg[[i]],
               data_raw_mg = data_mg[[i]],
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
               country = country[i],
               admin_unit = admin[i],
               state_check = 0)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_mgcorr_bulk_seas$status() #ittybitty_atlanticblackgoby - submitted 28 Feb 9:27am
#'amateurish_elver' - submitted 27 Feb 10:25am - [1] completed
nnp_mgcorr_bulk_seas$tasks[[1]]$log()
nnp_mgcorr_bulk_seas_270223 <- obj$task_bundle_get('amateurish_elver')
nnp_mgcorr_bulk_seas_270223$times()
nnp_mgcorr_bulk_bf_seas$status() #'chemophobic_macaque'
obj$unsubmit(nnp_mgcorr_bulk_seas$ids)
obj$login()
obj$cluster_load(TRUE)
nnp_mgcorr_bulk_bf_seas$tasks[[1]]$log()
obj$config

nnp_mgcorr_bulk_seas_single <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg,country,admin){
  run_pmcmc_mg(data_raw_pg = data_pg[[i]],
               data_raw_mg = data_mg[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 1,
               lag_rates = 10,
               seasonality_on = 1,
               country = country[i],
               admin_unit = admin[i],
               state_check = 0)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_mgcorr_bulk_seas_single$status() #'monotonous_betafish' - submitted 1 Mar 9:35am
nnp_mgcorr_bulk_seas_single$times()
nnp_mgcorr_bulk_seas_single <- obj$task_bundle_get('monotonous_betafish')
nnp_mgcorr_bulk_seas_results <- lapply(1:6, function(id){
  nnp_mgcorr_bulk_seas_single$tasks[[id]]$result()
})
obj$unsubmit(nnp_mgcorr_bulk_seas_single$ids[7:10])

nnp_ng_mgcorr_bulk_seas_single <- obj$enqueue_bulk(7:10, function(i,data_pg,data_mg,country,admin){
  run_pmcmc_mg(data_raw_pg = data_pg[[i]],
               data_raw_mg = data_mg[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 1,
               lag_rates = 10,
               seasonality_on = 1,
               country = country[i],
               admin_unit = admin[i],
               state_check = 0)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_ng_mgcorr_bulk_seas_single$status() #'antimonarchal_flee' submitted 2 Mar 10:40am
nnp_ng_mgcorr_bulk_seas_single <- obj$task_bundle_get('antimonarchal_flee')
nnp_ng_mgcorr_bulk_seas_results <- lapply(1:4, function(id){
  nnp_ng_mgcorr_bulk_seas_single$tasks[[id]]$result()
})
obj$task_bundle_list()
nnp_mgcorr_bulk_seas_results_update <- append(nnp_mgcorr_bulk_seas_results,nnp_ng_mgcorr_bulk_seas_results)


###### Run with sifter ######
##Configure cluster settings
ctx_sifter2 <- context::context_save("T:/jth/contexts.sift2",
                                    packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde','RecordLinkage'),
                                    package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate','jt-hicks/sifter')))
config_1 <- didehpc::didehpc_config(template = "24Core",cores =6, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
# config_dide <- didehpc::didehpc_config(template = "8Core",cores =1, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj <- didehpc::queue_didehpc(ctx_sifter2,config = config_1)
obj$cluster_load(TRUE)
obj$login()

nnp_pgmg_bulk_sifter_eir <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg,country,admin){
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
               comparison = 'pgmg')
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_pgmg_bulk_sifter_eir$status() #'creepy_oropendola' - submitted 17 July 10:40am
nnp_pgmg_bulk_sifter_eir$tasks[[7]]$log()
nnp_mgcorr_bulk_snnp_pgmg_bulk_sifter_eireas_single$times()
nnp_pgmg_bulk_sifter_eir <- obj$task_bundle_get('creepy_oropendola')
nnp_pgmg_bulk_sifter_eir$tasks[[7]] <- nnp_pgmg_bulk_sifter_eir_7.9$tasks[[1]]
nnp_pgmg_bulk_sifter_eir$tasks[[9]] <- nnp_pgmg_bulk_sifter_eir_7.9$tasks[[2]]

nnp_pgmg_bulk_sifter_eir_results <- lapply(c(1:10), function(id){
  nnp_pgmg_bulk_sifter_eir$tasks[[id]]$result()
})
names(nnp_pgmg_bulk_sifter_eir_results) <- names(nnp_pg_list)
results <- nnp_pgmg_bulk_sifter_eir_results
nnp_pgmg_bulk_sifter_eir_results[[1]]$history[,1,1]

nnp_pgmg_bulk_sifter_eir_7.9 <- obj$enqueue_bulk(c(7,9), function(i,data_pg,data_mg,country,admin){
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
                    comparison = 'pgmg',
                    seed = 1307)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_pgmg_bulk_sifter_eir_7.9$status() #floatable_gerbil (18 July 10:22am)

nnp_pgmg_bulk_sifter_eir_max2000 <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg,country,admin){
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
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_pgmg_bulk_sifter_eir_max2000$status() #'interuniversity_teledu' - submitted 21 July 10:19am
nnp_pgmg_bulk_sifter_eir_max2000_results <- lapply(c(1:10), function(id){
  nnp_pgmg_bulk_sifter_eir_max2000$tasks[[id]]$result()
})
names(nnp_pgmg_bulk_sifter_eir_max2000_results) <- c('Banfora','Gaoua','Orodara','Changara','Chemba','Guro','Asa','Ejigbo','Ife North','Moro')


### With mosquito emergence
##Some tests with the max betaa
nnp_pgmg_bulk_sifter_betaa_test1 <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg,country,admin){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 100,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_param=1000,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 50,
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
                    seed = 1707)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_pgmg_bulk_sifter_betaa_test1$status() #'anthropoid_paca' - submitted 20 July 10:26am
nnp_pgmg_bulk_sifter_betaa_test1$tasks[[1]]$log()
nnp_pgmg_bulk_sifter_betaa_test1_results <- lapply(c(1:10), function(id){
  nnp_pgmg_bulk_sifter_betaa_test1$tasks[[id]]$result()
})

nnp_pgmg_bulk_sifter_betaa_test2 <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg,country,admin){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 100,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_param=500,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 50, 
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
                    seed = 1707)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_pgmg_bulk_sifter_betaa_test2$status() #'melittological_tuna' - submitted 20 July 10:29am
nnp_pgmg_bulk_sifter_betaa_test2$tasks[[10]]$result()$history['betaa',50,]
nnp_pgmg_bulk_sifter_betaa_test2_results <- lapply(c(1:10), function(id){
  nnp_pgmg_bulk_sifter_betaa_test2$tasks[[id]]$result()
})

nnp_pgmg_bulk_sifter_betaa <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg,country,admin){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_param=62, #~1000/16
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
                    seed = 1707)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_pgmg_bulk_sifter_betaa$status() #'ununtrium_zebrafinch' - submitted 17 July 10:41am
nnp_pgmg_bulk_sifter_betaa$status() #'vwritten_johndory' - submitted 20 July 2:24pm

nnp_pgmg_bulk_sifter_betaa$times()
nnp_pgmg_bulk_sifter_betaa$tasks[['ab8339b81e36f767dd37e6e38aab436e']]$log()
nnp_pgmg_bulk_sifter_betaa$tasks[[7]]$log()
nnp_pgmg_bulk_sifter_betaa$tasks$e988e0d986be894d5cfb3b26f288ee35$log()
nnp_pgmg_bulk_sifter_betaa <- obj$task_bundle_get('ununtrium_zebrafinch')
nnp_pgmg_bulk_sifter_betaa$tasks[[7]] <- nnp_pgmg_bulk_sifter_betaa_7$tasks[[1]]

nnp_pgmg_bulk_sifter_betaa_results <- lapply(c(1:10), function(id){
  nnp_pgmg_bulk_sifter_betaa$tasks[[id]]$result()
})
names(nnp_pgmg_bulk_sifter_betaa_results) <- c('Banfora','Gaoua','Orodara','Changara','Chemba','Guro','Asa','Ejigbo','Ife North','Moro')

View(nnp_pgmg_bulk_sifter_betaa_results)
##Resubmit 7 and 10
nnp_pgmg_bulk_sifter_betaa_7 <- obj$enqueue_bulk(c(7), function(i,data_pg,data_mg,country,admin){
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
                    stoch_param = 'betaa',
                    comparison = 'pgmg',
                    seed=1307)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_pgmg_bulk_sifter_betaa_7$status() #'discoloured_easteuropeanshepherd'
nnp_pgmg_bulk_sifter_betaa_7 <- obj$task_bundle_get('discoloured_easteuropeanshepherd')
nnp_pgmg_bulk_sifter_betaa_7$tasks[[1]]$log()
asa_betaa <- nnp_pgmg_bulk_sifter_betaa_7$tasks[[1]]$result()
asa_betaa$history['betaa',500,-1]
asa_betaa$history[,500,1]
length(asa_betaa$history['betaa',525,-1])
plot(asa_betaa$history['betaa',500:1000,60:883],asa_betaa$history['EIR',500:1000,60:883])
plot(asa_betaa$pars$EIR_SD)
plot(rnorm(1000,1)*asa_betaa$pars$EIR_SD)

random_dev <- rnorm(1000,1)*asa_betaa$pars$EIR_SD
proposed_betaa <- sapply(60:883,function(x) exp(log(asa_betaa$history['betaa',500:1000,x]) + random_dev[500:1000]))
betaa_max <- log(30)
proposed_betaa_wmax <- sapply(60:883,function(x) exp(min(log(asa_betaa$history['betaa',500:1000,x]) + random_dev[500:1000],betaa_max)))

plot(exp(log(asa_betaa$history['betaa',500:1000,60:883])+random_dev(500:1000)))
plot(proposed_betaa,proposed_betaa_wmax)
plot(proposed_betaa_wmax)

obj$login()
nnp_pgmg_bulk_sifter_betaa_max125 <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg,country,admin){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_param=125, #~2000/16
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
                    seed = 1707)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_pgmg_bulk_sifter_betaa_max125 <- obj$task_bundle_get('webbed_sable')
nnp_pgmg_bulk_sifter_betaa_max125$status() #'webbed_sable' submitted 21Jul 10:16am
nnp_pgmg_bulk_sifter_betaa_max125_results <- lapply(c(1:10), function(id){
  nnp_pgmg_bulk_sifter_betaa_max125$tasks[[id]]$result()
})
names(nnp_pgmg_bulk_sifter_betaa_max125_results) <- c('Banfora','Gaoua','Orodara','Changara','Chemba','Guro','Asa','Ejigbo','Ife North','Moro')
proposal_matrix <- cov(nnp_pgmg_bulk_sifter_betaa_max125_results[[1]]$pars)
proposal_matrix_adj <- cov(nnp_pgmg_bulk_sifter_betaa_max125_results[[1]]$pars)
cov(nnp_pgmg_bulk_sifter_betaa_max125_results[[2]]$pars)
prop_mat_list <- lapply(c(1:10),function(i) cov(nnp_pgmg_bulk_sifter_betaa_max125_results[[i]]$pars))

obj$login()
nnp_pgmg_bulk_sifter_betaa_2407 <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg,country,admin,prop_mat_list){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = prop_mat_list[[i]],
                    max_param=125, #~2000/16
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
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin,prop_mat_list=prop_mat_list)
nnp_pgmg_bulk_sifter_betaa_2407$status() #'griffinish_killerwhale' submitted 24 July at 3:37pm
nnp_pgmg_bulk_sifter_betaa_2407$tasks[[1]]$log()
nnp_pgmg_bulk_sifter_betaa_2407 <- obj$task_bundle_get('griffinish_killerwhale')
nnp_pgmg_bulk_sifter_betaa_2407$tasks[[2]] <- nnp_pgmg_bulk_sifter_betaa_2407.2$tasks[[1]]
nnp_pgmg_bulk_sifter_betaa_2407$tasks[[4]] <- nnp_pgmg_bulk_sifter_betaa_2407.mz$tasks[[1]]
nnp_pgmg_bulk_sifter_betaa_2407$tasks[[5]] <- nnp_pgmg_bulk_sifter_betaa_2407.mz$tasks[[2]]
nnp_pgmg_bulk_sifter_betaa_2407$tasks[[6]] <- nnp_pgmg_bulk_sifter_betaa_2407.mz$tasks[[3]]
nnp_pgmg_bulk_sifter_betaa_2407$status()
nnp_pgmg_bulk_sifter_betaa_2407_results <- lapply(c(1:10), function(i){
  nnp_pgmg_bulk_sifter_betaa_2407$tasks[[i]]$result()
})
names(nnp_pgmg_bulk_sifter_betaa_2407_results) <- names(nnp_pg_list)
prop_mat_list_2407 <- lapply(c(1:10),function(i) cov(nnp_pgmg_bulk_sifter_betaa_2407_results[[i]]$pars))

nnp_pgmg_bulk_sifter_betaa_2407.2 <- obj$enqueue_bulk(2, function(i,data_pg,data_mg,country,admin,prop_mat_list){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = prop_mat_list[[i]],
                    max_param=125, #~2000/16
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
                    seed = 2507)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin,prop_mat_list=prop_mat_list)
nnp_pgmg_bulk_sifter_betaa_2407.2$status() #unsanitary_wireworm
nnp_pgmg_bulk_sifter_betaa_2407.2 <- obj$task_bundle_get('unsanitary_wireworm')

obj$login()
nnp_pgmg_bulk_sifter_betaa_2407.mz <- obj$enqueue_bulk(4:6, function(i,data_pg,data_mg,country,admin,prop_mat_list){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = prop_mat_list[[i]],
                    max_param=125, #~2000/16
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
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin,prop_mat_list=prop_mat_list)
nnp_pgmg_bulk_sifter_betaa_2407.mz$status()#ununbium_rottweiler
nnp_pgmg_bulk_sifter_betaa_2407.mz <- obj$task_bundle_get('ununbium_rottweiler')

nnp_pgmg_bulk_sifter_betaa_u5 <- obj$enqueue_bulk(1:10, function(i,data_pg,country,admin,prop_mat_list){
  sifter::run_pmcmc(data_raw = data_pg[[i]],
                    n_particles = 200,
                    proposal_matrix = prop_mat_list[[i]],
                    max_param=125, #~2000/16
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
                    comparison = 'u5',
                    seed = 2407)
},data_pg=nnp_pg_list,country=country,admin=admin,prop_mat_list=prop_mat_list_2407)
nnp_pgmg_bulk_sifter_betaa_u5$status() #unnecessary_gull
nnp_pgmg_bulk_sifter_betaa_u5 <- obj$task_bundle_get('unnecessary_gull')
nnp_pgmg_bulk_sifter_betaa_u5$times()
nnp_pgmg_bulk_sifter_betaa_u5$tasks[[1]]$log()
##Test workers
ctx_sifter.work <- context::context_save("T:/jth/contexts.workers",
                                         packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde','RecordLinkage','parallel'),
                                         package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate','jt-hicks/sifter')))
config_8 <- didehpc::didehpc_config(template = "16Core",cores =8, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
# config_dide <- didehpc::didehpc_config(template = "8Core",cores =1, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj_workers <- didehpc::queue_didehpc(ctx_sifter.work,config = config_8)
obj_workers$login()
nnp_pgmg_bulk_sifter_betaa_workertest <- obj_workers$enqueue_bulk(1, function(i,data_pg,data_mg,country,admin,prop_mat_list){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = prop_mat_list[[i]],
                    max_param=125, #~2000/16
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 8,
                    n_chains = 2,
                    n_workers = 1,
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
                    seed = 2707)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin,prop_mat_list=prop_mat_list)
nnp_pgmg_bulk_sifter_betaa_workertest$status() # overcooked_marmot submitted 28 Jul 11:08am
nnp_pgmg_bulk_sifter_betaa_workertest$status() #tardy_flyingsquirrel submitted 9 Aug 1:42pm
nnp_pgmg_bulk_sifter_betaa_workertest$tasks[[1]]$log()
obj_workers$unsubmit(nnp_pgmg_bulk_sifter_betaa_workertest$ids)
obj_workers$cluster_load(TRUE)
open_mp_check <- obj_workers$enqueue(dust::dust_openmp_support())
open_mp_check$status()
open_mp_check$result()
Sys.getenv("MC_CORES")
library(parallel)

##Run long runs
##Configure cluster settings
ctx_long <- context::context_save("T:/jth/contexts.long",
                                     packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde','RecordLinkage'),
                                     package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate','jt-hicks/sifter')))
config_long <- didehpc::didehpc_config(template = "16Core",cores =8, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
# config_dide <- didehpc::didehpc_config(template = "8Core",cores =1, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj_long <- didehpc::queue_didehpc(ctx_long,config = config_long)
obj_long$cluster_load(TRUE)
obj_long$login()

nnp_pgmg_bulk_sifter_betaa_1008 <- obj_long$enqueue_bulk(1:10, function(i,data_pg,data_mg,country,admin,prop_mat_list){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = prop_mat_list[[i]],
                    max_param=125, #~2000/16
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
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin,prop_mat_list=prop_mat_list_2407)
nnp_pgmg_bulk_sifter_betaa_1008$status() #'cardiac_archerfish' submitted 10 Aug at 2:02pm
# obj$unsubmit(nnp_pgmg_bulk_sifter_betaa_2407$ids)
nnp_pgmg_bulk_sifter_betaa_1008$tasks[[1]]$log()
nnp_pgmg_bulk_sifter_betaa_1008$tasks[[5]]$log()

obj_long$login()
nnp_pgmg_bulk_sifter_betaa_1008.2_5 <- obj_long$enqueue_bulk(c(2,5), function(i,data_pg,data_mg,country,admin,prop_mat_list){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = prop_mat_list[[i]],
                    max_param=125, #~2000/16
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
                    seed = 1108)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin,prop_mat_list=prop_mat_list_2407)
nnp_pgmg_bulk_sifter_betaa_1008.2_5$status() #'nepotistic_chick' submitted 11 Aug 2023

nnp_pgmg_bulk_sifter_betaa_1008.8 <- obj_long$enqueue_bulk(8, function(i,data_pg,data_mg,country,admin,prop_mat_list){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = prop_mat_list[[i]],
                    max_param=125, #~2000/16
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
                    seed = 1108)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin,prop_mat_list=prop_mat_list_2407)
nnp_pgmg_bulk_sifter_betaa_1008.8$status() #'glazed_boaconstrictor' submitted 12 Aug 2023

##Test
source('nnp/in_development/run_pmcmc_orig.R')

test_orig <- run_pmcmc_orig(data_raw = sifter::data_sim, #I've added data_sim to the package for an easy test
                       n_particles = 10,
                       proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                       max_EIR=1000,
                       max_steps = 1e7,
                       atol = 1e-5,
                       rtol = 1e-6,
                       n_steps = 10,
                       n_threads = 2,
                       lag_rates = 10,
                       country = 'Burkina Faso',
                       admin_unit = 'Cascades',
                       seasonality_on = 0,
                       state_check = 0,
                       seasonality_check = 0)

##Configure cluster settings
ctx_initial <- context::context_save("T:/jth/contexts.init",
                                  packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde','RecordLinkage'),
                                  package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate','jt-hicks/sifter@issue-10')))
config_init <- didehpc::didehpc_config(template = "16Core",cores=8, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
# config_dide <- didehpc::didehpc_config(template = "8Core",cores =1, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj_init <- didehpc::queue_didehpc(ctx_initial,config = config_init)
obj_long$cluster_load(TRUE)
obj_long$login()

nnp_pgmg_bulk_sifter_betaa_3008 <- obj_init$enqueue_bulk(1:10, function(i,data_pg,data_mg,country,admin,prop_mat_list){
  sifter::run_pmcmc(data_raw_pg = data_pg[[i]],
                    data_raw_mg = data_mg[[i]],
                    n_particles = 200,
                    proposal_matrix = prop_mat_list,
                    max_param=125, #~2000/16
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
                    seasonality_check = 0,
                    stoch_param = 'betaa',
                    comparison = 'pgmg',
                    seed = 1008)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin,prop_mat_list=proposal_mat)
nnp_pgmg_bulk_sifter_betaa_3008$status() #'seismoscopic_pitbull' submitted 31 Aug at 10:01am
obj_init$unsubmit(nnp_pgmg_bulk_sifter_betaa_3008$ids)
nnp_pgmg_bulk_sifter_betaa_3008$tasks[[1]]$log()
obj_init$login()
