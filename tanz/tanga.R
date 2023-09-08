library(foreign)
library(sas7bdat)
library(tidyverse)
library(readxl)
library(devtools)
devtools::install_github("mrc-ide/umbrella")
library(umbrella)
library(ggplot2)
library(lubridate)
library(zoo)
library(terra)
library(stringr)
library(geofacet)
library(reshape2)

theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))
#See tanz_read_data.R for how to pull and organize Tanga data
test_tanga <- sifter::run_pmcmc(data_raw = tanga_data_list_15to22[[1]], 
                                n_particles = 10,
                                proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                max_EIR=1000,
                                max_steps = 1e7,
                                atol = 1e-5,
                                rtol = 1e-6,
                                n_steps = 3,
                                n_threads = 2,
                                lag_rates = 10,
                                country = 'Tanzania',
                                admin_unit = 'Tanga',
                                seasonality_on = 1,
                                state_check = 0,
                                seasonality_check = 0,
                                preyears = 5,
                                start_pf_time = 30*4
                                )

##Configure cluster settings
ctx_sifter <- context::context_save("T:/jth/contexts.sift",
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde','RecordLinkage'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate','jt-hicks/sifter@issue-5')))
config_1 <- didehpc::didehpc_config(template = "24Core",cores =6, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
# config_dide <- didehpc::didehpc_config(template = "8Core",cores =1, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj_sift <- didehpc::queue_didehpc(ctx_sifter,config = config_1)
obj_sift$cluster_load(TRUE)
obj_sift$login()
##Submit test to cluster
test_tanga_clust <- obj_sift$enqueue(sifter::run_pmcmc(data_raw = tanga_data_list_15to22[[1]], 
                                n_particles = 10,
                                proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                max_EIR=1000,
                                max_steps = 1e7,
                                atol = 1e-5,
                                rtol = 1e-6,
                                n_steps = 3,
                                n_threads = 2,
                                lag_rates = 10,
                                country = 'Tanzania',
                                admin_unit = 'Tanga',
                                seasonality_on = 1,
                                state_check = 0,
                                seasonality_check = 0,
                                preyears = 5,
                                start_pf_time = 30*4
))
test_tanga_clust$status()
test_tanga_clust$log()

test_tanga_clust <- obj_sift$enqueue_bulk(1:11, function(i,data_list) {sifter::run_pmcmc(data_raw = data_list[[i]], 
                                                       n_particles = 10,
                                                       proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                                       max_EIR=1000,
                                                       max_steps = 1e7,
                                                       atol = 1e-5,
                                                       rtol = 1e-6,
                                                       n_steps = 3,
                                                       n_threads = 2,
                                                       lag_rates = 10,
                                                       country = 'Tanzania',
                                                       admin_unit = 'Tanga',
                                                       seasonality_on = 1,
                                                       state_check = 0,
                                                       seasonality_check = 0,
                                                       preyears = 5,
                                                       start_pf_time = 30*4
)},data_list = tanga_data_list_15to22)
test_tanga_clust$status()#'dramatisable_muskrat'

tanga_clust_290623 <- obj_sift$enqueue_bulk(1:11, function(i,data_list) {
  sifter::run_pmcmc(data_raw = data_list[[i]], 
                    n_particles = 200,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_EIR=1000,
                    max_steps = 1e7,
                    atol = 1e-5,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    country = 'Tanzania',
                    admin_unit = 'Tanga',
                    seasonality_on = 1,
                    state_check = 0,
                    seasonality_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4
  )},data_list = tanga_data_list_15to22)
tanga_clust_290623$status()#'chewable_bellfrog'
tanga_clust_290623 <- obj_sift$task_bundle_get('chewable_bellfrog')
tanga_clust_290623$tasks[[1]]$log()
tanga_clust_290623_results <- lapply(2:11, function(id){
  tanga_clust_290623$tasks[[id]]$result()
})
names(tanga_data_list_15to22)
names(tanga_clust_290623_results) <- names(tanga_data_list_15to22[2:11])

##re-run first district
tanga_clust_290623_redo1 <- obj_sift$enqueue_bulk(1, function(i,data_list) {
  sifter::run_pmcmc(data_raw = data_list[[i]], 
                    n_particles = 200,
                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                    max_EIR=1000,
                    max_steps = 1e7,
                    atol = 1e-5,
                    rtol = 1e-6,
                    n_steps = 1000,
                    n_threads = 4,
                    lag_rates = 10,
                    country = 'Tanzania',
                    admin_unit = 'Tanga',
                    seasonality_on = 1,
                    state_check = 0,
                    seasonality_check = 0,
                    preyears = 5,
                    start_pf_time = 30*4,
                    seed = 123
  )},data_list = tanga_data_list_15to22)
tanga_clust_290623_redo1$status()#'nosophobic_argali'
tanga_clust_290623_redo1 <- obj_sift$task_bundle_get('nosophobic_argali')
tanga_clust_290623$tasks[[1]] <- tanga_clust_290623_redo1$tasks[[1]]

tanga_clust_290623_results <- lapply(1:11, function(id){
  tanga_clust_290623$tasks[[id]]$result()
})

names(tanga_clust_290623_results) <- names(tanga_data_list_15to22[1:11])


###Summarize
rm(diagnostic_plots)
diagnostic_plots <- create_diag_figs_bulk(results=tanga_clust_290623_results)
windows(10,5)
diagnostic_plots$ar.plot
diagnostic_plots$ess.plot
diagnostic_plots$posterior.trace
diagnostic_plots$EIR_SD.trace
diagnostic_plots$init_EIR.trace
diagnostic_plots$EIR_SD.density
diagnostic_plots$init_EIR.density
tanga_rainfall <- read_csv('./tanz/Patrick/processed_inputs/rainfall_ad2.csv')%>%
  filter(NAME_1 == 'Tanga')%>%
  select(Month,NAME_1,NAME_2,ID,Rainfall)
tanga_itn <- read_csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Tanzania/Tanga/tanga_itn_access.csv')
names(tanga_itn) <- c('region','sites','date','itn_access')
source('./tanz/create_summary_plots.R')
tanga_plots <- create_summary_plots(results = tanga_clust_290623_results,
                                    data_list = tanga_data_list_15to22[1:11],
                                    level = 'Council',
                                    rainfall = tanga_rainfall,
                                    itn_access = tanga_itn,
                                    start_pf_time = 30*4,
                                    date_limits = c('2018-01-01',NA))
windows(10,10)
obs <- tanga_plots$obs_prev_plot
obs
ggsave('tanz/figures/obsprev_tanga_290623.pdf',plot = obs,width = 9,height=6)

est_prev <- tanga_plots$est_prev_plot
est_prev
ggsave('tanz/figures/estprev_tanga_290623.pdf',plot = est_prev,width = 9,height=6)

inc <- tanga_plots$inc.rainfall
inc
ggsave('tanz/figures/inc_tanga_290623.pdf',plot = inc,width = 9,height=6)

rel_inc <- tanga_plots$inc.rainfall.3
rel_inc
ggsave('tanz/figures/relinc_tanga_290623.pdf',plot = rel_inc,width = 9,height=6)

eir <- tanga_plots$eir.rainfall
eir
ggsave('tanz/figures/eir_tanga_290623.pdf',plot = eir,width = 9,height=6)

eir_sd <- tanga_plots$EIR_SD.density
eir_sd
ggsave('tanz/figures/eirsd_tanga_290623.pdf',plot = eir_sd,width = 9,height=6)

init_eir <- tanga_plots$init_EIR.density
init_eir
ggsave('tanz/figures/initEIR_tanga_290623.pdf',plot = init_eir,width = 9,height=6)
