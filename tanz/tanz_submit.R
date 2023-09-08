library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(didehpc)
library(pkgdepends)
library(dplyr)
library("coda")
library(binom)
library(ggplot2)
library(bayesplot)
library(zoo)
library(lubridate)
library(reshape2)

##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("shared/run_pmcmc.R",
             "shared/model_parameters.R",
             "shared/equilibrium-init-create-stripped.R",
             'shared/utils.R')
ctx <- context::context_save("T:/jth/contexts.temp", sources = sources,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate')))
ctx <- context::context_save("T:/jth/contexts.190423", sources = sources,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate')))
ctx <- context::context_save("T:/jth/contexts.230623", sources = sources,
                                    packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'),
                                    package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate')))
ctx <- context::context_save("T:/jth/contexts_may23", sources = sources,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate')))

config_1 <- didehpc::didehpc_config(template = "24Core",cores =6, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
# config_dide <- didehpc::didehpc_config(template = "8Core",cores =1, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj <- didehpc::queue_didehpc(ctx,config = config_1)
obj$login()
obj$cluster_load(TRUE)
names()

test_tanz_submit <- obj$enqueue_bulk(1:2, function(x,obs_list) {
  print(names(obs_list$lt20[x]))
  run_pmcmc(data = obs_list$lt20[[x]],
            n_particles = 10,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 100,
            n_threads = 4,
            lag_rates = 10,
            preyears = 5,
            seasonality_on = 0,
            state_check = 0,
            seasonality_check = 0)},obs_list=tanz_data_list_14to17)

test_tanz_submit$status()
test_tanz_submit$tasks$c788477447acd71ee8be7aba3a6ba47c$log()

tanz_lt20_2015to2017_submit <- obj$enqueue_bulk(1:26, function(x,obs_list) {
  print(names(obs_list$lt20[x]))
  run_pmcmc(data = obs_list$lt20[[x]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 4,
            lag_rates = 10,
            preyears = 5,
            seasonality_on = 0,
            state_check = 0,
            seasonality_check = 0)},obs_list=tanz_data_list_15to17)
tanz_lt20_2015to2017_submit$status()#'nonspiritous_lobo'
tanz_lt20_2015to2017_submit <- obj$task_bundle_get('nonspiritous_lobo')
test_result <- tanz_lt20_2015to2017_submit$tasks$`0519fdb81fdca79d98dcc1fa2b19cd3e`$result()

tanz_lt20_2015to2017_results <- lapply(1:26, function(id){
  tanz_lt20_2015to2017_submit$tasks[[id]]$result()
})

names(tanz_lt20_2015to2017_results) <- gsub(' Region','',names(tanz_data_list_15to17$lt20))


##Submit District level Lake Malawi
obj$login()
obj$cluster_load(TRUE)
obj$config

tanz_lakemal_lt20_2015to2017_submit <- obj$enqueue_bulk(1:14, function(x,obs_list) {
  print(names(obs_list$lt20[x]))
  run_pmcmc(data = obs_list$lt20[[x]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 4,
            lag_rates = 10,
            preyears = 5,
            seasonality_on = 0,
            state_check = 0,
            seasonality_check = 0)},obs_list=tanz_lakemal_list_15to17)

tanz_lakemal_lt20_2015to2017_submit$status() #'neutral_asianelephant'
tanz_lakemal_lt20_2015to2017_submit <- obj$task_bundle_get('neutral_asianelephant')
tanz_lakemal_lt20_2015to2017_submit$tasks$`754d1d793898e799ba97f07b42a06f73`$log()
tanz_lakemal_lt20_2015to2017_submit$times()
as.vector(tanz_lakemal_lt20_2015to2017_submit$status())

##Re-run tasks 8 and 13
obj$login()
obj$cluster_load(TRUE)
tanz_lakemal_lt20_2015to2017_submit_2 <- obj$enqueue_bulk(c(8,11), function(x,obs_list) {
  print(names(obs_list$lt20[x]))
  run_pmcmc(data = obs_list$lt20[[x]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 4,
            lag_rates = 10,
            preyears = 5,
            seasonality_on = 0,
            state_check = 0,
            seasonality_check = 0)},obs_list=tanz_lakemal_list_15to17)
obj$task_bundle_list()
bundle_info <- obj$task_bundle_info()
tanz_lakemal_lt20_2015to2017_submit_2$status()
tanz_lakemal_lt20_2015to2017_submit_2 <- obj$task_bundle_get('spiritless_germanshepherd')
tanz_lakemal_lt20_2015to2017_submit_2$status()
tanz_lakemal_lt20_2015to2017_submit_2$tasks[[2]]$log()
tanz_lakemal_lt20_2015to2017_submit$tasks[[8]] <- tanz_lakemal_lt20_2015to2017_submit_2$tasks[[1]]
tanz_lakemal_lt20_2015to2017_submit$tasks[[8]]$log()
tanz_lakemal_lt20_2015to2017_submit$tasks[[11]] <- tanz_lakemal_lt20_2015to2017_submit_2$tasks[[2]]
tanz_lakemal_lt20_2015to2017_submit$status()
tanz_lakemal_lt20_2015to2017_results <- lapply(1:14, function(id){
  tanz_lakemal_lt20_2015to2017_submit$tasks[[id]]$result()
})
names(tanz_lakemal_lt20_2015to2017_results) <- names(tanz_lakemal_list_15to17$lt20)

tanz_lakemal_lt20_2015to2017_results$`Ludewa District Council`$mcmc$log_prior

##Songea District Investigation
songeadist_test <- run_pmcmc(data = tanz_lakemal_list_15to17$lt20[[11]],
          n_particles = 200,
          proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
          max_EIR=1000,
          max_steps = 1e7,
          atol = 1e-5,
          rtol = 1e-6,
          n_steps = 100,
          n_threads = 4,
          lag_rates = 10,
          preyears = 5,
          seasonality_on = 0,
          state_check = 0,
          seasonality_check = 0)
tanz_lakemal_lt20_2015to2017_songeadsit_24 <- obj$enqueue_bulk(c(11), function(x,obs_list) {
  print(names(obs_list$lt20[x]))
  run_pmcmc(data = obs_list$lt20[[x]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 4,
            lag_rates = 10,
            preyears = 5,
            seasonality_on = 0,
            state_check = 0,
            seasonality_check = 0)},obs_list=tanz_lakemal_list_15to17)
#'waterbreathing_tadpole'
tanz_lakemal_lt20_2015to2017_songeadsit$status()
tanz_lakemal_lt20_2015to2017_songeadsit_24$status() #"cultivable_betafish"

tanz_lakemal_lt20_2015to2017_songeadsit_single <- obj$enqueue_bulk(c(11), function(x,obs_list) {
  print(names(obs_list$lt20[x]))
  run_pmcmc(data = obs_list$lt20[[x]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 4,
            lag_rates = 10,
            preyears = 5,
            seasonality_on = 0,
            state_check = 0,
            seasonality_check = 0,
            seed=2L)},obs_list=tanz_lakemal_list_15to17)
tanz_lakemal_lt20_2015to2017_songeadsit_single$status()#'megamalophilic_harborseal'
#'restful_okapi'
tanz_lakemal_lt20_2015to2017_songeadsit_single$tasks[[1]]$log()
obj$login()

##Submit new data
tanz_all_2015to2022_submit <- obj$enqueue_bulk(1:25, function(x,obs_list) {
  print(names(obs_list[x]))
  run_pmcmc(data = obs_list[[x]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 4,
            lag_rates = 10,
            preyears = 5,
            seasonality_on = 0,
            state_check = 0,
            seasonality_check = 0)},obs_list=tanz_data_list_15to22)
tanz_all_2015to2022_submit$status()#'uncrystallisable_africancivet'
tanz_all_2015to2022_submit$status()#'sharpwitted_mite' <- added incidence <5 yo
tanz_all_2015to2022_submit$tasks$`3a53a5cfb7eb3d51302a5fa94dd8454e`$log()
tanz_all_2015to2022_submit <- obj$task_bundle_get('sharpwitted_mite')
obj$task_bundle_list()
tanz_all_2015to2022_results <- lapply(1:25, function(id){
  tanz_all_2015to2022_submit$tasks[[id]]$result()
})
tanz_all_2015to2022_results[[1]]$history[,1,1]
names(tanz_all_2015to2022_results) <- gsub(' Region','',names(tanz_data_list_15to22))
zoomin_results <- tanz_all_2015to2022_results[c('Lindi','Mtwara','Ruvuma')]
mcmc_sample <- sample(c(201:1000), 100)
lindi_init_eir <- data.frame(init_EIR=exp(zoomin_results$Lindi$mcmc$log_init_EIR[mcmc_sample]),
                             sample_index = mcmc_sample)
lindi_eir_history <- data.frame(cbind(zoomin_results$Lindi$times[-1],t(zoomin_results$Lindi$history['EIR',mcmc_sample,-1])))
names(lindi_eir_history) <- c('t',mcmc_sample)
lindi_eir_history$date <- tanz_data_list_15to22$`Lindi Region`$month

mtwara_init_eir <- data.frame(init_EIR=exp(zoomin_results$Mtwara$mcmc$log_init_EIR[mcmc_sample]),
                             sample_index = mcmc_sample)
mtwara_eir_history <- data.frame(cbind(zoomin_results$Mtwara$times[-1],t(zoomin_results$Mtwara$history['EIR',mcmc_sample,-1])))
names(mtwara_eir_history) <- c('t',mcmc_sample)
mtwara_eir_history$date <- tanz_data_list_15to22$`Mtwara Region`$month

ruvuma_init_eir <- data.frame(init_EIR=exp(zoomin_results$Ruvuma$mcmc$log_init_EIR[mcmc_sample]),
                             sample_index = mcmc_sample)
ruvuma_eir_history <- data.frame(cbind(zoomin_results$Ruvuma$times[-1],t(zoomin_results$Ruvuma$history['EIR',mcmc_sample,-1])))
names(ruvuma_eir_history) <- c('t',mcmc_sample)
ruvuma_eir_history$date <- tanz_data_list_15to22$`Ruvuma Region`$month

saveRDS(list(init_EIR=lindi_init_eir,EIR_history=lindi_eir_history),'./tanz/lindi_eirs_4pw.rds')
saveRDS(list(init_EIR=mtwara_init_eir,EIR_history=mtwara_eir_history),'./tanz/mtwara_eirs_4pw.rds')
saveRDS(list(init_EIR=ruvuma_init_eir,EIR_history=ruvuma_eir_history),'./tanz/ruvuma_eirs_4pw.rds')

str(t(zoomin_results$Lindi$history['EIR',501:1000,-1]))
##Get den
run_pmcmc(data_raw = tanz_data_list_15to22[[1]],
          n_particles = 200,
          proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
          max_EIR=1000,
          max_steps = 1e7,
          atol = 1e-5,
          rtol = 1e-6,
          n_steps = 5,
          n_threads = 1,
          lag_rates = 10,
          preyears = 5,
          seasonality_on = 0,
          state_check = 1,
          seasonality_check = 0)
init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
prop_treated <- 0.4
het_brackets <- 5
mpl <- model_param_list_create(init_age = init_age,
                                  pro_treated = prop_treated,
                                  het_brackets = het_brackets)
state <- equilibrium_init_create_stripped(age_vector = mpl$init_age,
                                          init_EIR = 10,
                                          ft = 0.4,
                                          model_param_list = mpl,
                                          het_brackets = het_brackets)
state$den
plot(init_age,state$den)
sum(state$den)
age <- init_age * mpl$DY
na <- as.integer(length(age))  # number of age groups

age_rate <- age_width <- age_mid_point <- den <- c()
for (i in 1:(na-1))
{
  age_width[i] <- age[i+1] - age[i]
  age_rate[i] <- 1/(age[i + 1] - age[i])  # vector of rates at which people leave each age group (1/age group width)
  age_mid_point[i] <- 0.5 * (age[i] + age[i + 1])  # set age group vector to the midpoint of the group
  
}
age_rate[na] = 0
den <- 1/(1 + age_rate[1]/mpl$eta)
for (i in 1:(na-1))
{
  den[i+1] <- age_rate[i] * den[i]/(age_rate[i+1] + mpl$eta)  # proportion in each age_vector group
}
age_mid_point[length(age_mid_point)+1] <- 90*mpl$DY
plot(age_mid_point,den)
state$den
age_df <- data.frame(age_dy = age_mid_point,
                     age_yr = age_mid_point/mpl$DY,
                     init_age = init_age,
                     prop = den)
age_df <-  age_df %>%
  mutate(age_cat = case_when(
    init_age < 1 ~ '<1',
    init_age >= 1 & init_age < 5 ~ '1-5',
    init_age >= 5 & init_age < 10 ~ '5-10',
    init_age >= 10 & init_age < 20 ~ '10-20',
    init_age >= 20 & init_age < 30 ~ '20-30',
    init_age >= 30 & init_age < 40 ~ '30-40',
    init_age >= 40 & init_age < 50 ~ '40-50',
    init_age >= 50 & init_age < 60 ~ '50-60',
    init_age >= 60 & init_age < 70 ~ '60-70',
    init_age >= 70 & init_age < 80 ~ '70-80',
    init_age >= 80  ~ '80+',
    .default = NA
  ),
  age_cat_2 = factor(case_when(
    init_age < 1 ~ '<1',
    init_age >= 1 & init_age < 5 ~ '1-4',
    init_age >= 5 & init_age < 10 ~ '5-9',
    init_age >= 10 & init_age < 15 ~ '10-14',
    init_age >= 15 & init_age < 50 ~ '15-49',
    init_age >= 50 & init_age < 60 ~ '50-59',
    init_age >= 60 ~ '60+',
    .default = NA
  ),levels = c('<1','1-4','5-9','10-14','15-49','50-59','60+')))

age_cat_df <- age_df %>%
  group_by(age_cat)%>%
  summarise(age_prop = sum(prop))%>%
  mutate(age_factor = factor(age_cat,levels = c('<1','1-5','5-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80+')))
age_cat_2_df <- age_df %>%
  group_by(age_cat_2)%>%
  summarise(age_prop = sum(prop))
  
plot(age_cat_df$age_cat,age_cat_df$age_prop)
ggplot(age_cat_2_df,aes(x=age_cat_2, y=age_prop))+
  geom_bar(stat="identity")+
  scale_y_continuous(limits=c(0,0.6))
