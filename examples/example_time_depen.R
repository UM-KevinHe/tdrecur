## simulation
# basic information --------------------------------------------------------

require(mvtnorm) # used in simulating z
require(survival)
require(dplyr)

source("R/generate_data.R")
surv_dat_dep = generate_data()
p = 10

# result --------------------------------------------------------------------

require(Rcpp)
sourceCpp("src/main_dep.cpp")

Cresult14 = cox_breslow_dep(surv_dat_dep$events_time, surv_dat_dep$nevents, surv_dat_dep$t_start+1,
                            surv_dat_dep$t_end,
                            surv_dat_dep$z, surv_dat_dep$Z_tv, surv_dat_dep$t_tv, surv_dat_dep$max_v, 1.0e-6,
                            1.0e-6,
                            20, surv_dat_dep$facility, surv_dat_dep$num_facility, surv_dat_dep$ndays,
                            beta = rep(0, p),
                            parallel = F)
Cresult = cox_breslow_dep(surv_dat_dep$events_time, surv_dat_dep$nevents, surv_dat_dep$t_start+1,
                            surv_dat_dep$t_end,
                            surv_dat_dep$z, surv_dat_dep$Z_tv, surv_dat_dep$t_tv, surv_dat_dep$max_v, 1.0e-6,
                            1.0e-6,
                            20, surv_dat_dep$facility, surv_dat_dep$num_facility, surv_dat_dep$ndays,
                            beta = rep(0, p),
                            parallel = T, nthreads = 1)
