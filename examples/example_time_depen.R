## simulation
# basic information --------------------------------------------------------

require(mvtnorm) # used in simulating z
require(survival)
require(dplyr)

source("R/generate_data.R")
surv_dat_dep = generate_data()


# result --------------------------------------------------------------------

require(Rcpp)
sourceCpp("src/main_dep.cpp")

Cresult14 = cox_breslow_dep(surv_dat$events_time, surv_dat$nevents, surv_dat$t_start+1, surv_dat$t_end,
                            surv_dat$z, surv_dat$Z_tv, surv_dat$t_tv, surv_dat$, 1.0e-6, 1.0e-6,
                            20, surv_dat$facility, surv_dat$num_facility, surv_dat$ndays, beta = rep(0, p))

