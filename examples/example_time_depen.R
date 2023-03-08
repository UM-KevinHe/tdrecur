## simulation
# basic information --------------------------------------------------------

require(mvtnorm) # used in simulating z
require(survival)
require(dplyr)

source("R/generate_data_indep.R")
surv_dat_dep = generate_data()


# result --------------------------------------------------------------------

require(Rcpp)
sourceCpp("src/main_dep.cpp")
facility_result = compute_facility_idx(surv_dat_dep$facility, surv_dat_dep$f)
dresult = compute_d12(hosp_begin, max_d = O, D = ndays, z, Z_tv, t_tv, max_v,
                      id, facility, num_facility = f)
pp = p + 2
ddloglik_cpp6(L1 = rep(0, pp), L2 = matrix(0, pp, pp), t_start = T1+1, t_end = Tex,
              dresult$Sm, dresult$events_per_day_facility, id,
              z, Z_tv, t_tv, max_v, rep(0, pp), N, facility_result)
Cresult14 = cox_breslow_shr12(hosp_begin, max_d = O, t_start = T1+1, t_end = Tex, z, Z_tv, t_tv, max_v,
                            1e-6, 1e-6, 20,
                            facility, num_facility = f, D = ndays, beta = rep(0, p + 2))#, nthreads = 2)

