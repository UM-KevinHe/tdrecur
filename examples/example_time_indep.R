require(mvtnorm) # used in simulating z
require(survival)
require(Rcpp)
require(dplyr)

## simulation
source("R/generate_data.R")
surv_dat = generate_data(time_depen = F)

  p = 10
  # fixed-point -------------------------------------------------------------


  require(Rcpp)
  sourceCpp("src/shr_ind.cpp")
  t0 = proc.time()
  Cresult12 = cox_breslow_ind(surv_dat$events_time, surv_dat$nevents, surv_dat$t_start+1, surv_dat$t_end, surv_dat$z, 1.0e-6, 1.0e-6,
                                20, surv_dat$facility, surv_dat$num_facility, surv_dat$ndays, beta = rep(0, p))
  t12 = proc.time() - t0

  # coxph -------------------------------------------------------------------

  simdata = data.frame(surv_dat$facility, surv_dat$t_start, surv_dat$t_end, surv_dat$z, surv_dat$events_time)
  max_d = ncol(surv_dat$events_time)
  colnames(simdata)[1:(3 + p + max_d)] = c("facility","T1","Tex",paste0("z", 1:p), paste0("h", 1:max_d))
  simdata$id=1:surv_dat$N
  Ocox = tmerge(simdata[, c("id", "facility", paste0("z", 1:p))], simdata, id = id,tstart = T1, tstop = Tex)
  for(nevent in 1:(ncol(simdata) - 5-p)){
    Ocox = tmerge(Ocox, simdata, id = id, hosp_begin = event(simdata[, paste0("h", nevent)]))
  }

  # coxph strata ------------------------------------------------------------
  formula_coxph2 = paste("Surv(tstart, tstop, hosp_begin) ~ ", paste(paste0("z", 1:p, " + "), collapse = ""),
                         "strata(facility)")
  coxph_result2 = coxph(as.formula(formula_coxph2), ties = "breslow", data = Ocox)
  coxph_result2$coefficients
Cresult12$beta
