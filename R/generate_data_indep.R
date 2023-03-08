#'Generate Data with Time-independent Covariates
#'
#'The function generate_data_indep generates example recurrent events data with time-independent covariates only
#'
#'@importFrom mvtnorm rmvnorm

#'@param true_beta true parameters
#'@param F_pre a unique list of facilities
#'@param gamma true facility effects of each facility, notice length of gamma must be the same as length of F_pre
#'@param dmu0 a constant for baseline rate
#'@param ndays number of days followed
#'@param seed random seed
#'@param mean_sample_size expected sample size of each facility (stratum), a scalar same in all facilities
#'@param Sigma_z1 Var-Cov matrix of design matrix
#'@param time_depen indicator of whether the dataset contains time-dependent covariates. This function generates 2 time-dependent covariates for now
#'
#'@return a list of generated data:
#'\item{z}{design matrix}
#'\item{facility}{which facility each patient belongs to}
#'\item{true_beta}{true parameters}
#'\item{t_start}{enter time}
#'\item{t_end}{exit time}
#'\item{num_facility}{number of facilities (strata)}
#'\item{events_time}{observed event time for each patient}
#'\item{nevents}{number of events of each patient}
#'\item{ndays}{number of follow-up days}
#'\item{N}{number of patients}
#'
#'@examples
#'surv_dat = generate_data()
#'
#'@export
#'

generate_data = function(true_beta = c(c(0.5, -0.5, 1, -1, 1.5), rep(0, 5)), F_pre = 0:9, gamma = rep(1, 10),
                         dmu0 = 0.002,
                         ndays = 200, seed = 23021619,
                         mean_sample_size = 50, Sigma_z1 = 0.1 * diag(length(true_beta)),
                         time_depen = TRUE){
  f = length(F_pre)     # number of facilities
  p = length(true_beta)

  set.seed(seed)
  n_f = rpois(f, lambda = mean_sample_size)  #sample size for each facility
  N = sum(n_f)                               #number of people in total


  z= rmvnorm(N, mean=rep(0,p), sigma=Sigma_z1) # z~N(0, 0.1)


  alpha = log(gamma)       # log of facility effect
  alpha_subject=rep(alpha, n_f)   #extend alpha to each person


  facility = rep(F_pre, n_f)        # facility index for each person, start from 0


  # generate T
  T1 = floor(runif(N, 0, (ndays)))       # random start time
  T1[sample(1:N, round(N * 0.8))] = 0    # 80% start from 0
  t = ceiling(runif(N, 0, ndays)) # time at risk
  t[sample(1:N, round(N * 0.75))] = ndays
  t[t > ndays - T1] = ndays - T1[t > ndays - T1]
  Tex = t + T1

  if(time_depen){
    t_tv1 = matrix(sample(1:(2*ndays), 4*N, replace = T), nrow = N)
    t_tv1 = t(apply(t_tv1, 1, FUN = function(x){
      x[duplicated(x)] = ndays + 1
      return(x)}))
    t_tv1 = t(apply(t_tv1, 1, sort))
    t_tv1[t_tv1 > ndays] = NA

    t_tv2 = matrix(sample(1:(2*ndays), 2*N, replace = T), nrow = N)
    t_tv2 = t(apply(t_tv2, 1, FUN = function(x){
      x[duplicated(x)] = ndays + 1
      return(x)}))
    t_tv2 = t(apply(t_tv2, 1, sort))
    t_tv2[t_tv2 > ndays] = NA
    t_tv2 = cbind(t_tv2, NA, NA)

    t_tv = array(c(t_tv1, t_tv2), dim = c(N, 4, 2)) # N patients, 4 change times, and 2 variables

    z_tv1 = matrix(nrow = N, ncol = 4)
    z_tv1[!is.na(t_tv1)] = rnorm(sum(!is.na(t_tv1)))
    z_tv1 = cbind(rnorm(N), z_tv1)
    z_tv2 = cbind(rep(1, N), 0, NA, NA)
    z_tv2[is.na(t_tv2)] = NA
    z_tv2= cbind(0, z_tv2)
    Z_tv = array(c(z_tv1, z_tv2), dim = c(N, 5, 2))

    max_v = cbind(rowSums(!is.na(t_tv1)), rowSums(!is.na(t_tv2)))

    tempz1 = z_tv1[, 1]
    tempz2 = z_tv2[, 1]
  }
  O = rep(-1, N) # number of observed times of hospitalization
  X = T1 # time flag
  i = 0
  hosp_begin = NULL
  if(time_depen){
    while (sum(O == -1) > 0) {
      exp_zbeta1 <- exp(z[, 1:(p-2)] %*% true_beta[1:(p-2)] + tempz1 * true_beta[p-1] +
                          tempz2 * true_beta[p] + alpha_subject)
      truelambda1 = dmu0 * exp_zbeta1
      U = runif(N, 0, 1)
      X = ceiling((-log(U) / (truelambda1))) + X
      hosp_begin = cbind(hosp_begin, X)
      O[X > Tex & O == -1] = i
      i = i + 1

      # update z
      for(j in 1:4){
        tempz1[(X > t_tv1[, j]) & (!is.na(t_tv1[, j]))] = z_tv1[(X > t_tv1[, j]) & (!is.na(t_tv1[, j])), (j + 1)]
        tempz2[(X > t_tv2[, j]) & (!is.na(t_tv2[, j]))] = z_tv2[(X > t_tv2[, j]) & (!is.na(t_tv2[, j])), (j + 1)]
      }
    }
  }else{

    exp_zbeta <- exp(z %*% true_beta + alpha_subject)
    truelambda <- dmu0 * exp_zbeta
    while (sum(O == -1) > 0) {
      U = runif(N, 0, 1)
      X = ceiling((-log(U) / (truelambda))) + X
      hosp_begin = cbind(hosp_begin, X)
      O[X > Tex & O == -1] = i
      i = i + 1
    }
  }


  hosp_begin[hosp_begin > Tex] = NA
  max_d = ncol(hosp_begin)
  summary(hosp_begin)


  if(time_depen){
    simdata = list(z = z[, 1:(p-2)], Z_tv = Z_tv, t_tv = t_tv, facility = facility, num_facility = f, true_beta = true_beta, t_start = T1,
                   t_end = Tex, events_time = hosp_begin, nevents = O, ndays = ndays, N = N)
  }else{
    simdata = list(z = z, facility = facility, num_facility = f, true_beta = true_beta, t_start = T1,
                   t_end = Tex, events_time = hosp_begin, nevents = O, ndays = ndays, N = N)
  }

  return(simdata)
}
