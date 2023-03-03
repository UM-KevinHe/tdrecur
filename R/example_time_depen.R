rm(list = ls())
t0 = proc.time()
## simulation
# basic information --------------------------------------------------------

require(mvtnorm) # used in simulating z
require(survival)
require(dplyr)
f = 8#000     ###number of facility
F_pre = 0:(f - 1)
ndays = 36#5
truebeta=c(0.2, -0.1, 3.1, 0.9, 0.5)
p = length(truebeta)
Sigma_z1=diag(p) * 0.01
prob = seq(from = 0.2, to = 0.52, length.out = p)
n = 125
n1 = 125
loop = 106

set.seed(loop)
n_f = c(n1, rpois(f - 1, lambda = n))  #sample size for each facility
N = sum(n_f)                    #number of people in total

T1 = floor(runif(N, 0, (ndays)))       # random start time
T1[sample(1:N, round(N * 0.8))] = 0    # 80% start from 0
t = ceiling(runif(N, 0, ndays)) # time at risk
t[sample(1:N, round(N * 0.77))] = ndays
t[t > ndays - T1] = ndays - T1[t > ndays - T1]
Tex = t + T1
# imblc= sapply(1:ndays, FUN = function(t){sum(T1 < t & Tex >= t)})
# plot(1:ndays, imblc, main = "0.77")

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

z= rmvnorm(N, mean=rep(0,p), sigma=Sigma_z1) # z~N(0,I)
nh = matrix(rbinom(N * p, 1, prob), nrow = N, byrow = T)

alpha = rnorm(f, 0, 0.2)        # facility effect
alpha_subject=rep(alpha, n_f)   #extend alpha to each person
facility = rep(F_pre, n_f)      # facility index for each person, start from 0

# generate T
tempz1 = z_tv1[, 1]
tempz2 = z_tv2[, 1]


O = rep(-1, N) # number of observed times of hospitalization
X = T1 # time flag
i = 0
hosp_begin = NULL

while (sum(O == -1) > 0) {
  exp_zbeta1 <- exp(z %*% truebeta + nh %*% truebeta + tempz1 * 1 + tempz2 * 1 + alpha_subject)
  print(sum(O == -1))
  flag = X > 120 & X <= 302  # May to Oct.

  truelambda1 = 2e-5 * exp_zbeta1  # more events in winter
  truelambda1[flag] = 1e-5 * exp_zbeta1[flag]   # lambda(t) = 0.002
  U = runif(N, 0, 1)
  X = (ceiling((-log(U) / (truelambda1))) + X) # first event

  hosp_begin = cbind(hosp_begin, X)
  O[X > Tex & O == -1] = i
  i = i + 1

  # update z
  for(j in 1:4){
    tempz1[(X > t_tv1[, j]) & (!is.na(t_tv1[, j]))] = z_tv1[(X > t_tv1[, j]) & (!is.na(t_tv1[, j])), (j + 1)]
    tempz2[(X > t_tv2[, j]) & (!is.na(t_tv2[, j]))] = z_tv2[(X > t_tv2[, j]) & (!is.na(t_tv2[, j])), (j + 1)]
  }
}
hosp_begin[hosp_begin > Tex] = NA
max_d = ncol(hosp_begin)
summary(hosp_begin)

sum(hosp_begin, na.rm = T) / N

hosp_end = hosp_begin + 5
hosp_end[hosp_end > Tex] = NA
id = 0:(N - 1) # minus 1 so that the index starts from 0

t1 = proc.time() - t0
# result --------------------------------------------------------------------

require(Rcpp)
sourceCpp("src/shr_dep.cpp")
facility_result = compute_facility_idx(facility, f)
dresult = compute_d12(hosp_begin, max_d = O, D = ndays, z, Z_tv, t_tv, max_v,
                      id, facility, num_facility = f)
pp = p + 2
ddloglik_cpp6(L1 = rep(0, pp), L2 = matrix(0, pp, pp), t_start = T1+1, t_end = Tex,
              dresult$Sm, dresult$events_per_day_facility, id,
              z, Z_tv, t_tv, max_v, rep(0, pp), N, facility_result)
Cresult14 = cox_breslow_shr12(hosp_begin, max_d = O, t_start = T1+1, t_end = Tex, z, Z_tv, t_tv, max_v,
                            1e-6, 1e-6, 20,
                            facility, num_facility = f, D = ndays, beta = rep(0, p + 2))#, nthreads = 2)

