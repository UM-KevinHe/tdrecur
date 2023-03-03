#' Check Whether the Input is Valid
#'
#' Check absence of variables, missingness of variables, variation in covariates, correlation, and VIF.
#' If we fit a linear model of explanatory variable j vs. all other explanatory variables,
#' and calculate the model's R-square, then the VIF of variable j is an increasing function of
#' this R-square: $VIF_j=1/(1-R_j^2)$. In other words, VIF evaluates how much one explanatory variable can be
#'  explained by other explanatory variables in a linear relationship. The more a variable is linearly
#'  related to other variables, the higher the VIF of this variable. Usually, a VIF larger than 10 is
#'  considered as having high collinearity.

#' @param       data a data frame including response, provider ID, and covariates, with missing values imputed
#' @param    T1.char a character string as name of time to enter the study
#' @param    T2.char a character string as name of time to exit the study
#' @param    id.char a character string as name of patient ID
#' @param event.prefix a character string as the prefix of name of event time variable
#' @param  event.num maximum number of events for one patient
#' @param event.char a character string as the prefix of name of event time variable
#' @param nevents.char a character string as name of number of events of each patient
#' @param   Zti.char a vector of character strings as names of time-independent covariates
#' @param   Ztd.char a vector of character strings as names of time-dependent covariates
#' @param   ttd.char a vector of character strings as names of change points of time-dependent covariates
#' @param  prov.char a character string as name of variable consisting of provider IDs
#' @param         TD a Boolean indicating whether there are time-dependent covariates, default true
#' @param threshold.cor minimum pairwise correlation among time-independent covariates to trigger a warning
#' @param     cutoff an integer as cutoff of provider size with 10 as default. A provider will be filtered out if its
#' size is smaller or equal to cutoff
#' @param      check a Boolean (default FALSE) indicating whether checks of variation, correlation, collinearity are
#' needed
#' @return A data frame with non-missing patients and non-small facilities only and additional variables 'no.readm', which
#' equals 1 if the block has no events. Notice that the block effects of those with no events are truncated at -10.

fe.data.prep <- function(data, T1.char, T2.char, event.prefix = NULL, event.num = NULL, event.char = NULL,
                         nevents.char = NULL, Zti.char, Ztd.char = NULL, ttd.char = NULL,
                         prov.char, TD = TRUE, threshold.cor = 0.9, cutoff=10, check=FALSE) {
  #       data: a data frame including response, provider ID, and
  #             covariates, with missing values imputed
  #  prov.char: a character string as name of variable consisting of provider IDs
  #     cutoff: an integer as cutoff of provider size with 10 as default

    ## check absence of variables
    message("Checking absence of variables ... ")
    T1.ind <- match(T1.char, names(data))
    if (is.na(T1.ind)) {
      stop(paste("Time to start '", T1.char, "' NOT found!", sep=""),call.=F)
    }
    T2.ind <- match(T2.char, names(data))
    if (is.na(T2.ind)) {
      stop(paste("Time to exit '", T2.char, "' NOT found!", sep=""),call.=F)
    }
    # id.ind <- match(id.char, names(data))
    # if (is.na(id.ind)) {
    #   stop(paste("Patient ID '", id.char, "' NOT found!", sep=""),call.=F)
    # }
    if (length(event.char) == 0){     # no direct input of event.char
      if (length(event.prefix) == 0 | length(event.num) == 0){ # no input of event.prefix
        stop("Input for event time NOT complete", call.=F)
      }else{
        event.char = paste0(event.prefix, 1:event.num)
      }
    }else{
      if (length(event.prefix) > 0 & length(event.num) > 0){ # complete input of event.prefix
        event.char2 = paste0(event.prefix, 1:event.num)
        if (!all.equal(event.char, event.char2, check.attributes = F)){
          warning("event.char does NOT match string produced by event.prefix and event.num, the latter is ignored!",
                  immediate.=T,call.=F)
        }
      }
    }
    event.ind <- match(event.char, names(data))
    if (sum(is.na(event.ind))>0) {
      stop(paste("Event time variable '", paste(event.char[is.na(event.ind)], collapse="', '"), "' NOT found!", sep=""),call.=F)
    }
    if (length(nevents.char) == 0){     # no direct input of nevents.char
      nevents.char = "nevents"
      data[, nevents.char] = rowSums(data[, event.char], na.rm = T)
    }else{
      if(sum(data[, nevents.char] != rowSums(data[, event.char]), na.rm = T) >0){
        warning("Number of events in nevents.char does NOT match number produced by summing up events, the previous is ignored!",
                immediate.=T,call.=F)
      }
    }
    Zti.ind <- match(Zti.char, names(data))
    if (sum(is.na(Zti.ind)) > 0) {
      stop(paste("Covariate(s) '", paste(Zti.char[is.na(Zti.ind)], collapse="', '"), "' NOT found!", sep=""),call.=F)
    }
    if(TD){
      Ztd.ind <- match(Ztd.char, names(data))
      if (sum(is.na(Ztd.ind)) > 0) {
        stop(paste("Covariate(s) '", paste(Ztd.char[is.na(Ztd.ind)], collapse="', '"), "' NOT found!", sep=""),call.=F)
      }
      ttd.ind <- match(ttd.char, names(data))
      if (sum(is.na(ttd.ind)) > 0) {
        stop(paste("Changing point variable(s) '", paste(ttd.char[is.na(ttd.ind)], collapse="', '"), "' NOT found!", sep=""),call.=F)
      }
    }
    prov.ind <- match(prov.char, names(data))
    if (is.na(prov.ind)) {
      stop(paste("Provider ID '", prov.char, "' NOT found!", sep=""),call.=F)
    }
    message("Checking absence of variables completed!")

    ## check missingness of variables
    message("Checking missingness of variables ... ")
    complete_ind = complete.cases(data[,c(T1.char,T2.char,Zti.char,prov.char)])
    if (sum(complete_ind)==NROW(data)) {
      message("Missing values NOT found. Checking missingness of variables completed!")
    } else {
      check.na <- function(name) {
        if (sum(is.na(data[,name])) > 0) {
          warning(sum(is.na(data[,name]))," out of ",NROW(data[,name])," in '",name,"' missing!",immediate.=T,call.=F)
        }
      }
      invisible(sapply(c(T1.char,T2.char,Zti.char,prov.char), check.na))
      missingness <- (1 - sum(complete_ind) / NROW(data)) * 100
      warning(paste(round(missingness,2), "% of all observations are missing and will be ignored!",sep=""),immediate.=T,call.=F)
    }
    complete_ind = which(complete_ind)

  if (check) {
    ## check variation in time-independent covariates
    message("Checking variation in time-independent covariates ... ")
    time0 = proc.time()
    nzv <- caret::nearZeroVar(data[complete_ind,Zti.char], saveMetrics=T)
    if (sum(nzv$zeroVar==T) > 0) {
      stop("Covariate(s) '", paste(row.names(nzv[nzv$zeroVar==T,]), collapse="', '"),
           "' with zero variance(s)!", call.=F)
    } else if (sum(nzv$nzv==T) > 0) {
      warning("Covariate(s) '",paste(row.names(nzv[nzv$nzv==T,]), collapse="', '"),
              "' with near zero variance(s)!",immediate.=T,call.=F)
    }
    time1 = proc.time() - time0
    message(paste0("Checking variation in time-independent covariates completed in ", round(time1[3]),
                   " seconds!"))

    ## check correlation among time-independent covariates
    message("Checking pairwise correlation among time-independent covariates ... ")
    time0 = proc.time()
    cor <- cor(data[complete_ind,Zti.char])
    if (sum(abs(cor[upper.tri(cor)])>threshold.cor) > 0) {
      cor[lower.tri(cor,diag=T)] <- 0
      ind <- which(abs(cor)>threshold.cor)
      pairs <- sapply(ind, function(ind) c(rownames(cor)[ind%%NROW(cor)],
                                           colnames(cor)[ind%/%NROW(cor)+1]))
      warning("The following ", NCOL(pairs),
              " pair(s) of covariates are highly correlated (correlation > ",
              threshold.cor,"): ", immediate.=T, call.=F)
      invisible(apply(pairs,2,function(col) message('("',paste(col, collapse='", "'),'")')))
    }
    time1 = proc.time() - time0
    message(paste0("Checking pairwise correlation among time-independent covariates completed in ", round(time1[3]),
            " seconds!"))

    ## check collinearity of time-independent covariates
    message("Checking collinearity of time-independent covariates ... ")
    time0 = proc.time()
    # note: the outcome does not matter, so I just put in a non-missing variable
    m.lm <- lm(as.formula(paste(T1.char,"~",paste(Zti.char, collapse="+"))), data=data[complete_ind,])
    vif <- olsrr::ols_vif_tol(m.lm)
    if(sum(vif$VIF >= 10) > 0){
      warning("Covariate(s) '",
              paste(as.data.frame(vif)[vif$VIF>=10,"Variables"], collapse="', '"),
              "' are highly correlated with other covariates!",immediate.=T,call.=F)
    }
    time1 = proc.time() - time0
    message(paste0("Checking collinearity of time-independent covariates completed in ", round(time1[3]/60),
            " minutes!"))
  }

  ## remove small facilities
  prov.size <- as.integer(table(data[complete_ind,prov.char])) # provider sizes
  if(min(prov.size) <= cutoff){
    prov.size.long <- rep(prov.size,prov.size) # provider sizes assigned to patients
    not_small <- (prov.size.long > cutoff) # boolean of whether patients are not in small facilities
    warning(sum(prov.size<=cutoff)," out of ",length(prov.size),
            " providers considered small and filtered out!",immediate.=T,call.=F)
  }else{
    not_small = TRUE
  }
  prov.list <- unique(unlist(data[complete_ind[not_small],prov.char]))   # a reduced list of provider IDs

  # providers with no events
  prov.no.readm <-      # providers with no events
    prov.list[sapply(split(data[complete_ind[not_small],nevents.char],
                           data[complete_ind[not_small],prov.char]),sum)==0]
  data$no.readm <- 0
  data$no.readm[unlist(data[complete_ind[not_small],prov.char])%in%c(prov.no.readm)] <- 1
  message(paste(length(prov.no.readm),"out of",length(prov.list),
                "remaining providers with no events."))

  return(list(data[complete_ind[not_small], ], event.char, nevents.char))
}
#' Fit Recurrent Events Model
#'
#' Fit a recurrent events model with different options
#'
#' @param       data a data frame including response, provider ID, and covariates, with missing values imputed
#' @param    T1.char a character string as name of time to enter the study
#' @param    T2.char a character string as name of time to exit the study
#' @param    id.char a character string as name of patient ID
#' @param event.char a character string as the prefix of name of event time variable
#' @param nevents.char a character string as name of number of events of each patient
#' @param   Zti.char a vector of character strings as names of time-independent covariates
#' @param   Ztd.char a vector of character strings as names of time-dependent covariates
#' @param   ttd.char a vector of character strings as names of change points of time-dependent covariates
#' @param  prov.char a character string as name of variable consisting of provider IDs
#' @param       tolb a small positive number specifying stopping criterion of Newton-Raphson algorithm for beta: covariate effect
#' @param       tola a small positive number specifying stopping criterion of Newton-Raphson algorithm for alpha: block effect
#' @param      bound the maximum number that block effect can take
#' @param   max.iter maximum number of iterations of beta fitting
#' @param  backtrack a boolean indicating whether backtracking line search is implemented, defaulting to FALSE
#' @param  TimeDepen a Boolean indicating whether there are time-dependent covariates, default true
#' @param   parallel a Boolean indicating whether the parallel computing by OpenMP is used, default false
#' @return a data frame with non-missing patients and non-small facilities only and additional variables 'no.readm', which
#' equals 1 if the block has no events. Notice that the block effects of those with no events are truncated at -10
recur.fe.prov <- function(data, T1.char, T2.char, event.char, id.char,
                          nevents.char, Zti.char, Ztd.char = NULL, ttd.char = NULL,
                          prov.char, tolb=1e-6, tola=1e-6, bound = 10.0, max.iter = 10000, backtrack=FALSE, TimeDepen = TRUE, parallel=FALSE,
                          AUC=FALSE){
  # prov.char: a character string as name of variable consisting of provider IDs
  # backtrack: a boolean indicating whether backtracking line search is implemented, defaulting to FALSE

  if (!is.logical(backtrack)) stop("Argument 'backtrack' NOT as required!")
  n.prov <- sapply(split(data[, nevents.char], data[, prov.char]), length) # provider-specific number of patients
  n.readm.prov <- sapply(split(data[, nevents.char], data[, prov.char]), sum) # provider-specific number of events

  if(!TimeDepen){ # the simplest situation: no time dependent variables
    sourceCpp("shr12.cpp")
    Cresult = cox_breslow_ind(data[, event.char], data[, nevents.char], data[, T1.char], data[, T2.char],
                              data[, id.char], data[, Zti.char],
# 2/8 updated through here ----------------------------------------------

                                     1.0e-6, 20,
                                     sresult$facility, f, ndays, beta0)
  }



  Z <- as.matrix(data[,Z.char])
  gamma.prov <- rep(log(mean(data[,Y.char])/(1-mean(data[,Y.char]))), length(n.prov))
  beta <- rep(0, NCOL(Z))
  if (Rcpp) {
    ls <- logis_fe_prov(as.matrix(data[,Y.char]),Z,n.prov,gamma.prov,beta,backtrack,max.iter,bound,tol)
    gamma.prov <- as.numeric(ls$gamma); beta <- as.numeric(ls$beta)
  } else {
    iter <- 0
    beta.crit <- 100 # initialize stop criterion
    message("Implementing Newton-Raphson algorithm for fixed provider effects model ...")
    if (backtrack) {
      s <- 0.01; t <- 0.8 # initialize parameters for backtracking line search
      Loglkd <- function(gamma.obs, beta) {
        sum((gamma.obs+Z%*%beta)*data[,Y.char]-log(1+exp(gamma.obs+Z%*%beta)))
      }
      while (iter<=max.iter & beta.crit>=tol) {
        iter <- iter + 1
        # provider effect update
        gamma.obs <- rep(gamma.prov, n.prov)
        Z.beta <- Z%*%beta
        p <- c(plogis(gamma.obs+Z.beta)); pq <- p*(1-p)
        score.gamma.prov <- sapply(split(data[,Y.char]-p, data[,prov.char]), sum)
        d.gamma.prov <- score.gamma.prov / sapply(split(pq, data[,prov.char]), sum)
        v <- 1 # initialize step size
        loglkd <- Loglkd(rep(gamma.prov, n.prov), beta)
        d.loglkd <- Loglkd(rep(gamma.prov+v*d.gamma.prov, n.prov), beta) - loglkd
        lambda <- score.gamma.prov%*%d.gamma.prov
        while (d.loglkd < s*v*lambda) {
          v <- t * v
          d.loglkd <- Loglkd(rep(gamma.prov+v*d.gamma.prov, n.prov), beta) - loglkd
        }
        gamma.prov <- gamma.prov + v * d.gamma.prov
        gamma.prov <- pmin(pmax(gamma.prov, median(gamma.prov)-bound), median(gamma.prov)+bound)
        gamma.obs <- rep(gamma.prov, n.prov)

        # regression parameter update
        p <- c(plogis(gamma.obs+Z.beta)); pq <- p*(1-p)
        score.beta <- t(Z)%*%(data[,Y.char]-p)
        info.beta <- t(Z)%*%(c(pq)*Z)
        d.beta <- as.numeric(solve(info.beta)%*%score.beta)
        v <- 1 # initialize step size
        loglkd <- Loglkd(gamma.obs, beta)
        d.loglkd <- Loglkd(gamma.obs, beta+v*d.beta) - loglkd
        lambda <- t(score.beta)%*%d.beta
        while (d.loglkd < s*v*lambda) {
          v <- t * v
          d.loglkd <- Loglkd(gamma.obs, beta+v*d.beta) - loglkd
        }
        beta.new <- beta + v * d.beta
        beta.crit <- norm(matrix(beta-beta.new),"I") # stopping criterion
        beta <- beta.new
        cat(paste0("Iter ",iter,": Inf norm of running diff in est reg parm is ",
                   formatC(beta.crit,digits=3,format="e"),";\n"))
      }
    } else {
      while (iter<=max.iter & beta.crit>=tol) {
        iter <- iter + 1
        cat(paste0("\n Iter ",iter,":"))
        # provider effect update
        gamma.obs <- rep(gamma.prov, n.prov)
        Z.beta <- Z%*%beta
        p <- c(plogis(gamma.obs+Z.beta)); pq <- p*(1-p)
        gamma.prov <- sapply(split(data[,Y.char]-p, data[,prov.char]), sum) /
          sapply(split(pq, data[,prov.char]), sum) + gamma.prov
        gamma.prov <- pmin(pmax(gamma.prov, median(gamma.prov)-bound), median(gamma.prov)+bound)
        gamma.obs <- rep(gamma.prov, n.prov)

        # regression parameter update
        p <- c(plogis(gamma.obs+Z.beta)); pq <- p*(1-p)
        score.beta <- t(Z)%*%(data[,Y.char]-p)
        info.beta <- t(Z)%*%(c(pq)*Z)
        beta.new <- beta + as.numeric(solve(info.beta)%*%score.beta)
        beta.crit <- norm(matrix(beta-beta.new),"I") # stopping criterion
        beta <- beta.new
        cat(paste0(" Inf norm of running diff in est reg parm is ",formatC(beta.crit,digits=3,format="e"),";"))
      }

    }
    message("\n Newton-Raphson algorithm converged after ",iter," iterations!")
  }

  gamma.obs <- rep(gamma.prov, n.prov)
  gamma.null <- ifelse(null=="median", median(gamma.prov),
                       ifelse(class(null)=="numeric", null[1],
                              stop("Argument 'null' NOT as required!",call.=F)))
  Exp <- as.numeric(plogis(gamma.null+Z%*%beta)) # expected prob of readm within 30 days of discharge under null
  Pred <- as.numeric(plogis(gamma.obs+Z%*%beta))
  df.prov <- data.frame(Obs=sapply(split(data[,Y.char],data[,prov.char]),sum),
                        Exp=sapply(split(Exp,data[,prov.char]),sum))
  df.prov$SRR <- df.prov$Obs / df.prov$Exp
  df.prov$gamma <- gamma.prov
  neg2Loglkd <- -2*sum((gamma.obs+Z%*%beta)*data[,Y.char]-log(1+exp(gamma.obs+Z%*%beta)))
  AIC <- neg2Loglkd + 2 * (length(gamma.prov)+length(beta))
  BIC <- neg2Loglkd + log(nrow(data)) * (length(gamma.prov)+length(beta))
  gamma.prov[n.readm.prov==n.prov] <- Inf; gamma.prov[n.readm.prov==0] <- -Inf
  if (AUC) {
    AUC <- pROC::auc(data[,Y.char], Pred)
    return(list(beta=beta, Obs=data[, Y.char], Exp=Exp, df.prov=df.prov,
                neg2Loglkd=neg2Loglkd, AIC=AIC, BIC=BIC, AUC=AUC[1]))
  } else {
    return(list(beta=beta, Obs=data[, Y.char], Exp=Exp, df.prov=df.prov,
                neg2Loglkd=neg2Loglkd, AIC=AIC, BIC=BIC))
  }
  #       beta: a vector of fixed effect estimates
  #        Obs: a vector of responses for included providers
  #        Exp: a vector of expected probs of readmission within 30 days of discharge
  #    df.prov: a data frame of provider-level number of observed number of readmissions within 30 days
  #             expected number of readmissions within 30 days, a vector of SRRs, and a vector of
  #             provider effect estimates for included providers (considered as a fixed effect)
  # neg2Loglkd: minus two times log likelihood
  #        AIC: Akaike info criterion
  #        BIC: Bayesian info criterion
  #        AUC: area under the ROC curve
} # end of logis.fe.prov

