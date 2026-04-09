library(MASS)
library(survsim)
library(expm)
library(truncnorm)
library(tidyverse)
library(nlme)
library(DescTools)
library(survival)
library(tictoc)

SimData3 <- function(n=200, 
                     ZX=c('year','sindoes','doesW', 'trb1'), 
                     YX=NULL, sigma_lme=NULL, 
                     A=c(1.49, 0.74, 0.18, 0.18, 0.21, 0.06, 0.03), 
                     lambda = -5.51, sigma_nlme = 0.62,
                     Vtime=c(0,1,6,12,18,24,30),
                     endTrial=36,
                     alpha=c(-1.65, 0.15, 1.8, -0.05, 0.4,0.15), 
                     beta = c(2, 1, -0.3, 1.5, 0.5), delim_val=1.47, regime=1, #LLOQ number and type
                     rate=NULL, # only when regime = 2
                     Xi=c(0, -2.6), Asso = c( -4, 3.5), #(For surivival model) Xi -> parameter for fixed effect; Asso -> association parameter in random effect 
                     scale=205, shape=5, 
                     censor_type = "right", # or 'interval'
                     ran_p=c(1,1), # number of random effect for each model
                     SIGMA=diag(c(1, 1)) , 
                     cen_par=c(5,1000),
                     check=F,
                     freq=T,  # high frequency of check point
                     contaminated=F, # if there is outlier in lme
                     percent=0 #percent of outlier
                     ) { 
  
  # BEGIN OF THE FUNCTION （in terms of month）
  t <- seq(0, endTrial, by = 1/30)
  maxT <- max(t)        
  nVacc=length(Vtime)  # number of vaccination
  pVacc <- diff(c(Vtime, endTrial))   # (t_{i+1} - t_{i})  
  
  ni <- length(t) # number of time point
  sid <- rep(1:n, each=ni)  # sid: 1, 2, 3, ... (timepoint for each sid)
  time <- rep(t, n)
  base <- rep(rnorm(n, 0, 1), each=ni) # each sid has only 1 # from norm but repeated ni times
  

  dat <- data.frame(sid, time, base)
  dat$injectionNO <- findInterval(dat$time, Vtime, rightmost.closed = FALSE)
  dat$injectionNO[dat$injectionNO == 0] <- 1
  dat$injectionNO[dat$injectionNO > length(pVacc)] <- length(pVacc)
  
  interval_start <- Vtime[dat$injectionNO] # start time
  dat$perd <- pVacc[dat$injectionNO] # length of the interval
  dat$doest <- dat$time - interval_start # time of the interval
  
  
  dat$sindoes <- sin(dat$doest * pi / dat$perd)
  dat$month  <- dat$time
  dat$month2 <- dat$month^2
  dat$year  <- dat$time / 12
  dat$year2 <- dat$year^2
  dat$doesM <- dat$doest
  dat$doesW <- dat$doest * (52 / 12)   # 1 month ≈ 4.33 weeks
  dat$biweek <- dat$time * (52 / 24)   
  N <- nrow(dat)
  
  # -----------------------------------------------------------------------------
  # ------------------------- longitudinal models--------------------------------
  # -----------------------------------------------------------------------------
  
  
  ran_effects <- mvrnorm(n, mu=rep(0, sum(ran_p)), Sigma=SIGMA) # multivariate matrix  [1, 0]
                                                                #                      [0, 1]
  # in first model (NLME/LME)
  dat$trb1 <- rep(ran_effects[,1:ran_p[1]],each=ni) 
  # in GLME
  dat$trb2 <- rep(ran_effects[,(ran_p[1]+1):(sum(ran_p[1:2]))],each=ni)
  
  # NLME and LME should not coexist
  if (!is.null(A) && !is.null(YX)) {
    stop("A and YX cannot both be specified: NLME and LME are mutually exclusive.")
  }
  
  if (is.null(A) && is.null(YX)) {
    stop("Either A (for NLME) or YX (for LME) must be specified.")
  }
  
  if (!is.null(YX) && is.null(sigma_lme)) {
    stop("sigma_lme must be provided when YX is not NULL.")
  }
  
  # (1) generate exponential delay data
  if (!is.null(A)) {
    dose_time <- Vtime
    rate_i <- exp(lambda + dat$trb1) # subject-specific rate
    lag_mat <- sapply(dose_time, function(tk) pmax(dat$time - tk, 0)) #  lag: lag_k = max(time - tk, 0)
    ind_mat <- sapply(dose_time, function(tk) as.numeric(dat$time >= tk)) # # time >= tk
    
    colnames(lag_mat) <- paste0("lag_k", seq_along(dose_time))
    colnames(ind_mat) <- paste0("ind_k", seq_along(dose_time))
    
    dat <- cbind(dat, lag_mat, ind_mat)
    
    contrib_mat <- ind_mat * sweep(exp(-lag_mat * rate_i), 2, A, `*`)  # contrib_mat[,k] = ind_mat[,k] * A[k] * exp(-rate_i * lag_mat[,k])
    mu <- rowSums(contrib_mat)
    dat$y <- mu + rnorm(nrow(dat), 0, sigma_nlme) # error term
    
    # (3) generate longitudinal data with censoring data
  } else if (!is.null(YX)) {
    U <- as.matrix(cbind(1, dat[, YX], dat$trb1))
    dat$y <- U %*% beta + rnorm(N, mean = 0, sd = sigma_lme)
    
    if (contaminated == T) {
      outliers = sample(1:N, N * percent)
      dat$y[outliers] <- dat$y[outliers] + 5 * sd(dat$y)
    }
    if (regime == 0) {   # no censoring
      dat$c <- rep(0, N)
    }
    if (regime == 1) { #right censored with LLOQ
      dat$y[dat$y < delim_val] <- delim_val
      dat$c <- as.numeric(dat$y == delim_val)
    }
    
    if (regime == 2) { # LLOQ in random number from rnorm
      dat$y[dat$y < delim_val] <- delim_val
      if (rate == 0.2) {
        dat$stopT <- rep(rnorm(n, 10 + ran_effects[, 1:ran_p[1]] * 5, 1), each = ni) +
          (dat$injectionNO > 1) * 4.5
      } else if (rate == 0.4) {
        # with 40% censoring, where 34.8% from point mass, 5.2% from normal distribution
        dat$stopT <- rep(rnorm(n, 10, 3), each = ni) +
          (dat$injectionNO > 1) * 3
      } else if (rate == 0.25) {
        # with 40% censoring, where 25% from point mass, 15% from normal distribution
        dat$stopT <- rep(rnorm(n, 10, 3), each = ni) +
          (dat$injectionNO > 1) * 3.5
      } else {
        stop("For regime = 2, rate must be one of 0.2, 0.25, or 0.4.")
      }
      
      dat$y[dat$doest > dat$stopT] <- delim_val
      dat$c <- as.numeric(dat$y == delim_val) # mark observation outside LLOQ as delim_val
    }
  }
  
  # (2) generate binomial data
  X <- as.matrix(cbind(1, dat[, ZX], dat$trb2 * dat$month))
  logit.z <- X %*% alpha
  dat$p.z <- exp(logit.z) / (1 + exp(logit.z))
  dat$z <- rbinom(N, 1, prob = dat$p.z)
  
  
  
  
  
  
  # -----------------------------------------------------------------------------
  # ------------------------ generate survival time -----------------------------
  # -----------------------------------------------------------------------------
  
  if (!(censor_type %in% c("right", "interval"))) {
    stop("censor_type must be either 'right' or 'interval'.")
  }
  
  Xs <- cbind(1, dat[!duplicated(dat$sid, fromLast = TRUE), "base"])
  sexp <- exp(Xs %*% Xi + ran_effects %*% Asso)
  scale1 <- scale * sexp^(-1 / shape)
  
  S <- rweibull(n, shape = shape, scale = scale1)
  cen <- rweibull(n, shape = cen_par[1], scale = cen_par[2])
  
  visit_grid <- seq(0, maxT, by = 0.5)
  
  if (censor_type == "right") {
    
    event_time <- pmin(S, cen, maxT)
    event <- as.numeric(S <= cen & S <= maxT)
    
    surv.dat <- data.frame(
      sid = 1:n,
      event_time = event_time,
      event = event
    )
    
    # keep longitudinal observations before observed time
    dat <- merge(dat, surv.dat[, c("sid", "event_time")], by = "sid")
    dat <- dat[dat$time <= dat$event_time, ]
    
    # keep a unified cutoff variable for later subsampling
    dat$cutoff_time <- dat$event_time
    dat$event_time <- NULL
    
  } else if (censor_type == "interval") {
    
    L <- R <- rep(NA_real_, n)
    event <- rep(NA_real_, n)
    
    for (i in 1:n) {
      if (S[i] <= cen[i] && S[i] <= maxT) {
        # interval-censored event
        L[i] <- max(visit_grid[visit_grid < S[i]])
        R[i] <- min(visit_grid[visit_grid >= S[i]])
        event[i] <- 1
      } else {
        # right-censored by censoring or end of study
        t_cens <- min(cen[i], maxT)
        L[i] <- max(visit_grid[visit_grid <= t_cens])
        R[i] <- Inf
        event[i] <- 0
      }
    }
    
    surv.dat <- data.frame(
      sid = 1:n,
      L = L,
      R = R,
      event = event
    )
    
    # keep longitudinal observations before left endpoint
    dat <- merge(dat, surv.dat[, c("sid", "L")], by = "sid")
    dat <- dat[dat$time <= dat$L, ]
    
    # keep a unified cutoff variable for later subsampling
    dat$cutoff_time <- dat$L
    dat$L <- NULL
  }
  
  
  # -----------------------------------------------------------------------------
  # ----------------------------- subsampling -----------------------------------
  # -----------------------------------------------------------------------------
  
  ids <- unique(dat$sid)
  subsam_list <- vector("list", length(ids))
  
  for (i in seq_along(ids)) {
    subdat <- subset(dat, sid == ids[i])
    
    # subject-specific observed longitudinal end time
    maxt <- max(subdat$cutoff_time)
    
    if (freq == TRUE) {
      rtimes <- seq(0, maxt, by = 0.5)
    } else {
      rtimes <- c(0, 1, 6, 12, 18, 24, 30, 36)
      rtimes <- rtimes[rtimes <= maxt]
    }
    
    rtimes <- rtimes + runif(length(rtimes), -2/30, 2/30)
    rtimes <- pmax(0, pmin(rtimes, maxt))
    keept <- sort(unique(c(0, rtimes, maxt)))
    
    idx_keep <- sapply(keept, function(tt) {
      which.min(abs(subdat$time - tt))
    })
    
    idx_keep <- sort(unique(idx_keep))
    subsam_list[[i]] <- subdat[idx_keep, ]
  }
  
  finalDat <- do.call(rbind, subsam_list)
  row.names(finalDat) <- NULL
  
  # remove unified cutoff variable from longitudinal dataset
  finalDat$cutoff_time <- NULL
  
  base_df <- unique(finalDat[, c("sid", "base")])
  finalSdat <- merge(surv.dat, base_df, by = "sid", all.x = TRUE)
  
  dat <- finalDat
  
  # check the result
  if(check==T){
    hist(dat$p.z)       
    cat("mean of z:",mean(dat$z), '\n')
    hist(dat$y)
    if ("c" %in% names(dat)) {
      cat("mean of c", mean(dat$c), '\n')
    }
    
    #    par(mfrow=c(1,2))
    hist(S, main="true survival time")
    # hist(surv.dat$obs_time, main="observed event time")
    if (censor_type == "right") {
      hist(surv.dat$event_time[surv.dat$event == 1], main = "observed event time of interest")
      hist(surv.dat$event_time[surv.dat$event == 0], main = "observed censoring time")
    }
    
    if (censor_type == "interval") {
      hist(surv.dat$L[surv.dat$event == 1], main = "left endpoint of interval")
      hist(surv.dat$R[is.finite(surv.dat$R) & surv.dat$event == 1], main = "right endpoint of interval")
    }
  }
  return(list(dat=finalDat, surv.dat=finalSdat, ran_effects))}
  
  
  
run_one_sim <- function(seed = NULL, save_dir = "E:/UBC/Final Project/New folder/Final-Project/orin/HHJMs-master/example/simulationed_data_group") {
  if (!is.null(seed)) set.seed(seed)
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
 Sim_1 <-  SimData3(n=200, 
                      ZX=c('month','sindoes','doesW', 'trb1'), 
                      YX=NULL, sigma_lme=NULL, 
                      A=c(1.49, 0.74, 0.18, 0.18, 0.21, 0.06, 0.03), 
                      lambda = -5.51, sigma_nlme = 0.62,
                      Vtime=c(0,1,6,12,18,24,30),
                      endTrial=36,
                      alpha=c(-1.65, 0.11, 1.6, -0.02, 0.17, 0.09), 
                      beta = NULL, delim_val=1.47, regime=0, #LLOQ number and type
                      rate=NULL, # only when regime = 2
                      Xi=c(0, -0.26), Asso = c( 0.76, 0.7), #(For surivival model) Xi -> parameter for fixed effect; Asso -> association parameter in random effect 
                      scale=45, shape=5, 
                      censor_type = "right", # or 'interval'
                      ran_p=c(1,1), # number of random effect for each model
                      SIGMA=diag(c(1, 1)) , 
                      cen_par=c(5,30),
                      check=F,
                      freq=T,  # high frequency of check point
                      contaminated=F, # if there is outlier in lme
                      percent=0 #percent of outlier
 ) 
 
 
 Sim_1$dat
 
 
 # (1) Model 1: NLME of y
 fm1_sim <- y ~ 
   A1 * ind_k1 * exp(-lambda * lag_k1) +
   A2 * ind_k2 * exp(-lambda * lag_k2) +
   A3 * ind_k3 * exp(-lambda * lag_k3) +
   A4 * ind_k4 * exp(-lambda * lag_k4) +
   A5 * ind_k5 * exp(-lambda * lag_k5) +
   A6 * ind_k6 * exp(-lambda * lag_k6) +
   A7 * ind_k7 * exp(-lambda * lag_k7)
 
 
 start_nls <- list(
   A1=.5,A2=.6,A3=.6,A4=.6,A5=.5,A6=.5,A7=.5,
   lambda= 1
 )
 
 fit_nls <- nls(
   fm1_sim, data=Sim_1$dat, start=start_nls, algorithm="port",
   lower=c(A1=0,A2=0,A3=0,A4=0,A5=0,A6=0,A7=0, lambda=0.15),
   upper=c(A1=Inf,A2=Inf,A3=Inf,A4=Inf,A5=Inf,A6=Inf,A7=Inf, lambda=2),
   control=nls.control(maxiter=500, tol=1e-4, minFactor=1/2048)
 )
 
 
 
 fit_nlme_sim <- nlme(
   model = fm1_sim,
   data  = groupedData(y ~ month | sid, data = Sim_1$dat),
   fixed = A1 + A2 + A3 + A4 + A5 + A6 + A7 + lambda ~ 1,
   random =  list(sid = pdDiag(lambda  ~ 1)),
   start = as.numeric(coef(fit_nls)),
   control = nlmeControl(
     maxIter=2000, pnlsMaxIter=3000, msMaxIter=2000,
     tolerance = 1e-5, pnlsTol = 1e-5)
 )
 
 
 # get the estimated random effect by observation
 b_11_sim <- as.numeric(scale( ranef(fit_nlme_sim)[, "lambda"], center = TRUE, scale = TRUE))
 
 b_11_add_sim <- data.frame(
   sid = rownames(ranef(fit_nlme_sim) ),
   b_11 = b_11_sim)
 mydat_nlme_sim <- merge(Sim_1$dat, b_11_add_sim , by = "sid", all.x = TRUE)
 
 # (2) Model 2: a GLME model of C
 ## using the estimated random intercept from model 1 (i.e. estb11, scaled) as a covariate
 fm3_sim <- z ~ 1+month+sindoes+doesW+b_11+(month-1|sid)
 md3_sim <- glmer(fm3_sim, family="binomial", data=mydat_nlme_sim)
 
 # (3) Model 3: a survival model
 Sdata_nlme_sim <- Sim_1$surv.dat
 Sdata_nlme_sim$b_11 <- scale(ranef(fit_nlme_sim)[, "lambda"], center = TRUE, scale = TRUE)
 Sdata_nlme_sim$b_21 <- scale(ranef(md3_sim)$sid[,1], center=T, scale=T)
 
 
 # a Cox PH model
 fitCOX1_sim <- coxph(Surv(event_time, event) ~ base+ b_11 + b_21, data =  Sdata_nlme_sim)
 
 
 ###############################
 ######   Joint modeling   #####
 ###############################
 
 
 nlmeObject_sim <- list(
   type = "nlme",
   family = "EXP_Delay",
   fm = fm1_sim,
   sigma = "sigma",
   disp = "b_11",
   A_names = paste0("A",1:7),
   ind_names = paste0("ind_k",1:7),
   lag_names = paste0("lag_k",1:7),
   lambda_name = "lambda",
   ran_map = list(
     lambda = "b_11"
   ),
   str_val = c(as.list(fixef(fit_nlme_sim)),  sd(ranef(fit_nlme_sim)[["lambda"]])),      # A1..A7 + lambda
   CenObject = NULL
 )
 
 lambda_hat <- as.numeric(fixef(fit_nlme_sim)["lambda"])
 sd_lambda  <- as.numeric(sd(ranef(fit_nlme_sim)[["lambda"]]))
 
 nlmeObject_sim$str_val[[9]] <- log(lambda_hat / sd_lambda)
 nlmeObject_sim$str_val$lambda <- log(lambda_hat)
 
 
 # binary data
 glmeObject2_nlme_sim <- list(
   fm= z ~ 1+month+sindoes+doesW+b_11+month*b_21,
   family='binomial', par="alpha", ran.par=c('b_11','b_21'),
   sigma=NULL, 
   disp=c("alpha4", 'alpha5'),
   lower=c(-Inf,0),
   upper=rep(Inf,2),
   str_val=c(fixef(md3_sim),  sd(ranef(md3_sim)$sid[,1])),
   CenObject=NULL)
 
 glmeObject_nlme_sim <- list(nlmeObject_sim, glmeObject2_nlme_sim)
 
 
 survObject_nlme_sim <- list(
   fm= event_time ~ base+ b_11 + b_21,
   event="event", 
   par='xi',
   disp=NULL,
   lower=NULL, upper=NULL,
   distribution=NULL,
   str_val= summary(fitCOX1_sim)$coeff[,1])
 
 
 Sim_result <- JMfit(
   glmeObject = glmeObject_nlme_sim,
   survObject = survObject_nlme_sim,
   long.data = Sim_1$dat, 
   surv.data = Sim_1$surv.dat,
   idVar="sid", 
   eventTime="event_time",
   survFit= fitCOX1_sim,
   method = "h-likelihood", Silent=T)
 
 out <- JMsummary(Sim_result)
 
 ## -------- 2-step results --------
 step1_nlme_est <- fixef(fit_nlme_sim)
 step1_nlme_sd  <- sqrt(diag(vcov(fit_nlme_sim)))
 
 step2_glmer_est <- fixef(md3_sim)
 step2_glmer_sd  <- sqrt(diag(vcov(md3_sim)))
 
 step3_cox_est <- coef(fitCOX1_sim)
 step3_cox_sd  <- sqrt(diag(vcov(fitCOX1_sim)))
 
 saveRDS(
   list(
     seed = seed,
     sim_dat = Sim_1$dat,
     sim_surv_dat = Sim_1$surv.dat,
     mydat_nlme_sim = mydat_nlme_sim,
     Sdata_nlme_sim = Sdata_nlme_sim
   ),
   file = file.path(save_dir, paste0("sim_data_seed_", seed, ".rds"))
 )
 
 list(
   success = TRUE,
   seed = seed,
   data_file = file.path(save_dir, paste0("sim_data_seed_", seed, ".rds")),
   joint = list(
     estimate = out[, "Estimate"],
     se = out[, "Std.Error"],
     param_names = rownames(out),
     summary_result = out,
     fit = Sim_result
   ),
   twostep = list(
     nlme = list(
       estimate = step1_nlme_est,
       se = step1_nlme_sd
     ),
     glmer = list(
       estimate = step2_glmer_est,
       se = step2_glmer_sd
     ),
     cox = list(
       estimate = step3_cox_est,
       se = step3_cox_sd
     )
   )
 )
}

true_par <- c(
  A1 = 1.49,
  A2 = 0.74,
  A3 = 0.18,
  A4 = 0.18,
  A5 = 0.21,
  A6 = 0.06, # 0.12
  A7 = 0.03, # 0.1
  lambda = -5.51,
  alpha0 = -1.65,
  alpha1 = 0.11,
  alpha2 = 1.6,
  alpha3 = -0.02,
  xi0 = -0.26,
  xi1 = 0.76,
  xi2 = 0.7
)


 
# nsim <- 100
# res_list <- vector("list", nsim)


# for (i in 1:nsim) {
#   message("Simulation ", i, " start: ", Sys.time())
#   tic()
#   res_list[[i]] <- tryCatch(
#     run_one_sim(seed = 1000 + i),
#     error = function(e) {
#       message("Simulation ", i, " failed: ", e$message)
#       NULL
#     }
#   )
#   toc()
# }


# ------------------ 100 simulation for NLME case -------------------------

target_success <- 100
success_count <- 0
attempt <- 0

# res_new1 <- list()  # coverage lambda -> 40  sigme = 1
# res_new2 <- list()  # sigme = 0
res_new_for_test <- list() 
while (success_count < target_success) {
  
  attempt <- attempt + 1
  message("Attempt ", attempt, " start: ", Sys.time())
  
  tic()
  
  res <- tryCatch(
    run_one_sim(seed = 2000 + attempt),
    error = function(e) {
      message("Attempt ", attempt, " failed: ", e$message)
      NULL
    }
  )
  
  toc()
  
  if (!is.null(res)) {
    
    success_count <- success_count + 1
    res_new_for_test[[success_count]] <- res
    
    message("Success ", success_count, "/", target_success)
    
  }
}


# -------------------------------------------------------------------
#

# jmsummary_sd_naive <- vector("list", length(res_new_for_test))
jmsummary_sd_robust <- vector("list", length(res_new_for_test))

for (i in seq_along(res_new_for_test)) {
  cat("Running", i, "\n")
  
  fit_i <- res_new_for_test[[i]]$joint$fit
  
  jmsummary_sd_robust[[i]] <- tryCatch({
    sd_i <- JMsd_aGH(
      fit_i,
      ghsize = 4,
      srcpath = srcpath,
      paralle = TRUE
    )$sd_robust
    
    JMsummary(fit_i, newSD = sd_i)
    
  }, error = function(e) {
    message("Error at i = ", i, ": ", e$message)
    NULL
  })
}

# ----------------------- calculate SSE , BIAS, CONVERAGE, .... --------------
est_mat <- do.call(rbind, lapply(res_new_for_test, function(x) x$joint$estimate))
se_mat  <- do.call(rbind, lapply(res_new_for_test, function(x) x$joint$se))

# est_mat <- do.call(rbind, lapply(jmsummary_sd_robust, function(x) as.numeric(x[, "Estimate"])))
# se_mat  <- do.call(rbind, lapply(jmsummary_sd_robust, function(x) as.numeric(x[, "Std.Error"])))



param_names <- res_new_for_test[[1]]$joint$param_names
colnames(est_mat) <- param_names
colnames(se_mat)  <- param_names



keep <- complete.cases(est_mat) & complete.cases(se_mat)

est_mat <- est_mat[keep, , drop = FALSE]
se_mat  <- se_mat[keep, , drop = FALSE]





# n_total <- 100
# n_success <- sum(!sapply(res_list, is.null))
# success_rate <- n_success / n_total * 100
 
mean_Se <- colMeans(se_mat)
est_mean <- colMeans(est_mat)
SSE <- apply(est_mat, 2, sd)
rBias <- (est_mean - true_par) / true_par * 100

rMSE <- sqrt(colMeans((est_mat - matrix(true_par,
                                        nrow(est_mat),
                                        length(true_par),
                                        byrow=TRUE))^2)) /
  abs(true_par) * 100


lower <- est_mat - 1.96 * se_mat
upper <- est_mat + 1.96 * se_mat

true_mat <- matrix(true_par,
                   nrow(est_mat),
                   length(true_par),
                   byrow=TRUE)

coverage <- colMeans(
  (lower <= true_mat) & (true_mat <= upper)
) * 100
result_table <- data.frame(
  Par = param_names,
  True = true_par,
  Estimate = est_mean,
  SSE = SSE,
  meanse = mean_Se,
  rBias = rBias,
  rMSE = rMSE,
  Coverage = coverage
)

start <- which(names(result_table) == "Estimate")
end   <- which(names(result_table) == "Coverage")

result_table[, start:end] <- round(result_table[, start:end], 3)
# -------------------------------------------------------------------

# empirical MC SD




all_est <- lapply(1:length(res_new_for_test), function(b){
  tmp <- res_new_for_test[[b]]$joint$summary_result
  data.frame(
    sim = b,
    param = rownames(tmp),
    estimate = tmp[, "Estimate"]
  )
}) |> dplyr::bind_rows()



mc_summary <- all_est %>%
  group_by(param) %>%
  summarise(
    Estimate = mean(estimate, na.rm = TRUE),
    Empirical_SD = sd(estimate, na.rm = TRUE),
    Bias = Estimate - true_par[first(param)],
    rBias = 100 * Bias / abs(true_par[first(param)]),
    MSE = mean((estimate - true_par[first(param)])^2, na.rm = TRUE),
    SSE = sum((estimate - true_par[first(param)])^2, na.rm = TRUE),,
    rMSE = 100 * sqrt(MSE) / abs(true_par[first(param)])
  )

mc_sd_table <- all_est %>%
  group_by(param) %>%
  summarise(
    mc_sd = sd(estimate, na.rm = TRUE),
    .groups = "drop"
  )

all_est2 <- all_est %>%
  left_join(mc_sd_table, by = "param")



coverage_empirical <- all_est2 %>%
  group_by(param) %>%
  summarise(
    True = true_par[first(param)],
    Coverage_empirical = 100 * mean(
      estimate - 1.96 * mc_sd <= True &
        True <= estimate + 1.96 * mc_sd,
      na.rm = TRUE
    ),
    .groups = "drop"
  )







# ------------------------ two-step ---------------------------------
nlme_est <- do.call(rbind,
                    lapply(res_new1, function(x) x$twostep$nlme$estimate)
)

nlme_se <- do.call(rbind,
                   lapply(res_new1, function(x) x$twostep$nlme$se)
)

true_nlme <- c(
  A1 = 1.49,
  A2 = 0.74,
  A3 = 0.18,
  A4 = 0.18,
  A5 = 0.21,
  A6 = 0.06,
  A7 = 0.03,
  lambda = exp(-5.51)
)

lower <- nlme_est - 1.96 * nlme_se
upper <- nlme_est + 1.96 * nlme_se

coverage_nlme <- colMeans(
  sweep(lower, 2, true_nlme, `<=`) &
    sweep(upper, 2, true_nlme, `>=`),
  na.rm = TRUE
)
# -------------------------------------------------------------------







# 
# 
# # ----------------------- interval censored-------------------------
# 
# run_inverval_sim <- function(seed = NULL) {
#   if (!is.null(seed)) set.seed(seed)
#   
# 
#     
#   Sim_2 <-  SimData3(n=200, 
#                      ZX=c('month','sindoes','doesW', 'trb1'), 
#                      YX=c("month",'month2','sindoes'), sigma_lme=1.3, 
#                      A=NULL, 
#                      lambda = NULL, sigma_nlme = NULL,
#                      Vtime=c(0,1,6,12,18,24,30),
#                      endTrial=36,
#                      alpha=c(-1.2, 0.08, 0.9, -0.02, 0.15, 0.03), 
#                      beta = c(0.3, 1.8, 1.8, 0.97, -0.9), delim_val=1.47, regime=1, 
#                      rate=NULL, 
#                      Xi=c(0, 0.6), Asso = c(0.4, 0.05), #(For surivival model) Xi -> parameter for fixed effect; Asso -> association parameter in random effect 
#                      scale=50, shape=5, 
#                      censor_type = "right",
#                      ran_p=c(1,1), # number of random effect for each model
#                      SIGMA=diag(c(1, 1)) , 
#                      cen_par = c(5, 25),
#                      check=F,
#                      freq=T,  # high frequency of check point
#                      contaminated=F, # if there is outlier in lme
#                      percent=0 #percent of outlier
#   ) 
#   
#   
#   
#   # (1) Model 1: LME of y
#   fm1_sim <-  y ~ 1+ month+month2+sindoes+(1|sid)
#   md1_sim <- lmer(fm1_sim, data=Sim_2$dat)
#   
#   # get the estimated random effect by observation
#   estBi <- data.frame(row.names(ranef(md1_sim)$sid), 
#                       scale(ranef(md1_sim)$sid, center=T,scale=T))
#   names(estBi) <- c("sid", "estb11")
#   mydat_sim <- merge(Sim_2$dat, estBi, by = "sid", all.x = TRUE)
#   
#   
#   # (1.5) Model 1.5: a GLME model of C
#   # using the estimated random intercept from model 1 (i.e. estb11, scaled) as a covariate
#   fm2 <- c ~ 1+month+month2+sindoes+ estb11
#   md2 <- glm(fm2, family=binomial, data=mydat_sim)
#   
#   # (2) Model 2: a GLME model of C
#   ## using the estimated random intercept from model 1 (i.e. estb11, scaled) as a covariate
#   fm2_sim <- z ~ 1+month+sindoes+doesW+estb11+(month-1|sid)
#   md2_sim <- glmer(fm2_sim, family="binomial", data=mydat_sim)
#   
#   
#   
#   # (3) Model 3: a survival model
#   Sdata_sim <- Sim_2$surv.dat
#   Sdata_sim$estb11 <- scale(ranef(md1_sim)$sid[,1], center=T, scale=T)
#   Sdata_sim$estb21 <- scale(ranef(md2_sim)$sid[,1], center=T, scale=T)
#   
#   fitCOX2_sim <- coxph(Surv(event_time, event) ~ base+ estb11 + estb21, data =  Sdata_sim)
#   # a interval weibull model
#   # fitcox2_sim <- survreg(
#   #   Surv(time = L, time2 = R, type = "interval2") ~ base+ estb11+estb21,
#   #   data = Sdata_sim,
#   #   dist = "weibull"
#   # )
#   ###############################
#   ######   Joint modeling   #####
#   ###############################
#   LLOQ = 2 # lower limit of quantification of Y
#   
#   CenObject <- list(
#     fm= c ~ 1  + month +month2 +sindoes + b11,
#     family='binomial', par='eta', ran.par="b11",
#     disp="eta4",
#     lower= -Inf, upper=Inf,
#     str_val=coef(md2),
#     Cregime=1,
#     truncated=T, delim_val=LLOQ)
# 
#   
#   glmeObject1_sim <- list(
#     fm= y ~ 1 + month + month2 + sindoes + b11,
#     family='normal', par="beta", ran.par='b11', sigma='sigma',
#     disp='beta4',
#     lower=0, upper=Inf,   
#     str_val=c(fixef(md1_sim), sd(ranef(md1_sim)$sid[,1])),
#     CenObject=NULL)
#   
#   glmeObject2_sim <- list(
#     fm= z ~ 1 + month + sindoes +doesW + b11 + month*b21,
#     family='binomial', par="alpha", ran.par=c('b11','b21'),
#     sigma=NULL, 
#     disp=c("alpha4", 'alpha5'),
#     lower=c(-Inf,0),
#     upper=rep(Inf,2),
#     str_val=c(fixef(md2_sim), sd(ranef(md2_sim)$sid[,1])),
#     CenObject=NULL)
#   
#   glmeObject_sim <- list(glmeObject1_sim, glmeObject2_sim)
#   
#   
#   ## Weibull interval cencored
#   # survObject_intPH <- list( fm = Surv(L, R) ~ base+ b11+b21, 
#   #                             event = "event", 
#   #                           par = "lambda", # betas: lambda0, lambda1, ... 
#   #                           disp = NULL, 
#   #                           lower=c(0, 0), 
#   #                           upper=c(Inf, Inf), 
#   #                           distribution = "weibull_ph_interval", 
#   #                           str_val = -summary(fitcox2_sim)$coeff[-1]/summary(fitcox2_sim)$scale) 
#   # 
#   
#   survObject_sim <- list(
#     fm= event_time ~ base+ b11 + b21,
#     event="event", 
#     par='xi',
#     disp=NULL,
#     lower=NULL, upper=NULL,
#     distribution=NULL,
#     str_val= summary(fitCOX2_sim)$coeff[,1])
#   
#   
#   Sim_result <- JMfit(
#     glmeObject = glmeObject_sim,
#     survObject = survObject_sim,
#     long.data = Sim_2$dat, 
#     surv.data = Sim_2$surv.dat,
#     idVar="sid", 
#     eventTime="event_time",
#     survFit= fitCOX2_sim,
#     method = "h-likelihood", Silent=T)
#   
#   
#   
#   
#   
#   
#   
#   out <- JMsummary(Sim_result)
#   
#   ## -------- 2-step results --------
#   step1_lme_est <- fixef(md1_sim)
#   step1_lme_sd  <- sqrt(diag(vcov(md1_sim)))
#   
#   step2_glmer_est <- fixef(md2_sim)
#   step2_glmer_sd  <- sqrt(diag(vcov(md2_sim)))
#   
#   step3_weibull_est <- coef(fitcox2_sim)
#   step3_weibull_sd  <- sqrt(diag(vcov(fitcox2_sim)))
#   
#   list(
#     seed = seed,
#     joint = list(
#       estimate = out[, "Estimate"],
#       se = out[, "Std.Error"],
#       param_names = rownames(out),
#       raw = out
#     ),
#     twostep = list(
#       nlme = list(
#         estimate = step1_lme_est,
#         se = step1_lme_sd
#       ),
#       glmer = list(
#         estimate = step2_glmer_est,
#         se = step2_glmer_sd
#       ),
#       weibull = list(
#         estimate = step3_weibull_est,
#         se = step3_weibull_sd
#       )
#     )
#   )
# }
# 
# 
# 
# 
# 
# 
# 
# 
# # -------- run function -------------------
# target_success_2 <- 1
# success_count_2 <- 0
# attempt_2 <- 0
# 
# res_new00 <- list()
# 
# while (success_count_2 < target_success_2) {
#   
#   attempt_2 <- attempt_2 + 1
#   message("Attempt ", attempt_2, " start: ", Sys.time())
#   
#   tic()
#   
#   res <- tryCatch(
#     run_inverval_sim(seed = 1000 + attempt_2),
#     error = function(e) {
#       message("Attempt ", attempt_2, " failed: ", e$message)
#       NULL
#     }
#   )
#   
#   toc()
#   
#   if (!is.null(res)) {
#     
#     success_count <- success_count_2 + 1
#     res_new00[[success_count_2]] <- res
#     
#     message("Success ", success_count_2, "/", target_success_2)
#     
#   }
# }
# 
# 

lambda_vec <- sapply(res_new_for_test, function(x) x$joint$estimate[8])
lambda_se <- sapply(res_new_for_test, function(x) x$joint$fit$fixedsd[8])

exp(lambda_vec)*lambda_se
qqnorm(exp(lambda_vec)*lambda_se)
qqline(lambda_vec)

lower <- exp(lambda_vec - 1.96 * lambda_se)
upper <- exp(lambda_vec + 1.96 * lambda_se)

covered <- (exp(-5.51) >= lower) & (exp(-5.51) <= upper)
coverage <- mean(covered, na.rm = TRUE)
coverage <- mean(covered, na.rm = TRUE)
