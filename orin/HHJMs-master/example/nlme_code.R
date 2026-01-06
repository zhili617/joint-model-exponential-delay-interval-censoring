library(tidyverse)
library(nlme)
library(DescTools)
library(survival)

# dataset <- read_csv("E:/UBC/Final Project/vax004dataAb.txt")
dataset <- read_csv("E:/UBC/Final Project/vax004dataAb.txt")



# For simplification, if an sid has NA value in any variable, then delate all the rows of this sid
dataset_new <- dataset %>%
  group_by(sid) %>%
  filter(Vxind == 1, HIVinfectionind == 0) %>%      # 2815
  filter(n() > 1) %>%                             # Keep only sids with more than 1 record     2759
  dplyr::select(sid, sampledays, dosedays=doesdays, 
                dosenumberschedule=doesnumberschedule, 
                NAb, MNGNE8, riskscorecat, GNE8_CD4, MN_CD4) %>%
  filter(!any(is.na(across(everything())))) %>%  # Remove rows with NA
  ungroup()     #2235

dataset_new$sid <- as.factor(dataset_new$sid)

#----------------------------------------vaccine time to vaccine month-------------------
month_map <- c(`1` = 0, `2` = 1, `3` = 6, `4` = 12, `5` = 18, `6`= 24, `7` = 30)
dataset_new$vaccine_month <- month_map[as.character(dataset_new$dosenumberschedule)]


# method 1 : dosedays + (dose month[close] - dose month[k])
dataset_new <- dataset_new %>%
  mutate(
    dosedays = dplyr::coalesce(dosedays, dosedays),
    # Scheduled order of most recent immunization (1..7)
    m = as.integer(dosenumberschedule)  
  )

plan_months <- c(0, 1, 6, 12, 18, 24, 30)
mdays <- plan_months * 30.44

mk_cols <- paste0("k", 1:7)
for (k in 1:7) {
  lag_col <- paste0("lag_", mk_cols[k])
  ind_col <- paste0("ind_", mk_cols[k])
  dataset_new[[lag_col]] <- with(dataset_new, pmax((mdays[pmax(m,1)] - mdays[k]), 0))
  dataset_new[[ind_col]] <- with(dataset_new, ifelse(!is.na(m) & m >= k, 1L, 0L))
}



# method 2: sampledays - dose month
plan_months <- c(0, 1, 6, 12, 18, 24, 30)
mdays <- plan_months * 30.44   # dose month into dose day

for (k in 1:7) {
  lag_col <- paste0("lag_k", k)   
  ind_col <- paste0("ind_k", k)
  
  # ensure result >= 0
  dataset_new[[lag_col]] <- pmax(dataset_new$sampledays - mdays[k], 0)
  
  # make sure we choose sampledays > dose month
  dataset_new[[ind_col]] <- as.integer(dataset_new$sampledays + 1e-8 >= mdays[k])
}

# non-linear equation Sum_{k=1..7} A_k * I_{k <= m} * exp(-lambda * lag_k)
nl_formula <- NAb ~ 
  A1 * ind_k1 * exp(-lambda * lag_k1) +
  A2 * ind_k2 * exp(-lambda * lag_k2) +
  A3 * ind_k3 * exp(-lambda * lag_k3) +
  A4 * ind_k4 * exp(-lambda * lag_k4) +
  A5 * ind_k5 * exp(-lambda * lag_k5) +
  A6 * ind_k6 * exp(-lambda * lag_k6) +
  A7 * ind_k7 * exp(-lambda * lag_k7)

start_nls <- list(
  A1=.5,A2=.6,A3=.6,A4=.6,A5=.5,A6=.5,A7=.5,
  lambda=.30
)

fit_nls <- nls(
  nl_formula, data=dataset_new, start=start_nls, algorithm="port",
  lower=c(A1=0,A2=0,A3=0,A4=0,A5=0,A6=0,A7=0, lambda=0.05),
  upper=c(A1=Inf,A2=Inf,A3=Inf,A4=Inf,A5=Inf,A6=Inf,A7=Inf, lambda=2),
  control=nls.control(maxiter=500, tol=1e-4, minFactor=1/2048)
)



fit_nlme <- nlme(
  model = nl_formula,
  data  = groupedData(NAb ~ sampledays | sid, data = dataset_new),
  fixed = A1 + A2 + A3 + A4 + A5 + A6 + A7 + lambda ~ 1,
  random =  lambda ~ 1 | sid,
  start = as.numeric(coef(fit_nls)),
  control = nlmeControl(
    maxIter=200, pnlsMaxIter=300, msMaxIter=200,
    tolerance = 1e-5, pnlsTol = 1e-5)
)

summary(fit_nlme)







# Power Law

for (k in 1:7) {
  lag_col <- paste0("lag_k", k)
  ind_col <- paste0("ind_k", k)
  
  raw_lag <- dataset_new$sampledays - mdays[k]   
  
  dataset_new[[ind_col]] <- as.integer(raw_lag > 0)       
  dataset_new[[lag_col]] <- ifelse(raw_lag > 0, raw_lag, 1)  
}

# Sum_{k=1..7} C_k * I_{k<=m} * (lag_k)^(-m)
nl_formula_pl <- NAb ~
  C1 * ind_k1 * (lag_k1 )^(-m) +
  C2 * ind_k2 * (lag_k2 )^(-m) +
  C3 * ind_k3 * (lag_k3 )^(-m) +
  C4 * ind_k4 * (lag_k4 )^(-m) +
  C5 * ind_k5 * (lag_k5 )^(-m) +
  C6 * ind_k6 * (lag_k6 )^(-m) +
  C7 * ind_k7 * (lag_k7 )^(-m)

start_vals_pl <- c(C1=.3, C2=.3, C3=.3, C4=.3, C5=.3, C6=.3, C7=.3, m=0.6)

fit_nlme_pl <- nlme(
  model  = nl_formula_pl,
  data   = groupedData(NAb ~ sampledays | sid, data = dataset_new),
  fixed  = C1 + C2 + C3 + C4 + C5 + C6 + C7 + m ~ 1,
  random = m ~ 1 | sid,
  start  = start_vals_pl,
  control = nlmeControl(
    maxIter     = 200, pnlsMaxIter = 300, msMaxIter = 200,
    tolerance   = 1e-6, pnlsTol    = 1e-6, msTol   = 1e-6,
    returnObject = TRUE
  )
)
summary(fit_nlme)
summary(fit_nlme_pl)
AIC(fit_nlme_pl); BIC(fit_nlme_pl)
intervals(fit_nlme_pl)
AIC(fit_nlme)













set.seed(2) 
selected_sids <- sample(unique(dataset_new$sid), 4)

dataset_plot <- dataset_new %>%
  mutate(
    observed = NAb,
    fit = predict(fit_nlme, level = 0)  
  )

plot_data <- dataset_plot %>%
  filter(sid %in% selected_sids) %>%
  pivot_longer(cols = c(observed, fit),
               names_to = "type",
               values_to = "value")

ggplot(plot_data, aes(x = sampledays, y = value, color = type)) +
  geom_point(data = subset(plot_data, type == "observed")) +
  geom_line(data = subset(plot_data, type == "fit")) +
  facet_wrap(~ sid, scales = "free_x") +
  labs(
    title = paste("Randomly selected sids:", paste(selected_sids, collapse = ", ")),
    x = "Days", 
    y = "NAb value", 
    color = "Legend"
  ) +
  theme_minimal()






# ------------predicted graph in total --------------------------
dataset_plot <- dataset_new %>%
  mutate(
    observed = NAb,
    fit = predict(fit_nlme, level = 1)   
  )

plot_data <- dataset_plot %>%
  pivot_longer(cols = c(observed, fit),
               names_to = "type",
               values_to = "value")

ggplot(plot_data, aes(x = sampledays, y = value, color = type, group = interaction(sid, type))) +
  geom_point(data = subset(plot_data, type == "observed"), alpha = 0.5) +
  geom_line(data = subset(plot_data, type == "fit"), alpha = 0.7) +
  labs(
    title = "Observed (points) and Predicted (lines) for all sids",
    x = "Days", 
    y = "NAb value", 
    color = "Legend"
  ) +
  theme_minimal()





# ----------------Original NAb graph----------------------------
ggplot(dataset_new, aes(x = sampledays, y = NAb, group = sid)) +
  geom_line(alpha = 0.6, color = "steelblue") +
  geom_point(alpha = 0.4, color = "steelblue") +
  labs(
    title = "Original observed NAb trajectories (all sids)",
    x = "Days",
    y = "NAb value"
  ) +
  theme_minimal()


# random 5 cases
set.seed(222)  
sid_sample <- sample(unique(dataset_new$sid), 5)

random_NAb_ori <- dataset_new %>%
  dplyr::filter(sid %in% sid_sample)

ggplot(random_NAb_ori, aes(x = sampledays, y = NAb, color = sid, group = sid)) +
  geom_line(alpha = 0.8, linewidth = 1) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = "Randomly Selected 5 NAb Trajectories",
    x = "Days",
    y = "NAb value"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())

#----------------goodness of fit(compare fitted value v.s. model curve)


# ----------------------------------------
# Model 1
# ----------------------------------------

set.seed(10)
sids_sample <- sample(unique(dataset_new$sid), 3) #randomly choose 3 sid

tgrid <- seq(min(dataset_new$sampledays),
             max(dataset_new$sampledays), length.out=200)

# define model formula
f_model1 <- function(t, A, lambda, plan_months=c(0,1,6,12,18,24,30)) {
  mdays <- plan_months * 30.44
  val <- rep(0, length(t))
  for (k in 1:7) {
    Ak <- A[k]
    val <- val + Ak * as.integer(t >= mdays[k]) *
      exp(-lambda * pmax(t - mdays[k], 0))
  }
  return(val)
}

par(mfrow=c(3,1))  # put 3 graph into same line
for (sid_i in sids_sample) {
  
  dat_i <- dataset_new %>% filter(sid == sid_i)
  
  # orginal point (from data)
  plot(dat_i$sampledays, dat_i$NAb, pch=16, col="black",
       main=paste("SID =", sid_i), xlab="Days", ylab="NAb")
  
  # predicted value 
  pred_points <- predict(fit_nlme, dat_i, level=1)
  points(dat_i$sampledays, pred_points, col="red", pch=19)
  
  # model curve (based on fixed effects parameter)
  coefs <- fixef(fit_nlme)
  A <- coefs[paste0("A", 1:7)]
  lambda <- coefs["lambda"]
  curve_model <- f_model1(tgrid, A, lambda)
  lines(tgrid, curve_model, col="blue", lwd=2)
  
  legend("topright", legend=c("Observed","nlme fit","Formula curve"),
         col=c("black","red","blue"), pch=c(16,19,NA),
         lty=c(NA,NA,1), bty="n")
}




#-----------------------------------------------
# Model 3
#-----------------------------------------------
# define power-law function
f_model2 <- function(t, C1, C2, C3, C4, C5, C6, C7, m,
                     plan_months = c(0,1,6,12,18,24,30)) {

  mdays <- plan_months * 30.44

  # pure lag (no repair)
  lag1 <- t - mdays[1]
  lag2 <- t - mdays[2]
  lag3 <- t - mdays[3]
  lag4 <- t - mdays[4]
  lag5 <- t - mdays[5]
  lag6 <- t - mdays[6]
  lag7 <- t - mdays[7]

  # effect only when lag > 0
  term1 <- ifelse(lag1 > 0, C1 * lag1^(-m), 0)
  term2 <- ifelse(lag2 > 0, C2 * lag2^(-m), 0)
  term3 <- ifelse(lag3 > 0, C3 * lag3^(-m), 0)
  term4 <- ifelse(lag4 > 0, C4 * lag4^(-m), 0)
  term5 <- ifelse(lag5 > 0, C5 * lag5^(-m), 0)
  term6 <- ifelse(lag6 > 0, C6 * lag6^(-m), 0)
  term7 <- ifelse(lag7 > 0, C7 * lag7^(-m), 0)

  return(term1 + term2 + term3 + term4 + term5 + term6 + term7)
}





tgrid <- seq(min(dataset_new$sampledays),
             max(dataset_new$sampledays), length.out=200)

par(mfrow=c(3,1))
for (sid_i in sids_sample) {
  
  dat_i <- dataset_new %>% filter(sid == sid_i)
  
  plot(dat_i$sampledays, dat_i$NAb, pch=16, col="black",
       main=paste("SID =", sid_i), xlab="Days", ylab="NAb")
  
  pred_points <- predict(fit_nlme_pl, dat_i, level=1)
  points(dat_i$sampledays, pred_points, col="red", pch=19)
  
  # find fixed effects parameter
  coefs_pl <- fixef(fit_nlme_pl)
  C <- coefs_pl[paste0("C", 1:7)]
  m <- coefs_pl["m"]
  
  # model curve
  curve_model2 <- f_model2(
    tgrid,
    C1 = coefs_pl["C1"], C2 = coefs_pl["C2"], C3 = coefs_pl["C3"],
    C4 = coefs_pl["C4"], C5 = coefs_pl["C5"], C6 = coefs_pl["C6"], C7 = coefs_pl["C7"],
    m  = coefs_pl["m"]
  )
  
  lines(tgrid, curve_model2, col="blue", lwd=2)
  
  legend("topright", legend=c("Observed","nlme fit","Formula curve"),
         col=c("black","red","blue"), pch=c(16,19,NA),
         lty=c(NA,NA,1), bty="n")
}


#------------------------ two step method -------------------

# =========================
# Improved two-step (bootstrap) for:
#   Step 1: NLME -> get subject-specific random effects (lambda_i)
#   Step 2: Cox  -> treat lambda_i as covariate
#   Repeat B times (cluster bootstrap by subject) to get robust SE/CI
# =========================

library(nlme)
library(survival)
library(dplyr)
library(tibble)

# ---------- helper: build a bootstrap sample with unique boot_sid ----------


make_boot_sample <- function(long_data, surv_data, idVar = "sid", seed = NULL) { 
  if (!is.null(seed)) set.seed(seed) 
  sids <- sort(unique(long_data[[idVar]])) 
  n_sid <- length(sids) 
  sid_star <- sample(sids, size = n_sid, replace = TRUE) # cut by sid
  long_list <- split(long_data, long_data[[idVar]]) 
  surv_list <- split(surv_data, surv_data[[idVar]]) 
  long_b_list <- vector("list", n_sid) 
  surv_b_list <- vector("list", n_sid) 
  for (k in seq_len(n_sid)) { 
    s <- as.character(sid_star[k]) 
    boot_sid <- k 
    tmp_long <- long_list[[s]] 
    tmp_surv <- surv_list[[s]] 
    tmp_long[[idVar]] <- boot_sid 
    tmp_surv[[idVar]] <- boot_sid 
    long_b_list[[k]] <- tmp_long 
    surv_b_list[[k]] <- tmp_surv } 
  long_b <- dplyr::bind_rows(long_b_list) 
  surv_b <- dplyr::bind_rows(surv_b_list) 
  list(long_b = long_b, surv_b = surv_b) }

# ---------- helper: fit nlme safely ----------
fit_nlme_safe <- function(long_b,
                          nl_formula,
                          start,
                          idVar = "sid",
                          control = nlmeControl(maxIter=200, pnlsMaxIter=300, msMaxIter=200,
                                                tolerance=1e-5, pnlsTol=1e-5,
                                                returnObject=TRUE)) {
  # groupedData needs the same id variable name
  gd <- groupedData(formula = as.formula(paste0("NAb ~ sampledays | ", idVar)), data = long_b)
  
  fit <- tryCatch(
    nlme(
      model = nl_formula,
      data  = gd,
      fixed = A1 + A2 + A3 + A4 + A5 + A6 + A7 + lambda ~ 1,
      random = lambda ~ 1 | sid,
      start = start,
      control = control,
      na.action = na.exclude
    ),
    error = function(e) e
  )
  fit
}

# ---------- helper: extract lambda_i and merge to survival ----------
add_lambda_to_surv <- function(fit_nlme_b, surv_b, idVar = "sid",
                               use = c("lambda_i","b_lambda"),
                               standardize = TRUE) {
  use <- match.arg(use)
  
  # ranef() gives BLUPs per subject
  re <- ranef(fit_nlme_b) %>%
    rownames_to_column(idVar) %>%
    mutate(!!idVar := as.integer(.data[[idVar]]))
  
  # fixed effect lambda
  lam_fix <- fixef(fit_nlme_b)["lambda"]
  
  re <- re %>%
    rename(b_lambda = lambda) %>%
    mutate(lambda_i = as.numeric(lam_fix + b_lambda))
  
  out <- surv_b %>% left_join(re %>% dplyr::select(all_of(idVar), b_lambda, lambda_i), by = idVar)
  
  if (use == "lambda_i") {
    out$lambda_cov <- out$lambda_i
  } else {
    out$lambda_cov <- out$b_lambda
  }
  
  if (standardize) {
    out$lambda_cov <- as.numeric(scale(out$lambda_cov))
  }
  
  out
}

# ---------- helper: fit Cox safely ----------
fit_cox_safe <- function(surv_df,
                         timeVar = "L",
                         eventVar = "dropped_out_right",
                         covars = c("GNE8_CD4","risk1","risk2","lambda_cov")) {
  fml <- as.formula(
    paste0("Surv(", timeVar, ", ", eventVar, ") ~ ", paste(covars, collapse = " + "))
  )
  fit <- tryCatch(coxph(fml, data = surv_df, x = TRUE), error = function(e) e)
  fit
}

# ---------- MAIN: improved two-step via bootstrap ----------
two_step_bootstrap <- function(long_data,
                               surv_data_final,
                               nl_formula,
                               fit_nlme_ref,          # your original fit_nlme (for start/control)
                               B = 300,
                               seed = 1,
                               idVar = "sid",
                               timeVar = "L",
                               eventVar = "dropped_out_right",
                               lambda_use = c("lambda_i","b_lambda"),
                               standardize_lambda = TRUE,
                               verbose_every = 50) {
  lambda_use <- match.arg(lambda_use)
  set.seed(seed)
  
  # starting values: reuse reference fit
  start <- fixef(fit_nlme_ref)  # includes A1..A7 and lambda
  ctrl  <- fit_nlme_ref$control
  
  # storage
  cox_names <- c("GNE8_CD4","risk1","risk2","lambda_cov")
  beta_mat <- matrix(NA_real_, nrow = B, ncol = length(cox_names),
                     dimnames = list(NULL, cox_names))
  ok_nlme <- logical(B)
  ok_cox  <- logical(B)
  err_msg <- character(B)
  
  for (b in 1:B) {
    if (!is.null(verbose_every) && b %% verbose_every == 0) {
      message("Bootstrap ", b, "/", B)
    }
    
    samp <- make_boot_sample(long_data, surv_data_final, idVar = idVar)
    long_b <- samp$long_b
    surv_b <- samp$surv_b
    
    # Step 1: NLME
    fit_b <- fit_nlme_safe(long_b, nl_formula = nl_formula, start = start, idVar = idVar, control = ctrl)
    if (inherits(fit_b, "error")) {
      ok_nlme[b] <- FALSE
      err_msg[b] <- paste0("NLME error: ", fit_b$message)
      next
    }
    ok_nlme[b] <- TRUE
    
    # Build survival covariate (lambda)
    surv_b2 <- add_lambda_to_surv(fit_b, surv_b, idVar = idVar,
                                  use = lambda_use, standardize = standardize_lambda)
    
    # Step 2: Cox
    cox_b <- fit_cox_safe(surv_b2, timeVar = timeVar, eventVar = eventVar, covars = cox_names)
    if (inherits(cox_b, "error")) {
      ok_cox[b] <- FALSE
      err_msg[b] <- paste0("Cox error: ", cox_b$message)
      next
    }
    ok_cox[b] <- TRUE
    
    # store betas
    bb <- coef(cox_b)
    # ensure all exist
    beta_mat[b, names(bb)] <- as.numeric(bb)
  }
  
  beta_df <- as.data.frame(beta_mat) %>%
    mutate(iter = 1:B, ok_nlme = ok_nlme, ok_cox = ok_cox, err = err_msg) %>%
    relocate(iter, ok_nlme, ok_cox, err)
  
  # summary on successful runs
  keep <- which(ok_nlme & ok_cox)
  beta_keep <- beta_mat[keep, , drop = FALSE]
  
  summ <- tibble(
    term = colnames(beta_keep),
    mean = apply(beta_keep, 2, mean, na.rm = TRUE),
    se   = apply(beta_keep, 2, sd,   na.rm = TRUE),
    q025 = apply(beta_keep, 2, quantile, probs = 0.025, na.rm = TRUE),
    q975 = apply(beta_keep, 2, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
    mutate(
      HR_mean = exp(mean),
      HR_q025 = exp(q025),
      HR_q975 = exp(q975),
      n_success = length(keep),
      B_total = B,
      lambda_used = lambda_use,
      lambda_standardized = standardize_lambda
    )
  
  list(
    beta_each = beta_df,
    summary = summ
  )
}
res <- two_step_bootstrap(
  long_data = dataset_new,
  surv_data_final = surv_data_final,
  nl_formula = nl_formula,
  fit_nlme_ref = fit_nlme,
  B = 300,
  seed = 1,
  idVar = "sid",
  timeVar = "L",
  eventVar = "dropped_out_right",
  lambda_use = "lambda_i",          # 或 "b_lambda"
  standardize_lambda = TRUE,
  verbose_every = 50
)

res$summary
head(res$beta_each)
