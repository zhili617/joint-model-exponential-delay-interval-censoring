library(tidyverse)
library(nlme)
library(DescTools)

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
  dataset_new[[lag_col]] <- with(dataset_new, pmax(dosedays + (mdays[pmax(m,1)] - mdays[k]), 0))
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
  random =  list(sid = pdDiag(A1 + A2 + A3 + A4 + A5 + A6 + A7 + lambda~ 1)),
  start = as.numeric(coef(fit_nls)),
  control = nlmeControl(
    maxIter=200, pnlsMaxIter=300, msMaxIter=200,
    tolerance = 1e-5, pnlsTol = 1e-5)
)

summary(fit_nlme)

# equation 2
nls_formula <- NAb ~ 
  A1*ind_k1*(-expm1(-beta*lag_k1))*exp(-lambda1*lag_k1) +
  A2*ind_k2*(-expm1(-beta*lag_k2))*exp(-lambda2*lag_k2) +
  A3*ind_k3*(-expm1(-beta*lag_k3))*exp(-lambda3*lag_k3) +
  A4*ind_k4*(-expm1(-beta*lag_k4))*exp(-lambda4*lag_k4) +
  A5*ind_k5*(-expm1(-beta*lag_k5))*exp(-lambda5*lag_k5) +
  A6*ind_k6*(-expm1(-beta*lag_k6))*exp(-lambda6*lag_k6) +
  A7*ind_k7*(-expm1(-beta*lag_k7))*exp(-lambda7*lag_k7)

start_nls <- list(
  A1=.5,A2=.6,A3=.6,A4=.6,A5=.5,A6=.5,A7=.5,
  lambda=.30,  
  beta=.80
)

fit_nls <- nls(
  nls_formula, data=dataset_new, start=start_nls, algorithm="port",
  lower=c(A1=0,A2=0,A3=0,A4=0,A5=0,A6=0,A7=0, lambda=0.05, beta=0.05),
  upper=c(A1=Inf,A2=Inf,A3=Inf,A4=Inf,A5=Inf,A6=Inf,A7=Inf, lambda=2,    beta=2),
  control=nls.control(maxiter=500, tol=1e-4, minFactor=1/2048)
)


nlme_formula <- nls_formula




fit_nlme <- nlme(
  model = nlme_formula,
  data  = groupedData(NAb ~ sampledays | sid, data = dataset_new),
  fixed = A1 + A2 + A3 + A4 + A5 + A6 + A7 + beta + lambda ~ 1,
  random = list(sid = pdDiag(A1 + A2 + A3 + A4 + A5 + A6 + A7 + beta + lambda ~ 1)),  
  start  = as.numeric(coef(fit_nls)),
  control = nlmeControl(
    pnlsMaxIter = 100, msMaxIter = 500, maxIter = 300,
    pnlsTol     = 1e-4,  
    tolerance   = 1e-5,  
    msVerbose   = TRUE, returnObject = TRUE
  )
)









# model 2 

for (k in 1:7) {
  L <- paste0("lag_k", k)
  dataset_new[[L]] <- pmax(dataset_new[[L]], 0) + eps
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
    maxIter     = 200, pnlsMaxIter = 30, msMaxIter = 200,
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

#----------------------value compare (md1 & method 1)------------------------------------------

df <- dataset_new %>%
  mutate(
    pred_nlme = tryCatch(predict(fit_nlme, dataset_new, level = 1), error = function(e) NA),
    pred_nls  = tryCatch(predict(fit_nls,  dataset_new),             error = function(e) NA),
    pred_md1 = tryCatch(predict(md1), error = function(e) NA)
  )



rmse_nlme  <- sqrt(mean((df$pred_nlme - df$NAb)^2, na.rm = TRUE))
mae_nlme   <- mean(abs(df$pred_nlme - df$NAb),      na.rm = TRUE)
smape_nlme <- mean(2*abs(df$pred_nlme - df$NAb) / (abs(df$pred_nlme) + abs(df$NAb) + 1e-8), na.rm = TRUE) * 100
bias_nlme  <- mean(df$pred_nlme - df$NAb, na.rm = TRUE)
r2_nlme    <- cor(df$pred_nlme, df$NAb, use = "complete.obs")^2
ccc_nlme   <- CCC(df$pred_nlme, df$NAb)$rho.c[1]  


rmse_nls  <- sqrt(mean((df$pred_nls - df$NAb)^2, na.rm = TRUE))
mae_nls   <- mean(abs(df$pred_nls - df$NAb),      na.rm = TRUE)
smape_nls <- mean(2*abs(df$pred_nls - df$NAb) / (abs(df$pred_nls) + abs(df$NAb) + 1e-8), na.rm = TRUE) * 100
bias_nls  <- mean(df$pred_nls - df$NAb, na.rm = TRUE)
r2_nls    <- cor(df$pred_nls, df$NAb, use = "complete.obs")^2
ccc_nls   <- CCC(df$pred_nls, df$NAb)$rho.c[1]  


rmse_md1  <- sqrt(mean((df$pred_md1 - df$NAb)^2, na.rm = TRUE))
mae_md1   <- mean(abs(df$pred_md1 - df$NAb),      na.rm = TRUE)
smape_md1 <- mean(2*abs(df$pred_md1 - df$NAb) / (abs(df$pred_md1) + abs(df$NAb) + 1e-8), na.rm = TRUE) * 100
bias_md1  <- mean(df$pred_md1 - df$NAb, na.rm = TRUE)
r2_md1    <- cor(df$pred_md1, df$NAb, use = "complete.obs")^2
ccc_md1   <- CCC(df$pred_md1, df$NAb)$rho.c[1]  




cal_int  <- c(
  coef(lm(NAb ~ pred_nlme, data = df))[1],
  coef(lm(NAb ~ pred_nls,  data = df))[1],
  coef(lm(NAb ~ pred_md1,      data = df))[1]
)
cal_slope <- c(
  coef(lm(NAb ~ pred_nlme, data = df))[2],
  coef(lm(NAb ~ pred_nls,  data = df))[2],
  coef(lm(NAb ~ pred_md1,      data = df))[2]
)

results <- tibble(
  model      = c("nlme", "nls", "md1"),
  RMSE       = c(rmse_nlme,  rmse_nls,  rmse_md1),
  MAE        = c(mae_nlme,   mae_nls,   mae_md1),
  SMAPE_pct  = c(smape_nlme, smape_nls, smape_md1),
  Bias       = c(bias_nlme,  bias_nls,  bias_md1),
  R2         = c(r2_nlme,    r2_nls,    r2_md1),
  CCC        = as.numeric(c(ccc_nlme,   ccc_nls,   ccc_md1)),
  Cal_intercept = cal_int,
  Cal_slope     = cal_slope
) %>%
  mutate(across(-model, ~ round(., 4))) %>%
  arrange(RMSE)

results







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

par(mfrow=c(1,3))  # put 3 graph into same line
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
f_model2 <- function(t, C1,C2,C3,C4,C5,C6,C7,m, plan_months=c(0,1,6,12,18,24,30)) {
  mdays <- plan_months * 30.44
  
  lag1 <- pmax(t - mdays[1], 1e-8)
  lag2 <- pmax(t - mdays[2], 1e-8)
  lag3 <- pmax(t - mdays[3], 1e-8)
  lag4 <- pmax(t - mdays[4], 1e-8)
  lag5 <- pmax(t - mdays[5], 1e-8)
  lag6 <- pmax(t - mdays[6], 1e-8)
  lag7 <- pmax(t - mdays[7], 1e-8)
  
  val <- C1*(t>=mdays[1])*(lag1^(-m)) +
    C2*(t>=mdays[2])*(lag2^(-m)) +
    C3*(t>=mdays[3])*(lag3^(-m)) +
    C4*(t>=mdays[4])*(lag4^(-m)) +
    C5*(t>=mdays[5])*(lag5^(-m)) +
    C6*(t>=mdays[6])*(lag6^(-m)) +
    C7*(t>=mdays[7])*(lag7^(-m))
  
  return(val)
}




set.seed(10)
sids_sample <- sample(unique(dataset_new$sid), 3)
tgrid <- seq(min(dataset_new$sampledays),
             max(dataset_new$sampledays), length.out=200)

par(mfrow=c(1,3))
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


