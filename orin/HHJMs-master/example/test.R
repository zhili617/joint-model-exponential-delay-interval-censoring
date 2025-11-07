library(tidyverse)
library(readr)
library(lme4)
library(tictoc)
library(survival)
library(MASS)
library(glmmTMB)
library(Matrix)
library(ggplot2)


## source all the R code
srcpath = "E:/UBC/Final Project/orin/orin/HHJMs-master/R"
setwd(srcpath)
(file.sources = list.files(pattern="\\.[rR]$"))
sapply(file.sources,source,.GlobalEnv)

--------------# decide interval censored as drop out---------------

dataset <- read_csv("E:/UBC/Final Project/vax004dataAb.txt")




# For simplification, if an sid has NA value in any variable, then delate all the rows of this sid
dataset_new <- dataset %>%
  group_by(sid) %>%
  filter(Vxind == 1, HIVinfectionind == 0) %>%      # 2815
  filter(n() > 1) %>%                             # Keep only sids with more than 1 record     2759
  dplyr::select(sid, sampledays, doesdays, doesnumberschedule, NAb, MNGNE8, riskscorecat, GNE8_CD4, MN_CD4) %>%
  filter(!any(is.na(across(everything())))) %>%  # Remove rows with NA
  ungroup()     #2235




# find average interval of sampledays for each sample sid
avg_interval_df <- dataset_new %>%  
  group_by(sid) %>%
  arrange(sampledays) %>%
  summarise(
    n = n(),
    last_sampleday = max(sampledays, na.rm = TRUE),
    avg_interval = round(mean(diff(sampledays), na.rm = TRUE)))


# take 90% of the max sampledays as cutoff day and 80% of the itnerval day as average interval day
cutoff_day <- round(quantile(avg_interval_df$last_sampleday, probs = 0.9))  
average_interval <- round(quantile(avg_interval_df$avg_interval, probs = 0.8))


# dropped_out_interval as interval censored result
# dropped_out_right as right censored result
surv_data <- avg_interval_df %>%
  mutate(
    dropout_L = last_sampleday,
    dropout_R = case_when(
        last_sampleday > cutoff_day ~ last_sampleday,
       (last_sampleday + average_interval) <= cutoff_day ~ last_sampleday + average_interval, 
       (last_sampleday + average_interval) > cutoff_day ~ last_sampleday + average_interval,
    ),
    dropped_out_interval = ifelse( dropout_R <= cutoff_day, 1, 0),
    dropped_out_right = ifelse( last_sampleday <= cutoff_day, 1, 0)
  )


#----------------------------------------------------------------------------------------------------

# -----------------------------build long_data and surv_data----------------------------


#step 1: based on information, match the dose number and vaccine time period
vaccine_time_map <- c(0, 1, 6, 12, 18, 24, 30)

vac_schedule <- data.frame(
  dosenumberschedule = 1:7,
  vac_month = vaccine_time_map
)


dataset_new <- dataset_new %>%
  mutate(
    sample_month = sampledays / 30
  )



# step 2 : create a function to compute delta
get_deltas <- function(df) {
  vac_times <- vaccine_time_map
  last_vac_time <- tail(vac_times, 1)
  valid_delta_list <- list()
  
  for (j in 1:(length(vac_times) - 1)) {
    start <- vac_times[j]
    end <- vac_times[j + 1]
    
    # check if sample falls in between the interval of a vaccine range
    # if it is, save the dose number and delta_ij into valid_delta_list 
    
    if (any(df$sample_month >= start & df$sample_month < end)) {
      
      sel <- df$sample_month >= start & df$sample_month < end
      temp <- df[sel, ]
      valid_delta_list[[length(valid_delta_list) + 1]] <- 
        data.frame(
          sid = temp$sid,
          start_doesnumber = j,
          end_doesnumber = j + 1,
          sampledays = temp$sampledays,
          delta_ij = end - start
        )
    }}
    
    # for sample month that larger than the last dose month, we decide 
    # the time between the final vaccination and the final measurement time as delta
    sel <- df$sample_month >= last_vac_time
    if (any(sel)) {
      temp <- df[sel, ]
      valid_delta_list[[length(valid_delta_list) + 1]] <- 
        data.frame(
          sid = temp$sid,
          start_doesnumber = length(vac_times),  # the 7th dose number
          end_doesnumber = NA,                   # no more dose
          sampledays = temp$sampledays,
          delta_ij = (temp$sampledays - (last_vac_time * 30))/30 
        )
    }
    
    do.call(rbind, valid_delta_list)
  }


# step 3 : run function 
delta_info <- dataset_new %>%
  group_by(sid) %>%
  group_split() %>%
  lapply(get_deltas) %>%
  bind_rows()


# step 4 : add delta into orginal longitutional dataset
long_data <- dataset_new %>%
  left_join(delta_info, by = c("sid", "sampledays"))



# Step 5: other variables
long_data <- long_data %>%
  mutate(
    tdij = doesdays /30,
    t = sampledays/30,
    t_dij_star = tdij*30/7,
    t_star =  t/12,
    t_star2 = t_star^2,
    # for those delta = 0, we use 0.01 instead to compute sin term, such that we will not have NaNs produced
    sin_term = ifelse(delta_ij != 0, sin(pi * tdij / delta_ij), sin(pi * tdij / 0.01)), 
    risk1 = as.numeric(riskscorecat == 1),
    risk2 = as.numeric(riskscorecat == 2)
  )




# step 6 : add the nessary variable from long_data to surv_data (choose the baseline value)
baseline_vars <- long_data %>%
  group_by(sid) %>%
  slice_min(sampledays, with_ties = FALSE) %>%
  dplyr::select(sid, NAb, MNGNE8, riskscorecat, GNE8_CD4, MN_CD4) %>%
  mutate(
    risk1 = as.numeric(riskscorecat == 1),
    risk2 = as.numeric(riskscorecat == 2)
  )

surv_data_final <- surv_data %>%
  dplyr::select(sid, dropout_L = last_sampleday, dropout_R, dropped_out_interval, dropped_out_right) %>%
  mutate(
    L = dropout_L,
    R = dropout_R
  ) %>%
  left_join(baseline_vars, by = "sid")



# step 7: crete Cij and Zij value
threshold_nab <- quantile(long_data$NAb, 0.27, na.rm = TRUE)  # set up the censored value
long_data <- long_data %>%
  mutate(
    Cij = as.numeric(NAb <= threshold_nab)
  )

long_data <- long_data %>% # define Zij based on paper
  mutate(Zij = as.integer(MNGNE8 > 0.57))

#------------  Set up models--------------------


# (1) LME model 1
md1 <- lmer(
  NAb ~ 1 + t_star + sin_term + risk1 + risk2 + (1 | sid),
  data =long_data 
)


# get the estimated random effect by observation
estBi <- data.frame(row.names(ranef(md1)$sid), 
                    scale(ranef(md1)$sid, center=T,scale=T))
names(estBi) <- c("sid", "b11")
mydat <- merge(long_data, estBi, by='sid',all=T)


# (2) GLME model for NAb

fm2 <- Cij ~ 1 + t_star + sin_term + risk1 + risk2 + b11
md2 <- glm(fm2, data = mydat, family = binomial())

# (3) GLME model for Cij

fm3 <-  Zij ~ 1 + t + sin_term + t_dij_star + b11 + ( t - 1 | sid)
md3 <- glmer(fm3, data = mydat, family = binomial())
summary(md3)
Sdata <- surv_data_final
Sdata$b11 <- scale(ranef(md1)$sid[,1], center=T, scale=T)
Sdata$b21 <- scale(ranef(md3)$sid[,1], center=T, scale=T)



#---------------------right censored survival model----------------------
# a Cox PH model
fitCOX1 <- coxph(Surv(L, dropped_out_right) ~ GNE8_CD4 + risk1 + risk2 + b11 + b21, data = Sdata)   

# a Weibull model
fitCOX2 <- survreg(Surv(L, dropped_out_right) ~ GNE8_CD4 + risk1 + risk2 + b11 + b21, data = Sdata, dist='weibull')


############################################
######   Joint modeling  right censored#####
############################################



#-------------(1) Create objects for longitudinal models---------------
LLOQ = 2
# right censored part
CenObject <- list(
  fm = Cij ~ t_star + sin_term + risk1 + risk2 + b11,
  family='binomial', par='eta', ran.par="b11",
  disp="eta5",
  lower = -Inf,
  upper = Inf,
  str_val=coef(md2),
  Cregime=1,
  truncated=T, delim_val=LLOQ)



# continuous data
glmeObject1 <- list(
  fm = NAb ~ t_star + sin_term + risk1 + risk2 + b11,
  family='normal', par="beta", ran.par='b11', sigma='sigma',
  disp='beta5',
  lower = 0,
  upper = Inf,
  str_val=c(fixef(md1), sd(ranef(md1)$sid[,1])),
  CenObject=CenObject)





# binary data
glmeObject2 <- list(
  fm= Zij ~ 1 + t + sin_term + t_dij_star + b11 + b21*t,
  family='binomial', par="alpha", ran.par=c('b11','b21'),
  sigma=NULL, 
  disp=c("alpha4", 'alpha5'),
  lower=c(-Inf,0),
  upper=rep(Inf,2),
  str_val=c(fixef(md3), sd(ranef(md3)$sid[,1])),
  CenObject=NULL)


glmeObject <- list(glmeObject1, glmeObject2)




### (2) Create object for survival model
## case 1: if a Cox PH model is fit for survival data
survObject1 <- list(
  fm= L ~ GNE8_CD4 + risk1 + risk2 + b11 + b21,
  event="dropped_out_right", 
  par='lambda',
  disp=NULL,
  lower=NULL, upper=NULL,
  distribution=NULL,
  str_val= summary(fitCOX1)$coeff[,1])

## case 2: if a Weibull model is fit for survival data
survObject2 <- list(
  fm= L ~ GNE8_CD4 + risk1 + risk2 + b11 + b21,
  event="dropped_out_right", 
  par='lambda',
  disp=NULL,
  lower=c(0, 0), upper=c(Inf, Inf),
  distribution='weibull',
  str_val= summary(fitCOX2)$coeff[-1])


# ------------------ Cox model h-likelihood ----------------------
tic()
testjm1 <- try(JMfit(glmeObject, survObject1, 
                     long_data, surv_data_final,
                     idVar="sid", eventTime="L",
                     survFit=fitCOX1,
                     method = "h-likelihood"), silent=T)

new_sd1 = JMsd_aGH(testjm1, ghsize=4, srcpath=srcpath, paralle=T)
ptm <- toc()
(ptm$toc-ptm$tic)/60   # takes 3.33 min
# return coefficient table of Fixed effects 
JMsummary(testjm1)
JMsummary(testjm1, newSD=new_sd1)

#-----------------without drop out---------------------

#censored right eta0-eta4
md_eta <- glmer(Cij ~ t_star + sin_term + risk1 + risk2  + (1 | sid),
                data = mydat, family = binomial())
fixef(md_eta) 

# continuous beta0 - beta 4
md_beta <- lmer(NAb ~ t_star + sin_term + risk1 + risk2 + (1 | sid),
                data = mydat)
fixef(md_beta)

# binary alpha0 - alpha 3
md_bin <- glmer(Zij ~ 1 + t + sin_term + t_dij_star + (1 + t | sid),
             data = mydat, family = binomial())
fixef(md_bin)

# ----------------Weibull H-likelihood--------------------------------

tic()
testjm2 <- try(JMfit(glmeObject, survObject2, 
                     long_data, surv_data_final,
                     idVar="sid", eventTime="L",
                     survFit=fitCOX2,
                     method = "h-likelihood", Silent = F), silent=T)

new_sd2 = JMsd_aGH(testjm2, ghsize=4, srcpath=srcpath, paralle=T)
ptm <- toc()
(ptm$toc-ptm$tic)/60   # takes 3.33 min
# return coefficient table of Fixed effects 
JMsummary(testjm2)
JMsummary(testjm2, newSD=new_sd2)




#------------------Weibull Interval Censored model--------------------


surv_model <- survreg(
  Surv(time = L, time2 = R, type = "interval2") ~ GNE8_CD4 + risk1 + risk2 + b11 + b21,
  data = Sdata,
  dist = "weibull"
)





survObject <- list(
  fm = Surv(time = L, time2 = R, type = "interval2") ~ GNE8_CD4 + risk1 + risk2 + b11 + b21,
  event = "interval",             
  par = 'lambda',  # coefficients
  disp = c("lambda_base", "phi"),  # true scale & shape
  lower = c(0, 0),
  upper = rep(Inf, 2),
  distribution = "weibull_interval",
  str_val = str_val_surv,
  timegrid = timegrid_vec,
  rnames = paste0("r", 1:(length(timegrid_vec) - 1))
)



#-----------------------interval censored survival model---------------------------------
surv_model <- survreg(
  Surv(time = L, time2 = R, type = "interval2") ~ GNE8_CD4 + risk1 + risk2 + b11 + b21,
  data = Sdata,
  dist = "weibull"
)




#-------------(1) Create objects for longitudinal models---------------

# right censored part
CenObject <- list(
  fm = Cij ~ t_star + t_star2 + sin_term + risk1 + risk2 + b11,
  family = "binomial",
  par = "eta",
  ran.par = "b11",
  disp = "eta_disp",
  lower = -Inf,
  upper = Inf,
  str_val = coef(md2),
  Cregime = 1,             # right-censored
  truncated = TRUE,        # censored model
  delim_val = threshold_nab
)




# continuous data
glmeObject1 <- list(
  fm = NAb ~ t_star + t_star2 + sin_term + risk1 + risk2 + b11,
  family = "normal",
  par = "beta",
  ran.par = "b11",
  sigma = "sigma",
  disp = "beta_disp",
  lower = rep(-Inf, 6),
  upper = rep(Inf, 6),
  str_val = c(fixef(md1), sd(ranef(md1)$sid[, 1])),
  CenObject = CenObject
)





# binary data
glmeObject2 <- list(
  fm= Zij ~ 1 + t + sin_term + t_dij_star + b11 + b21,
  family='binomial', 
  par="alpha", 
  sigma=NULL, 
  disp=c("alpha4", 'alpha5'),
  lower=c(-Inf,0),
  upper=rep(Inf,2),
  str_val=c(fixef(md3), sd(ranef(md3)$sid[,1])),
  CenObject=NULL)


glmeObject <- list(glmeObject1, glmeObject2)


#-------------(1) Create objects for survival models---------------
timegrid_vec <- seq(min(surv_data$L, na.rm = TRUE), max(surv_data$R, na.rm = TRUE), by = 100)


for (i in 1:(length(timegrid_vec) - 1)) {
  varname <- paste0("r", i)
  Sdata[[varname]] <- as.numeric(Sdata$L <= timegrid_vec[i] & Sdata$R > timegrid_vec[i])
}


str_val_surv <- c(
  -summary(surv_model)$coeff[-1] / summary(surv_model)$scale,
  exp(summary(surv_model)$coeff[1]),
  1 / summary(surv_model)$scale
)
names(str_val_surv) <- c("lambda0", "lambda1", "lambda2", "lambda3", "lambda4", "lambda_base", "phi")




survObject <- list(
  fm = Surv(time = L, time2 = R, type = "interval2") ~ GNE8_CD4 + risk1 + risk2 + b11 + b21,
  event = "interval",             
  par = 'lambda',  # coefficients
  disp = c("lambda_base", "phi"),  # true scale & shape
  lower = c(0, 0),
  upper = rep(Inf, 2),
  distribution = "weibull_interval",
  str_val = str_val_surv,
  timegrid = timegrid_vec,
  rnames = paste0("r", 1:(length(timegrid_vec) - 1))
)




#----------------Joint modeling using h-likelihood method -------------------

long_data <- long_data %>%
  left_join(b11_df, by = "sid") %>%
  left_join(b21_df, by = "sid")



long_data <- long_data %>%
  filter(!is.na(b11), !is.na(b21))

Sdata <- Sdata %>%
  filter(!is.na(b11), !is.na(b21))



debug(estRaneff)

jm_result <- try(JMfit(
  glmeObject = glmeObject,
  survObject = survObject,
  long.data = long_data,
  surv.data = Sdata,
  idVar = "sid",
  eventTime = "R",             # ?
  survFit = surv_model,      
  method = "h-likelihood"
), silent=T)





bad_sids <- c()
unique_ids <- unique(long_data$sid)

for (sid_i in unique_ids) {
  subdat <- subset(long_data, sid == sid_i)
  subsurv.dat <- subset(Sdata, sid == sid_i)
  
  result <- try({
    estRaneff(
      RespLog = RespLog, 
      raneff = Jraneff, 
      long.data = subdat,
      surv.dat = subsurv.dat,
      invSIGMA = invSIGMA0,
      sigma = sigma0,
      ParVal = fixedest0,
      Silent = TRUE
    )
  }, silent = TRUE)
  
  if (inherits(result, "try-error")) {
    cat("❌ Error with sid:", sid_i, "\n")
    bad_sids <- c(bad_sids, sid_i)
  } else {
    cat("✅ OK for sid:", sid_i, "\n")
  }
}
clean_long_data <- long_data[!(long_data$sid %in% bad_sids), ]
clean_Sdata <- Sdata[!(Sdata$sid %in% bad_sids), ]



#--------------simulation-----------------------




# Step 1: set parameter
n_sims_large <- 10
n_sims <- 1000
n_subjects <- 100
param_names <- rownames(JMsummary(testjm1))
estimates <- matrix(NA, nrow = length(rownames(JMsummary(testjm1))), ncol = n_sims_large,
                    dimnames = list(param_names, NULL))
std_errors <- matrix(NA, nrow = length(rownames(JMsummary(testjm1))), ncol = n_sims_large,
                     dimnames = list(param_names, NULL))
error_count <- 0


# Cox version parameter
lambda_pars <- c(2.237, 2.094, 1.176, -3.933, -1.430)  # GNE8_CD4, risk1, risk2, b11, b21
beta_pars   <- c(0.421, 1.903, 0.622, 0.141, -0.272)   # intercept, t_star, sin, risk1, risk2
eta_pars    <- c(-1.894, 0.943, 0.509, 1.337, -0.002)  # intercept, t_star, sin, risk1, risk2, b11
alpha_pars  <- c(-1.762, 0.122, 1.773, -0.026)         # intercept, t, sin, t_dij_star, b11
lambda_pars <- lambda_pars*7

# get shape from real data
# fit_wb <- survreg(Surv(L, dropped_out_right) ~ GNE8_CD4 + risk1 + risk2 + b11 + b21,
#                  data = Sdata, dist = "weibull")

# Weibull -> scale = 1 / shape
# shape_phi <- 1 / fit_wb$scale # Weibull shape
# logscale <- log(fit_wb$scale) # Weibull scale
 shape_phi <- 100
 logscale <- 	-450


# Weibull survival model
simulate_weibull <- function(lp, shape, logscale, n) {
  # lp = X^T beta
  # shape = k, logscale = log(lambda)
  # S(t) = exp( - (lambda * t)^k ) = exp( - exp(logscale + log(t^k)) )
  # inverse: T = ( -log(U) / lambda )^(1/k)
  scale <- exp(logscale)  # lambda
  u <- runif(n)
  T <- (-log(u) / (scale * exp(lp)))^(1 / shape)
  return(T)
}


for (i in 1:n_sims_large) {
  cat("Simulation #", i, "\n")
  
  
  
for (i in 1:n_sims) {
  
 
  
  
  # ---- Step 1: simulate longitudinal data (1–15th random number of t for each sid) ----
  long_data_simulation <- purrr::map_dfr(1:n_subjects, function(sid) {
    n_obs <- sample(2:15, 1)
    real_t_density <- density(long_data$t)
    positive_idx <- real_t_density$x >= 0
    t <- sort(sample(
      real_t_density$x[positive_idx],
      size = n_obs,
      prob = real_t_density$y[positive_idx],
      replace = F
    ))
    
    # t <-sort(round(60 * rbeta(n_obs, shape1 = 2.1, shape2 = 5.5)))
    
   
    real_d_density <- density(long_data$t_dij_star)
    positive_dx <- real_d_density$x >= 0
    
    
    dose_times <- c(0,1,6,12,18,24,30)
    A_ik <- 0.05 + 0.05 * (1:length(dose_times))
    lambda_i <- runif(1, 0.1, 0.15)               
    t_rise <- 2  
    
    
    GNE8_CD4 <- sapply(t, function(ti) {
      sum(sapply(1:length(dose_times), function(k) {
        if (ti >= dose_times[k]) {
          A_ik[k] * exp(-lambda_i * (ti - dose_times[k]))
        } else {
          0
        }
      })) + rnorm(1, 0, 0.02)
    })
    
    
#    t_dij_star <- sort(sample(real_d_density$x, size = n_obs, prob = real_d_density$y, replace = TRUE))
    data.frame(
      t = t,
      sid = sid,
      t_star = t / 12,
      t_dij_star = sort(sample(real_d_density$x[positive_dx], size = n_obs, prob = real_d_density$y[positive_dx], replace = TRUE)),
#      t_dij_star = sort(2 * rbeta(n_obs, shape1 = 1.5, shape2 = 5)),
      risk1 = rbinom(n_obs, 1, 0.4),
      risk2 = rbinom(n_obs, 1, 0.02),
      GNE8_CD4 =  GNE8_CD4
    )
  })
  
  
  # get sin_term
  
  
  get_deltas <- function(df) {
    vac_times <- c(0, 1, 6, 12, 18, 24, 30)
    last_vac_time <- tail(vac_times, 1)
    valid_delta_list <- list()
    
    for (j in 1:(length(vac_times) - 1)) {
      start <- vac_times[j]
      end <- vac_times[j + 1]
      if (any(df$t >= start & df$t < end)) {
        temp <- df[df$t >= start & df$t < end, ]
        valid_delta_list[[length(valid_delta_list) + 1]] <- 
          data.frame(
            sid = temp$sid,
            t = temp$t,
            delta_ij = end - start
          )
      }
    }
    
    sel <- df$t >= last_vac_time
    if (any(sel)) {
      temp <- df[sel, ]
      valid_delta_list[[length(valid_delta_list) + 1]] <- 
        data.frame(
          sid = temp$sid,
          t = temp$t,
          delta_ij = temp$t - last_vac_time
        )
    }
    
    do.call(rbind, valid_delta_list)
  }
  
  delta_df <- long_data_simulation %>%
    group_by(sid) %>%
    group_split() %>%
    lapply(get_deltas) %>%
    bind_rows()
  
  long_data_simulation <- long_data_simulation %>%
    left_join(delta_df, by = c("sid", "t"))
  
  
  
  long_data_simulation <- long_data_simulation %>%
    mutate(
      sin_term = ifelse(delta_ij != 0,
                        sin(pi * t_dij_star*7/30 / delta_ij),
                        sin(pi * t_dij_star*7/30 / 0.01))
    )
  
  
  # ---- Step 2: add random effect b11, b21 ----
  bi <- mvrnorm(n_subjects, mu = c(0, 0), Sigma = matrix(c(0.8, 0.3, 0.3, 0.8), 2))
  ranef_df <- data.frame(sid = 1:n_subjects, b11 = bi[, 1], b21 = bi[, 2])
  long_data_simulation <- long_data_simulation %>% left_join(ranef_df, by = "sid")
  
  # ---- Step 3: simulate event time from Weibull ----
  surv_data_simulation <- long_data_simulation %>%
    group_by(sid) %>%
    summarise(
      GNE8_CD4 = first(GNE8_CD4),
      risk1 = first(risk1),
      risk2 = first(risk2),
      b11 = first(b11),
      b21 = first(b21)
    )
  
  
  lp_surv <- with(surv_data_simulation,
                  lambda_pars[1]*GNE8_CD4 + lambda_pars[2]*risk1 + lambda_pars[3]*risk2 +
                    lambda_pars[4]*b11 + lambda_pars[5]*b21
  )
  
  T_surv <- simulate_weibull(lp_surv, shape_phi, logscale, n_subjects)
  
  
  # ---- Step 4: assign true event time ----
  surv_data_simulation <- surv_data_simulation %>%
    mutate(eventtime = T_surv)
  
  # ---- Step 5: set administrative censoring cutoff (e.g. 90th percentile) ----
  cutoff <- quantile(T_surv, 0.9)
  
  # ---- Step 6: define observed time L and censoring indicator ----
  surv_data_simulation <- surv_data_simulation %>%
    mutate(
      L = pmin(eventtime, cutoff),  # observed follow-up time (may be censored)
      dropped_out_right = as.integer(eventtime >= cutoff)  # 1 = censored, 0 = event
    )
  

  
  long_data_simulation <- long_data_simulation %>%
    left_join(surv_data_simulation %>% dplyr::select(sid, L), by = "sid") 
  

  
  # ---- Step 8: simulate outcome（NAb, Cij, Zij） ----
  X_long <- model.matrix(~ t_star + sin_term + risk1 + risk2 , data = long_data_simulation)
  long_data_simulation$NAb <- as.numeric(X_long %*% beta_pars + long_data_simulation$b11 + rnorm(nrow(X_long), sd = 0.5))
  
  X_c <- model.matrix(~ t_star + sin_term + risk1 + risk2, data = long_data_simulation)
  linpred_c <- as.numeric(X_c %*% eta_pars + long_data_simulation$b11)  # same as random intercept
  long_data_simulation$Cij <- rbinom(nrow(X_c), 1, plogis(linpred_c))
  
  X_z <- model.matrix(~ t + sin_term + t_dij_star, data = long_data_simulation)
  linpred_z <- as.numeric(X_z %*% alpha_pars + long_data_simulation$b11 + long_data_simulation$t * long_data_simulation$b21)
  long_data_simulation$Zij <- rbinom(nrow(X_z), 1, plogis(linpred_z))
  
  # ---- Step 6: save result ----
  write_csv(long_data_simulation, paste0("sim_data_long_", i, ".csv"))
  write_csv(surv_data_simulation, paste0("sim_data_surv_", i, ".csv"))
}



#-------------------------testing----------------



# (1) LME model 1
md1_new <- lmer(
  NAb ~ 1 + t_star + sin_term + risk1 + risk2 + (1 | sid),
  data =long_data_simulation 
)


# get the estimated random effect by observation
estBi <- data.frame(row.names(ranef(md1_new)$sid), 
                    scale(ranef(md1_new)$sid, center=T,scale=T))
names(estBi) <- c("sid", "esb11")
mydat_new <- merge(long_data_simulation, estBi, by='sid',all=T)


# (2) GLME model for NAb

fm2_new <- Cij ~ 1 + t_star + sin_term + risk1 + risk2 + esb11
md2_new <- glm(fm2_new, data = mydat_new, family = binomial())

# (3) GLME model for Cij

fm3_new <-  Zij ~ 1 + t + sin_term + t_dij_star + esb11 + ( t - 1 | sid)
md3_new <- glmer(fm3_new, data = mydat_new, family = binomial())
summary(md3_new)
Sdata_new <- surv_data_simulation
Sdata_new$esb11 <- scale(ranef(md1_new)$sid[,1], center=T, scale=T)
Sdata_new$esb21 <- scale(ranef(md3_new)$sid[,1], center=T, scale=T)



#---------------------right censored survival model----------------------
# a Cox PH model
fitCOX1_new <- coxph(Surv(L, dropped_out_right) ~ GNE8_CD4 + risk1 + risk2 + esb11 + esb21, data = Sdata_new)   

# a Weibull model
fitCOX2_new <- survreg(Surv(L, dropped_out_right) ~ GNE8_CD4 + risk1 + risk2 + esb11 + esb21, data = Sdata_new, dist='weibull')


############################################
######   Joint modeling  right censored#####
############################################



#-------------(1) Create objects for longitudinal models---------------
LLOQ = 2
# right censored part
CenObject_new <- list(
  fm = Cij ~ t_star + sin_term + risk1 + risk2 + esb11,
  family='binomial', par='eta', ran.par="esb11",
  disp="eta5",
  lower = -Inf,
  upper = Inf,
  str_val=coef(md2_new),
  Cregime=1,
  truncated=T, delim_val=LLOQ)



# continuous data
glmeObject1_new <- list(
  fm = NAb ~ t_star + sin_term + risk1 + risk2 + esb11,
  family='normal', par="beta", ran.par='esb11', sigma='sigma',
  disp='beta5',
  lower = 0,
  upper = Inf,
  str_val=c(fixef(md1_new), sd(ranef(md1_new)$sid[,1])),
  CenObject=CenObject_new)





# binary data
glmeObject2_new <- list(
  fm= Zij ~ 1 + t + sin_term + t_dij_star + esb11 + esb21*t,
  family='binomial', par="alpha", ran.par=c('esb11','esb21'),
  sigma=NULL, 
  disp=c("alpha4", 'alpha5'),
  lower=c(-Inf,0),
  upper=rep(Inf,2),
  str_val=c(fixef(md3_new), sd(ranef(md3_new)$sid[,1])),
  CenObject=NULL)


glmeObject_new <- list(glmeObject1_new, glmeObject2_new)




### (2) Create object for survival model
## case 1: if a Cox PH model is fit for survival data
survObject1_new <- list(
  fm= L ~ GNE8_CD4 + risk1 + risk2 + esb11 + esb21,
  event="dropped_out_right", 
  par='lambda',
  disp=NULL,
  lower=NULL, upper=NULL,
  distribution=NULL,
  str_val= summary(fitCOX1_new)$coeff[,1])

## case 2: if a Weibull model is fit for survival data
survObject2_new <- list(
  fm= L ~ GNE8_CD4 + risk1 + risk2 + esb11 + esb21,
  event="dropped_out_right", 
  par='lambda',
  disp=NULL,
  lower=c(0, 0), upper=c(Inf, Inf),
  distribution='weibull',
  str_val= summary(fitCOX2_new)$coeff[-1])


# ------------------ Cox model h-likelihood ----------------------
tic()
testjm1_new <- try(JMfit(glmeObject_new, survObject1_new, 
                     long_data_simulation, surv_data_simulation,
                     idVar="sid", eventTime="L",
                     survFit=fitCOX1_new,
                     method = "h-likelihood", Silent = F), silent=T)



#---------------------------

if (inherits(testjm1_new, "try-error")) {
  error_count <- error_count + 1
  next
}

coef_vec <- try(JMsummary(testjm1_new), silent = TRUE)
se_vec <- try(JMsd_aGH(testjm1_new, ghsize=4, srcpath=srcpath, paralle=TRUE), silent = TRUE)

if (inherits(coef_vec, "try-error") || inherits(se_vec, "try-error")) {
  error_count <- error_count + 1
  next
}
for (j in param_names) {
  if (!is.na(coef_vec[j]) && !is.na(se_vec[j])) {
    estimates[j, i] <- coef_vec[j]
    std_errors[j, i] <- se_vec[j]
  }
}
}

#--------------Plot-----------------
ggplot(long_data, aes(x = t, y =  GNE8_CD4, group = sid)) +
  geom_line(color = "black", alpha = 0.5) +   
  geom_point(color = "black", size = 1) + 
  labs(title = "GNE8_CD4 longitudinal trajectories", y = "GNE8_CD4", x = "month") +
  theme_minimal()



ggplot(long_data_simulation, aes(x = t, y = GNE8_CD4, group = sid)) +
  geom_line(color = "black", alpha = 0.5) +   
  geom_point(color = "black", size = 1) + 
  geom_hline(yintercept = 2, linetype = "solid") +
  annotate("text", x = 10, y = 2 + 0.1, 
           label = "Lower limit of quantification (LLOQ)", hjust = 0)  +
  labs(title = "GNE8_CD4 longitudinal trajectories", y = "GNE8_CD4", x = "month") +
  theme_minimal()




long_data$type <- "Real"
long_data_simulation$type <- "simulation"

combined_data <- bind_rows(long_data, long_data_simulation)

combined_data <- combined_data %>%
  arrange(type, t_dij_star) %>%
  group_by(type) %>%
  mutate(rank = row_number())

ggplot(combined_data, aes(x = t_dij_star, fill = type)) +
  geom_density(alpha = 0.5) +
  xlim(0, 100)+
  labs(title = "Comparison of t (mesurment time)",
       x = "Participant Index (Sorted)", y = "t (days)",
       color = "Dataset") +
  theme_minimal()



surv_data_final$type <- "Real"
surv_data_simulation$type <- "simulation"

combined_data <- bind_rows(surv_data_final, surv_data_simulation)

combined_data <- combined_data %>%
  arrange(type, L) %>%
  group_by(type) %>%
  mutate(rank = row_number())

ggplot(combined_data, aes(x = L, fill = type)) +
  geom_density(alpha = 0.5) +
  xlim(0, 5000) +
  labs(title = "Comparison of L (Observation Time)",
       x = "L", y = "density",
       color = "Dataset") +
  theme_minimal()







ggplot(surv_data, aes(x = last_sampleday/30)) +        
  geom_histogram(aes(y = ..density..),
                 binwidth = 1,
                 fill = "steelblue", color = "white", alpha = 0.7) +
  geom_density(color = "red", size = 1) +

  geom_vline(xintercept = cutoff_day/30,
             linetype = "dashed", color = "darkred", size = 1) +
  labs(x = "Last-visit time (months)",
       y = "Density",
       title = "Empirical distribution of last-visit times\nwith cutoff day") +
  theme_minimal()

