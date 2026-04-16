library(tidyverse)
library(readr)
library(lme4)
library(tictoc)
library(survival)
library(MASS)
library(glmmTMB)
library(Matrix)
library(ggplot2)
library(flexsurv)

## source all the R code
srcpath = "E:/UBC/Final Project/New folder/Final-Project/orin/HHJMs-master/R"
setwd(srcpath)
(file.sources = list.files(pattern="\\.[rR]$"))

sapply(file.sources,source,.GlobalEnv)

# -------------- decide interval censored as drop out---------------

dataset <- read_csv("E:/UBC/Final Project/vax004dataAb.txt")




# For simplification, if an sid has NA value in any variable, then delate all the rows of this sid
dataset_new <- dataset %>%
  group_by(sid) %>%
  filter(Vxind == 1, HIVinfectionind == 0) %>%      # 2815
  filter(n() > 1) %>%                             # Keep only sids with more than 1 record     2759
  dplyr::select(sid, sampledays, doesdays,peakvalley, doesnumberschedule, NAb, MNGNE8, riskscorecat, GNE8_CD4, MN_CD4) %>%
  filter(!any(is.na(across(everything())))) %>%  # Remove rows with NA
  ungroup()     #2235


# -----------------------------------
# valley-only planned days (used to impute R)
# doesnumberschedule 1~7 corresponds to Month 0,1,6,12,18,24,30
# -----------------------------------
valley_grid <- c(0, 30, 180, 360, 540, 720, 900)  # unit: days

# -----------------------------------
# population-level: mean actual visit day per doesnumberschedule (valley only)
# -----------------------------------
dose_avg_days <- dataset_new %>%
  filter(peakvalley == 0) %>%
  group_by(doesnumberschedule) %>%
  summarise(avg_sampleday = mean(sampledays, na.rm = TRUE), .groups = "drop")

# -----------------------------------
# subject-level: systematic visit bias
# (individual actual day - population mean for that dose)
# -----------------------------------
individual_bias <- dataset_new %>%
  filter(peakvalley == 0) %>%
  left_join(dose_avg_days, by = "doesnumberschedule") %>%
  group_by(sid) %>%
  summarise(
    visit_bias = mean(sampledays - avg_sampleday, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------
# dropout flag + L
# dropout defined as: did not complete all 7 scheduled valley visits
# L = last observed visit day (peak or valley)
# -----------------------------------
avg_interval_df <- dataset_new %>%
  group_by(sid) %>%
  arrange(sampledays) %>%
  summarise(
    n = n(),
    last_sampleday = max(sampledays, na.rm = TRUE),
    avg_interval = round(mean(diff(sampledays), na.rm = TRUE)),
    max_dose_valley = max(doesnumberschedule[peakvalley == 0], na.rm = TRUE),
    .groups = "drop"
  )

cutoff_day <- round(quantile(avg_interval_df$last_sampleday, probs = 0.4))  
# -----------------------------------
# construct survival data
# -----------------------------------
surv_data <- avg_interval_df %>%
  left_join(individual_bias, by = "sid") %>%
  mutate(
    dropout_L = last_sampleday,
    next_dose = max_dose_valley + 1,
    R_base = case_when(
      next_dose <= 7 ~ dose_avg_days$avg_sampleday[next_dose],
      TRUE ~ Inf
    ),
    dropout_R = ifelse(is.finite(R_base), R_base + visit_bias, Inf),
    dropout_R = pmax(dropout_R, dropout_L + 1),
    dropout_R = ifelse(is.finite(R_base), ceiling(dropout_R), Inf),
    
    dropped_out_right    = ifelse(last_sampleday <= cutoff_day, 1, 0),
    dropped_out_interval = dropped_out_right,
    
    # for people who complete the survey but event=1（R=Inf）use individual average time
    dropout_R = case_when(
      dropped_out_right == 1 & is.finite(dropout_R)  ~ dropout_R,
      dropped_out_right == 1 & !is.finite(dropout_R) ~ as.numeric(last_sampleday +
                                                                    avg_interval),
      TRUE ~ Inf
    )
  )

#----------------------------------------------------------------------------------------------------

# -----------------------------build long_data and surv_data----------------------------


#step 1: based on information, match the dose number and vaccine time period
vaccine_time_map <- c(0, 1, 6, 12, 18, 24, 30)

vac_schedule <- data.frame(
  doesnumberschedule = 1:7,
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
        end_doesnumber = 8,                   # no more dose
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
    L = dropout_L/30,
    R = dropout_R/30
  ) %>%
  left_join(baseline_vars, by = "sid")



# step 7: crete Cij and Zij value
quantile(long_data$NAb, 0.27, na.rm = TRUE)
threshold_nab <- min(long_data$NAb)
long_data <- long_data %>%
  mutate(
    Cij = as.numeric(NAb <= threshold_nab)
  )

long_data <- long_data %>% # define Zij based on paper
  mutate(Zij = as.integer(MNGNE8 > 0.57))




#------------  Set up models--------------------


# (1) LME model 
md1 <- lmer(
  NAb ~  1 + t_star + sin_term + risk1 + risk2 + (1 | sid),
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

fm3 <-  Zij ~ 1 + t  + sin_term + t_dij_star + b11 + ( t - 1 | sid)
md3 <- glmer(fm3, data = mydat, family = binomial())

Sdata <- surv_data_final
Sdata$b11 <- scale(ranef(md1)$sid[,1], center=T, scale=T)
Sdata$b21 <- scale(ranef(md3)$sid[,1], center=T, scale=T)




#---------------------right censored survival model----------------------


# a Cox PH model
fitCOX1 <- coxph(Surv(L, dropped_out_right) ~ GNE8_CD4 + risk1 + risk2 + b11 + b21, data = Sdata)   

# a Weibull model
fitCOX2 <- survreg(Surv(L, dropped_out_right) ~ GNE8_CD4  + risk1 + risk2 + b11 + b21, data = Sdata, dist='weibull')


############################################
######   Joint modeling  right censored#####
############################################


#-------------(1) Create objects for longitudinal models---------------
LLOQ = threshold_nab
# right censored part
CenObject <- list(
  fm = Cij ~ 1 + t_star  + sin_term + risk1 + risk2 + b11,
  family='binomial', par='eta', ran.par="b11",
  disp="eta5",
  lower = -Inf,
  upper = Inf,
  str_val=coef(md2),
  Cregime=1,
  truncated=T, delim_val=LLOQ)



# continuous data
glmeObject1 <- list(
  fm = NAb ~ 1 + t_star  + sin_term + risk1 + risk2 + b11,
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
  lower=NULL, upper=NULL,
  distribution='weibull',
  str_val= -summary(fitCOX2)$coeff[-1]/summary(fitCOX2)$scale)


# ------------------ Cox model h-likelihood ----------------------
tic()

# joint model without LLOQ part
set.seed(3)
testjm_without <- try(JMfit(glmeObject, survObject1, # the glmeObject1 with CenObject=NULL
                            long_data, surv_data_final,
                            idVar="sid", eventTime="L",
                            survFit=fitCOX1,
                            method = "h-likelihood", Silent=T), silent=F)
JMsummary(testjm_without)
new_sd1 = JMsd_aGH(testjm_without, ghsize=4, srcpath=srcpath, paralle=T)
JMsummary(testjm_without, new_sd1$sd_naive)

set.seed(3)
testjm12  <- try(JMfit(glmeObject, survObject1, 
                      long_data, surv_data_final,
                      idVar="sid", eventTime="L",
                      survFit=fitCOX1,
                      method = "h-likelihood", Silent=T), silent=F) 


testjm1$RespLog

# return coefficient table of Fixed effects 
JMsummary(testjm1)
JMsummary(testjm1, newSD=new_sd1)

# random effect
testjm1$sigma

testjm1$covBi



# ----------------Weibull H-likelihood--------------------------------

tic()
set.seed(1)
testjm2  <- JMfit(glmeObject, survObject2, 
                  long_data, surv_data_final,
                  idVar="sid", eventTime="L",
                  survFit=fitCOX2,
                  method = "h-likelihood", Silent = T) #itertol = 2e-2

new_sd2 = JMsd_aGH(testjm2, ghsize=4, srcpath=srcpath, paralle=T)
ptm <- toc()
(ptm$toc-ptm$tic)/60   # takes 3.33 min
# return coefficient table of Fixed effects 
JMsummary(testjm2)
JMsummary(testjm2, newSD=new_sd2)



#---------------------------------------------------------------------
#------------------Weibull Interval Censored model--------------------
#---------------------------------------------------------------------

#-----------------------interval censored survival model---------------------------------

# AFT Weibull Model
surv_model <- survreg(
  Surv(time = L, time2 = R, type = "interval2") ~ GNE8_CD4 + risk1 + risk2 + b11 + b21,
  data = Sdata,
  dist = "weibull"
)


# PH Weibull Model
surv_model_ph <- flexsurvreg(
  Surv(L, R, type = "interval2") ~ GNE8_CD4 + risk1 + risk2 + b11 + b21,
  data = Sdata,
  dist = "weibullPH"   
)


# set a relatively large number for R = inf
R_f <- surv_data_final$R[is.finite(surv_data_final$R)]
R_big <- max(R_f, na.rm=TRUE) + 6
surv_data_final$R[is.infinite(surv_data_final$R)] <- R_big


#-------------(1) Create objects for survival models---------------


survObject_intPH <- list( fm = Surv(L, R) ~ GNE8_CD4 + risk1 + risk2 + b11 + b21, 
                          event = "dropped_out_interval", 
                          par = "lambda", # betas: lambda0, lambda1, ... 
                          disp = NULL, 
                          lower=c(0, 0), 
                          upper=c(Inf, Inf), 
                          distribution = "weibull_ph_interval", 
                          str_val = -summary(surv_model)$coeff[-1]/summary(surv_model)$scale) 

set.seed(2)
K1<- JMfit(
  glmeObject = glmeObject, survObject_intPH,
  long.data = long_data,
  surv.data = surv_data_final,
  idVar = "sid",
  eventTime = NULL,      # no need
  survFit = surv_model,    
  method = "h-likelihood", Silent = T
)
JMsummary(K1)


# ----------------NLME + Right censored ---------------
long_df <- as.data.frame(dataset_new_NLME)
long_df$Zij <- long_data$Zij
long_df$sin_term <- long_data$sin_term
long_df$t_dij_star <- long_data$t_dij_star

# add b_lambda into long data 
b_11 <- as.numeric(scale( ranef(fit_nlme)[, "lambda"], center = TRUE, scale = TRUE))

b_11_add <- data.frame(
  sid = rownames(ranef(fit_nlme) ),
  b_11 = b_11)
mydat_nlme <- merge(long_df, b_11_add , by = "sid", all.x = TRUE)

# model 2
fm3_same <- Zij ~ 1 + t + sin_term + t_dij_star + b_11 + (t - 1 | sid)
md3_same <- glmer(fm3_same, data = mydat_nlme, family = binomial())


# model from nlme_code
Sdata_nlme <- Sdata
Sdata_nlme$b_11 <- scale(ranef(fit_nlme)[, "lambda"], center = TRUE, scale = TRUE)
# Sdata_nlme$b_21 <- scale(ranef(fit_nlme)[, "A2"], center = TRUE, scale = TRUE)
Sdata_nlme$b_21 <- scale(ranef(md3_same)$sid[,1], center=T, scale=T)



fitCOX1_nlme <- coxph(Surv(L, dropped_out_right) ~ GNE8_CD4 + risk1 + risk2 + b_11 +b_21, data = Sdata_nlme)   



nlmeObject <- list(
  type = "nlme",
  family = "EXP_Delay",
  fm = nl_formula,
  sigma = "sigma",
  disp = "sigma_b11",                          # dispersion of standardised b_11
  lower = 0,                                   # sigma_b11 >= 0
  upper = Inf,
  A_names = paste0("A",1:7),
  ind_names = paste0("ind_k",1:7),
  lag_names = paste0("lag_k",1:7),
  lambda_name = "lambda",
  ran_map = list(
    lambda = "b_11"
  ),
  re_disp = list(
    b_11 = "sigma_b11"                         # lambda_i = exp(lambda + sigma_b11 * b_11)
  ),
  str_val = c(as.list(fixef(fit_nlme)),
              sd(ranef(fit_nlme)[["lambda"]]),  # sigma_b11 starting value (~0.0038)
              fit_nlme$sigma),                  # sigma (NAb measurement error) starting value
  CenObject = NULL
)


nlmeObject$str_val$lambda <- log(nlmeObject$str_val$lambda)



# binary data
glmeObject2_nlme <- list(
  fm= Zij ~ 1 + t + sin_term + t_dij_star + b_11 + b_21*t,
  family='binomial', par="alpha", ran.par=c('b_11','b_21'),
  sigma=NULL, 
  disp=c("alpha4", 'alpha5'),
  lower=c(-Inf,0),
  upper=rep(Inf,2),
  str_val=c(fixef(md3_same),  sd(ranef(md3_same)$sid[,1])),
  CenObject=NULL)


glmeObject_nlme <- list(nlmeObject, glmeObject2_nlme)


survObject_nlme <- list(
  fm= L ~ GNE8_CD4 + risk1 + risk2 + b_11 +b_21,
  event="dropped_out_right", 
  par='xi',
  disp=NULL,
  lower=NULL, upper=NULL,
  distribution=NULL,
  str_val= summary(fitCOX1_nlme)$coeff[,1])

set.seed(2)
H <- JMfit(
  glmeObject = glmeObject_nlme,
  survObject = survObject_nlme,
  long.data = long_df, 
  surv.data = surv_data_final,
  idVar="sid", 
  eventTime="L",
  survFit= fitCOX1_nlme,
  method = "h-likelihood", Silent=T)


JMsummary(H)
new_sd2 = JMsd_aGH(H, ghsize=4, srcpath=srcpath, paralle=T)
JMsummary(H, new_sd2$sd_robust)
