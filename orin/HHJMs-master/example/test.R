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

#----------------------
# ---- save result ----
#----------------------

out_dir <- "E:/UBC/Final Project/New folder/Final-Project/orin/HHJMs-master/example"

write_csv(
  long_data,
  file.path(out_dir, paste0("long_data", ".csv"))
)

write_csv(
  surv_data_final,
  file.path(out_dir, paste0("surv_data", ".csv"))
)


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

set.seed(1)
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

# random effect
testjm1$sigma

testjm1$covBi
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
testjm2 <- JMfit(glmeObject, survObject2, 
                     long_data, surv_data_final,
                     idVar="sid", eventTime="L",
                     survFit=fitCOX2,
                     method = "h-likelihood", Silent = T)

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


#-------------(1) Create objects for survival models---------------

# rename the coefficients
co <- coef(surv_model_ph)

phi <- exp(co["shape"])
lambda_base <- exp(co["scale"])
lambda_slopes <- co[c("GNE8_CD4", "risk1", "risk2", "b11", "b21")]
str_val_surv <- c(unname(lambda_slopes), lambda_base, phi)

names(str_val_surv) <- c(
  "lambda0", "lambda1", "lambda2", "lambda3", "lambda4",
  "lambda_base", "phi"
)





survObject <- list(
  fm   = Surv(L, R, type = "interval2") ~ GNE8_CD4 + risk1 + risk2 + b11 + b21,
  event = "interval",
  par   = "lambda",                 
  disp  = c("lambda_base", "phi"),  
  lower = c(0, 0),                  
  upper = rep(Inf, 2),
  distribution = "weibull_interval",
  str_val = str_val_surv
)


testjm_interval <- try(JMfit(
  glmeObject = glmeObject,
  survObject = survObject,   
  long.data  = long_data,
  surv.data  = Sdata,
  idVar      = "sid",
  eventTime  = NULL,        
  survFit    = surv_model_ph,         
  method     = "h-likelihood",
  Silent     = FALSE
), silent = TRUE)





#--------------simulation-----------------------

1111
#-----------------------------------
#--------------Plot-----------------
#-----------------------------------



# Oringal GNE8_CD4
ggplot(long_data, aes(x = t, y =  sin_term, group = sid)) +
  geom_line(color = "black", alpha = 0.5) +   
  geom_point(color = "black", size = 1) + 
  labs(title = "GNE8_CD4 longitudinal trajectories", y = "GNE8_CD4", x = "month") +
  theme_minimal()

ggplot(long_data_simulation, aes(x = t, y =  t_dij_star, group = sid)) +
  geom_line(color = "black", alpha = 0.5) +   
  geom_point(color = "black", size = 1) + 
  labs(title = "GNE8_CD4 longitudinal trajectories", y = "GNE8_CD4", x = "month") +
  theme_minimal()

# randomly pick 5 sid,and draw their line plot (GEN8_CD4)
set.seed(2020) # GEN -> 2025, NAb -> 2020
subset_ids <- sample(unique(long_data$sid), 5)

long_data %>% 
  filter(sid %in% subset_ids) %>%
  ggplot(aes(x = t, y = GNE8_CD4, group = sid, color = factor(sid))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = paste("GNE8_CD4 trajectories for 5 random subjects"),
    x = "Month", y = "GNE8_CD4", color = "Subject ID"
  ) +
  theme_minimal()


# Original NAb
ggplot(long_data, aes(x = t, y =  NAb, group = sid)) +
  geom_line(color = "black", alpha = 0.5) +   
  geom_point(color = "black", size = 1) + 
  labs(title = "NAb longitudinal trajectories", y = "NAb", x = "month") +
  theme_minimal()


# randomly pick 5 sid,and draw their line plot (NAb)

long_data %>% 
  filter(sid %in% subset_ids) %>%
  ggplot(aes(x = t, y = NAb, group = sid, color = factor(sid))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = paste("NAb trajectories for 5 random subjects"),
    x = "Month", y = "NAb", color = "Subject ID"
  ) +
  theme_minimal()






# Simulated GEN8_CD4
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

