library(readr)
library(tidyverse)
library(survival)
data <- read_csv("E:/UBC/Final Project/vax004dataAb.txt")
View(data)
long.data <- read_csv("example/Longdata.csv")
View(Longdata)
surv.data <- read_csv("example/Survdata.csv")

#------------survival data create--------------
infected <- data %>%
  filter(HIVinfectionind == 1, Vxind == 1)
surv.data <- infected %>%
  group_by(sid) %>%
  summarize(
    L = min(sampledays),  # the earliest observed time
    R = max(sampledays),  # latest observed time
    NAb = first(NAb),     # fixed covariates
    riskscorecat = first(riskscorecat),
    .groups = "drop"
  ) %>%
  mutate(event = 3)


S <- with(surv.data, Surv(L, R, event, type = "interval2"))






# use group_modify，keep each ID's history
surv_data <- df_vx %>%
  group_by(sid) %>%
  group_modify(~{
    df <- .x
    times <- df$sampledays
    infect <- df$HIVinfectionind
    
    if (all(infect == 0)) {
      # never infect right-censored
      return(tibble(
        left = max(times),
        right = Inf,
        event = 0
      ))
    }
    
    # right part
    first_infect_time <- min(times[infect == 1])
    
    # check if the infection during the test period
    prior_neg <- times[infect == 0 & times < first_infect_time]
    
    if (length(prior_neg) > 0) {
      # interval-censored 
      return(tibble(
        left = max(prior_neg),
        right = first_infect_time,
        event = 3
      ))
    } else {
      # no prior negative -> left-censored
      return(tibble(
        left = 0,
        right = first_infect_time,
        event = 2
      ))
    }
  }) %>%
  ungroup()

df_vx %>% filter(sid ==104716) %>% arrange(sampledays) %>% select(sid, sampledays, HIVinfectionind)



df_vx %>%
  arrange(sid, sampledays) %>%
  group_by(sid) %>%
  summarise(
    has_neg = any(HIVinfectionind == 0),
    has_pos = any(HIVinfectionind == 1),
    first_pos_time = if (any(HIVinfectionind == 1)) min(sampledays[HIVinfectionind == 1]) else NA_real_,
    has_neg_before_pos = if (any(HIVinfectionind == 1)) {
      any(HIVinfectionind == 0 & sampledays < min(sampledays[HIVinfectionind == 1]))
    } else {
      FALSE
    }
  ) %>% filter(has_neg_before_pos == TRUE)



long.data %>%
  group_by(sid) %>%
  summarise(
    both = any(HIVinfectionind == 1) &  any(HIVinfectionind == 0)
 ) %>%
  count(both)

sum(df_vx$has_neg_before_pos == TRUE)

















cutoff <- 2.0

interval_data <- data %>%
  arrange(sid, sampledays) %>%
  group_by(sid) %>%
  mutate(event = NAb > cutoff) %>%
  summarise(
    has_event = any(event),
    event_idx = which(event)[1],
    R = ifelse(has_event, sampledays[event_idx], 2000),
    L = ifelse(has_event & event_idx > 1, sampledays[event_idx - 1], 0),
    status = case_when(
      !has_event ~ 0,              # right-censored
      has_event & is.na(L) ~ 1,    # left-censored
      has_event ~ 2                # interval-censored
    )
  ) 


interval_data$new_id <- seq_len(nrow(interval_data))

unique_ids <- unique(data$sid)

# give unique sid a base value from N(0, 1)
set.seed(123)
base_values <- data.frame(
  sid = unique_ids,
  base = rnorm(length(unique_ids), mean = 0, sd = 1)
)

# add base into long and survival
long <- merge(data, base_values, by = "sid")
survival <- merge(interval_data, base_values, by = "sid")


#-----------------------------

long_data <- data %>%
  mutate(
    t_ij = sampledays / 30,           # month
    t_star_ij = t_ij / 12             # year
  )


long_data <- long_data %>%
  group_by(sid) %>%
  mutate(
    last_vaccine_day = max(sampledays[dosenumberschedule == max(dosenumberschedule)]), # 最后接种时间
    d_ij = sampledays - last_vaccine_day
  ) %>% ungroup()


long_data <- long_data %>%
  mutate(
    risk1 = as.numeric(riskscorecat == 1),
    risk2 = as.numeric(riskscorecat == 2)
  )


long_data <- long_data %>%
  mutate(
    GNE8_V3_bin = as.numeric(GNE8_V3 > log10(100))  # 你已有代码
  )


# -----------------------------------longitudinal models (no change)---------------------------
LLOQ <- 1.477  # log10(30) from Gilbert et al. (2005)

long_data <- long %>%
  mutate(year = sampledays / 365.25, year2 = (sampledays / 365.25)^2)
long_data$GNE8_V3_bin <- as.numeric(long_data$GNE8_V3 > log10(100))


md_NAb <- lmer(NAb ~ t_star_ij + I(t_star_ij^2) + sin(pi / Delta_ij * d_ij) + risk1 + risk2 + (1 | sid), data = long_data)
md_bin <- glmer(GNE8_V3_bin ~ t_ij + sin(pi / Delta_ij * d_ij) + t_star_ij + risk1 + risk2 + (1 + t_ij | sid),
                data = long_data, family = binomial)


ranef_NAb <- ranef(md_NAb)$sid
ranef_NAb$sid <- rownames(ranef_NAb)
colnames(ranef_NAb)[1] <- "b11"

ranef_bin <- ranef(md_bin)$sid
ranef_bin$sid <- rownames(ranef_bin)
colnames(ranef_bin)[1:2] <- c("b11_bin", "b21_bin")


long_data <- merge(long_data, ranef_NAb[, c("sid", "b11")], by = "sid", all.x = TRUE)
long_data <- merge(long_data, ranef_bin[, c("sid", "b11_bin", "b21_bin")], by = "sid", all.x = TRUE)





# longitutional part


glmeObject1 <- list(
  fm = NAb ~ t_star_ij + I(t_star_ij^2) + sin(pi / Delta_ij * d_ij) + risk1 + risk2 + b11,
  family = "normal",
  par = c("Intercept", "t_star_ij", "I_t_star_ij^2", "sin_term", "risk1", "risk2", "b11"),
  ran.par = "b11",
  sigma = "sigma",
  disp = "beta4",
  lower = rep(-Inf, 7),
  upper = rep(Inf, 7),
  str_val = c(fixef(md_NAb), sd(ranef(md_NAb)$sid[,1])),
  CenObject = NULL
)





glmeObject2 <- list(
  fm = GNE8_V3_bin ~ t_ij + sin(pi / Delta_ij * d_ij) + t_star_ij + risk1 + risk2 + b11_bin + b21_bin,
  family = "binomial",
  par = c("Intercept", "t_ij", "sin_term", "t_star_ij", "risk1", "risk2", "b11_bin", "b21_bin"),
  ran.par = c("b11_bin", "b21_bin"),
  disp = "alpha5",
  lower = rep(-Inf, 8),
  upper = rep(Inf, 8),
  str_val = c(fixef(md_bin), sqrt(diag(as.matrix(VarCorr(md_bin)$sid)))),
  CenObject = NULL
)






glmeObject <- list(glmeObject1, glmeObject2)



#--------------------------------------------------------

timegrid <- sort(unique(long$doesnumberschedule))  
timegrid <- timegrid[!is.na(timegrid)]             
K <- length(timegrid) - 1                           
rnames <- paste0("r", 1:K)  


survObject_weibull_interval <- list(
  resp = rnames,
  linear_pred = "base + b11",   
  timegrid = timegrid,
  rnames = rnames,
  par = c("lambda", "phi", "gamma0", "gamma1", "gamma2", "gamma3", "gamma4"),
  sigma = c("lambda", "phi"),
  lower = c(0, 0, -Inf, -Inf, -Inf, -Inf, -Inf),
  upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf),
  str_val = c(0.01, 1.2, 0, 0, 0, 0, 0),
  distribution = "weibull_interval",
  disp = NULL
)




fit_weibull <- survreg(Surv(L, R, type = "interval2") ~ base + risk1 + risk2, data = survival, dist = "weibull")


fit <- JMfit(
  glmeObject = list(glmeObject1, glmeObject2),
  survObject = survObject_weibull_interval,
  long.data = long_data,
  surv.data = survival,
  idVar = "sid",
  survFit = fit_weibull,
  method = "h-likelihood"
)






