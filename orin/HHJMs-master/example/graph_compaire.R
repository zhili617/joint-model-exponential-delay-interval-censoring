library(survival)
library(survivalROC)
library(dplyr)
library(pROC)

graph_data <- dataset_new_NLME

graph_data <- graph_data %>%
  mutate(sid = as.character(sid))

re_df <- data.frame(
  sid = as.character(testjm1$uniqueID),
  b11 = testjm1$Bi[, 1],
  b21 = testjm1$Bi[, 2]
)

# re_df <- data.frame(
#   sid = as.character(testjm_without$uniqueID),
#   b11 = testjm_without$Bi[, 1],
#   b21 = testjm_without$Bi[, 2]
# )
graph_data_T <- graph_data %>%
  left_join(
    re_df %>% dplyr::select(sid, b11),
    by = "sid"
  )

re_H <- data.frame(
  sid = as.character(H$uniqueID),
  b11_nlme = H$Bi[, 1],
  b21_nlme = H$Bi[, 2]
)

long_new <- data.frame(
  sid = as.character(long_data$sid),
  risk1 = long_data$risk1,
  risk2 = long_data$risk2,
  sin_term = long_data$sin_term,
  t_star = long_data$t_star,
  tdij = long_data$tdij,
  delta_ij = long_data$delta_ij,
  t = long_data$t
)

graph_data_A <- graph_data_T  %>%
  left_join(
    long_new %>% dplyr::select(sid, sin_term, risk1, t_star, risk2, tdij, delta_ij, t),
    by = c("sid", "t")
  )


graph_data_N <- graph_data_A %>%
  left_join(
    re_H %>% dplyr::select(sid, b11_nlme),
    by = "sid"
  )

graph_data_N <- graph_data_N %>%
  distinct(sid, sampledays, .keep_all = TRUE)


graph_data_N <- graph_data_N %>%
  mutate(actual_dose_day = sampledays - dosedays)


# ----------------------------------------------------
# f_model1_smooth <- function(t, A, lambda,
#                             plan_months = c(0,1,6,12,18,24,30),
#                             kappa = 20) {
#   
#   smooth_step <- function(x) {
#     z <- kappa * x
#     ifelse(z >= 0,
#            1 / (1 + exp(-z)),
#            exp(z) / (1 + exp(z)))
#   }
#   
#   softplus <- function(x) {
#     z <- kappa * x
#     ifelse(z > 30,
#            x,
#            log1p(exp(z)) / kappa)
#   }
#   
#   val <- rep(0, length(t))
#   for (k in seq_along(plan_months)) {
#     x <- t - plan_months[k]
#     val <- val + A[k] * smooth_step(x) *
#       exp(-exp(lambda) * softplus(x))
#   }
#   
#   val
# }
# ----------------------------------------------------






get_delta_grid <- function(tgrid, vac_times = vaccine_time_map) {
  last_vac_time <- tail(vac_times, 1)
  delta <- numeric(length(tgrid))
  
  for (i in seq_along(tgrid)) {
    ti <- tgrid[i]
    assigned <- FALSE
    
    for (j in 1:(length(vac_times) - 1)) {
      start <- vac_times[j]
      end   <- vac_times[j + 1]
      
      if (ti >= start && ti < end) {
        delta[i] <- end - start
        assigned <- TRUE
        break
      }
    }
    
    if (!assigned) {
      if (ti >= last_vac_time) {
        delta[i] <- ti - last_vac_time
      } else {
        delta[i] <- 0
      }
    }
  }
  
  delta
}






set.seed(20) #2
sids_sample <- sample(unique(graph_data_N$sid), 1)   

# # without LLOQ section model
# beta_vec <- c(
#   beta0 = as.numeric(testjm_without$fixedest["beta0"]),
#   beta1 = as.numeric(testjm_without$fixedest["beta1"]) ,
#   beta2 = as.numeric(testjm_without$fixedest["beta2"]),
#   beta3 = as.numeric(testjm_without$fixedest["beta3"]),
#   beta4 = as.numeric(testjm_without$fixedest["beta4"]),
#   unlist(testjm_without$sigma)["beta5"]
# )

# # with LLOQ section model
# beta_vec <- c(
#   beta0 = as.numeric(testjm1$fixedest["beta0"]),
#   beta1 = as.numeric(testjm1$fixedest["beta1"]) ,
#   beta2 = as.numeric(testjm1$fixedest["beta2"]),
#   beta3 = as.numeric(testjm1$fixedest["beta3"]),
#   beta4 = as.numeric(testjm1$fixedest["beta4"]),
#    unlist(testjm1$sigma)["beta5"]
# )

# nonlinear model function
f_model1 <- function(t, A, lambda, plan_months = c(0,1,6,12,18,24,30)) {
  val <- rep(0, length(t))
  for (k in 1:7) {
    val <- val + A[k] * as.integer(t >= plan_months[k]) *
      exp(-exp(lambda) * pmax(t - plan_months[k], 0))
  }
  val
}

par(mfrow = c(length(sids_sample), 1), mar = c(4,4,3,1))

for (sid_i in sids_sample) {
  
  dat_i <- graph_data_N %>% filter(sid == sid_i)
  
  tgrid <- seq(min(dat_i$t), max(dat_i$t), length.out = 500)
  
  # observed points
  plot(dat_i$t, dat_i$NAb,
       pch = 16, col = "black",
       main = paste("SID =", sid_i),
       xlab = "Months", ylab = "NAb", ylim = c(0, 8))
  
  ## ======================================================
  ## 1) H / fit_nlme subject-specific nonlinear curve
  ## ======================================================
  
  b11_nlme <- dat_i$b11_nlme[1]

  # Get lambda and c_b11 from converged model
  # lambda is in fixedest (log-scale decay rate)
  # c_b11 is in sigma (dispersion) if nlmeObject$disp = "c_b11"
  lambda_fixed <- H$fixedest["lambda"]

  coefs_sid <- c(
    A1 = as.numeric(H$fixedest["A1"]),
    A2 = as.numeric(H$fixedest["A2"]),
    A3 = as.numeric(H$fixedest["A3"]),
    A4 = as.numeric(H$fixedest["A4"]),
    A5 = as.numeric(H$fixedest["A5"]),
    A6 = as.numeric(H$fixedest["A6"]),
    A7 = as.numeric(H$fixedest["A7"]),
    lambda = as.numeric(lambda_fixed) +  b11_nlme
  )
  
  A_names <- paste0("A", 1:7)
  A_vals  <- as.numeric(coefs_sid[A_names])
  lambda_val <- as.numeric(coefs_sid["lambda"])   # subject-specific log-rate
  
  curve_H <- f_model1(tgrid, A = A_vals, lambda = lambda_val)
  lines(tgrid, curve_H, col = "blue", lwd = 2)
  
  ## optional: nlme fitted points at observed times
  pred_points <- f_model1(dat_i$t, A = A_vals, lambda = lambda_val)
  points(dat_i$t, pred_points, col = "blue", pch = 16)
  
  ## ======================================================
  ## 2) testjm1 subject-specific linear predictor curve
  ## NAb ~ t_star + sin_term + risk1 + risk2 + b11
  ## ======================================================
  risk1_i <- dat_i$risk1[1]
  risk2_i <- dat_i$risk2[1]
  b11_i   <- dat_i$b11[1]
  
  t_star_grid   <- tgrid / 12
  
  actual_imm <- sort(unique(round((dat_i$sampledays - dat_i$dosedays) / 30, 3)))
  
  tdij_grid <- sapply(tgrid, function(tt) {
    prev_candidates <- actual_imm[actual_imm <= tt]
    if (length(prev_candidates) == 0) {
      return(NA_real_)
    } else {
      return(tt - max(prev_candidates))
    }
  })
  
  delta_grid <- get_delta_grid(tgrid, vaccine_time_map)
  
  sin_term_grid <- ifelse(delta_grid != 0,
                          sin(pi * tdij_grid / delta_grid),
                          sin(pi * tdij_grid / 0.01))
  
  curve_beta <- beta_vec["beta0"] +
    beta_vec["beta1"] * (tgrid / 12) +
    beta_vec["beta2"] * sin_term_grid +
    beta_vec["beta3"] * risk1_i +
    beta_vec["beta4"] * risk2_i +
    beta_vec["beta5"] * b11_i
  
  lines(tgrid, curve_beta, col = "red", lwd = 2, lty = 2)
  
  pred_lme_points <- beta_vec["beta0"] +
    beta_vec["beta1"] * dat_i$t_star +
    beta_vec["beta2"] * dat_i$sin_term +
    beta_vec["beta3"] * dat_i$risk1 +
    beta_vec["beta4"] * dat_i$risk2 +
    beta_vec["beta5"] * dat_i$b11
  
  points(dat_i$t, pred_lme_points, col = "red", pch = 16)
  
  dose_time_map <- c(`1` = 0, `2` = 1, `3` = 6, `4` = 12, `5` = 18, `6` = 24, `7` = 30)
  
  dose_nums <- sort(unique(dat_i$dosenumberschedule))
  dose_times <- dose_time_map[as.character(dose_nums)]
  
  abline(v = dose_times, col = "darkgreen", lty = 3)
  
  
  text(
    x = dose_times,
    y = 0,   
    labels = paste0("dose", dose_nums),
    col = "darkgreen",
    cex = 0.8,
    pos = 3  
  )
  
  legend("topright",
         inset = c(-0.2, -0.05),
         legend = c("Observed NAb",
                    "NLME curve",
                    "NLME fitted points",
                    "LME curve",
                    "LME fitted points",
                    "Immunization"),
         col = c("black", "blue", "blue", "red", "red", "darkgreen"),
         pch = c(16, NA, 16, NA,16, NA),
         lty = c(NA, 1, NA, 2, NA, 3),
         bty = "n")
  
  # check
  tdij_obs_rebuild <- sapply(dat_i$t, function(tt) {
    prev_candidates <- actual_imm[actual_imm <= tt]
    tt - max(prev_candidates)
  })
  
  delta_obs_rebuild <- get_delta_grid(dat_i$t, vaccine_time_map)
  
  sin_obs_rebuild <- ifelse(delta_obs_rebuild != 0,
                            sin(pi * tdij_obs_rebuild / delta_obs_rebuild),
                            sin(pi * tdij_obs_rebuild / 0.01))
  
  check_compare_new <- data.frame(
    t = dat_i$t,
    tdij_data = dat_i$tdij,
    tdij_rebuild = tdij_obs_rebuild,
    delta_data = dat_i$delta_ij,
    delta_rebuild = delta_obs_rebuild,
    sin_data = dat_i$sin_term,
    sin_rebuild = sin_obs_rebuild
  )
  
  print(check_compare_new)
}




# ----------------------------------------------------------------


# ------------------- Obseved v.s. fitted NAb (MSE and SSE)-------------



## -----------------------------
##  fitted points
## -----------------------------
all_sid <- unique(graph_data_N$sid)

pred_list <- vector("list", length(all_sid))

for (ii in seq_along(all_sid)) {
  
  sid_com <- all_sid[ii]
  dat_com <- graph_data_N %>% filter(sid == sid_com)
  
  ## ===== NLME fitted points =====
  b11_nlme <- dat_com$b11_nlme[1]
  
   lambda_fixed <- H$fixedest["lambda"]

  coefs_sid <- c(
    A1 = as.numeric(H$fixedest["A1"]),
    A2 = as.numeric(H$fixedest["A2"]),
    A3 = as.numeric(H$fixedest["A3"]),
    A4 = as.numeric(H$fixedest["A4"]),
    A5 = as.numeric(H$fixedest["A5"]),
    A6 = as.numeric(H$fixedest["A6"]),
    A7 = as.numeric(H$fixedest["A7"]),
    lambda = as.numeric(lambda_fixed) +  b11_nlme
  )
  A_vals <- as.numeric(coefs_sid[paste0("A", 1:7)])
  lambda_val <- as.numeric(coefs_sid["lambda"])
  
  pred_nlme_points <- f_model1(dat_com$t, A = A_vals, lambda = lambda_val)
  
  ## ===== LME fitted points =====
  pred_lme_points <- beta_vec["beta0"] +
    beta_vec["beta1"] * dat_com$t_star +
    beta_vec["beta2"] * dat_com$sin_term +
    beta_vec["beta3"] * dat_com$risk1 +
    beta_vec["beta4"] * dat_com$risk2 +
    beta_vec["beta5"] * dat_com$b11
  
  ## ===== save =====
  pred_list[[ii]] <- dat_com %>%
    mutate(
      pred_nlme = pred_nlme_points,
      pred_lme  = pred_lme_points
    ) %>%
    dplyr::select(sid, t, NAb, pred_nlme, pred_lme)
}

pred_df <- bind_rows(pred_list)


# --------- Compute SSE, MSE -------------------

SSE_nlme <- sum((pred_df$NAb - pred_df$pred_nlme)^2, na.rm = TRUE)
SSE_lme  <- sum((pred_df$NAb - pred_df$pred_lme)^2, na.rm = TRUE)

MSE_nlme <- mean((pred_df$NAb - pred_df$pred_nlme)^2, na.rm = TRUE)
MSE_lme  <- mean((pred_df$NAb - pred_df$pred_lme)^2, na.rm = TRUE)

data.frame(
  Model = c("NLME", "LME"),
  SSE = c(SSE_nlme, SSE_lme),
  MSE = c(MSE_nlme, MSE_lme)
)