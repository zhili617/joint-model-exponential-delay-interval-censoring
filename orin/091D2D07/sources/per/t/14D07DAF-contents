#' Create interval-censored Weibull log-likelihood function
#' For use in joint model survival component

weibull_interval_loglike <- function(survObject) {
  # Extract required elements from the survObject
  model <- survObject$fm
  timegrid <- survObject$timegrid
  distribution <- survObject$distribution
  lower <- survObject$lower
  upper <- survObject$upper
  sigma <- survObject$sigma
  str_val <- survObject$str_val
  
  # Extract variable names from model formula
  model_vars <- all.vars(model)
  resp_vars <- model_vars[1:2]  # L and R
  covariates <- model_vars[-(1:2)]  # fixed effects
  
  # Build parameter names based on 'par' field
  n_cov <- length(covariates)
  par_names <- paste0(survObject$par, 0:(n_cov - 1))
  names(str_val)[1:n_cov] <- par_names
  
  # Construct linear predictor string
  linear_pred <- paste(par_names, "*", covariates, collapse = " + ")
  
  # Extract scale and shape parameter names
  par_lambda <- survObject$disp[1]  # usually "lambda"
  par_phi <- survObject$disp[2]     # usually "phi"
  
  # Determine number of intervals
  n_interval <- length(timegrid) - 1
  loglike_parts <- vector("list", n_interval)
  
  # Build log-likelihood expression
  for (j in 1:n_interval) {
    t_j1 <- paste0("(", timegrid[j], ")^", par_phi)
    t_j <- paste0("(", timegrid[j + 1], ")^", par_phi)
    delta_j <- paste0("(", t_j, " - ", t_j1, ")")
    
    haz_expr <- paste0(par_lambda, " * exp(", linear_pred, ") * ", delta_j)
    pi_ij <- paste0("1 - exp(-(", haz_expr, "))")
    
    if (j == 1) {
      inner_sum <- ""
    } else {
      inner_terms <- c()
      for (l in 1:(j - 1)) {
        t_l <- paste0("(", timegrid[l + 1], ")^", par_phi)
        t_l1 <- paste0("(", timegrid[l], ")^", par_phi)
        delta_l <- paste0("(", t_l, " - ", t_l1, ")")
        term_l <- paste0("-(", par_lambda, " * exp(", linear_pred, ") * ", delta_l, ")")
        inner_terms <- c(inner_terms, term_l)
      }
      inner_sum <- paste(inner_terms, collapse = " + ")
    }
    
    rname <- paste0("r", j)
    loglike_j <- paste0("(", rname, ") * (", inner_sum, if (inner_sum != "") " + ", "log(", pi_ij, "))")
    loglike_parts[[j]] <- loglike_j
  }
  
  loglike_expr <- paste(loglike_parts, collapse = " + ")
  
  # Display the final expression
  message("=== Interval-censored Weibull Log-Likelihood Expression ===")
  cat(loglike_expr, "\n")
  
  return(list(
    resp = resp_vars,
    loglike = loglike_expr,
    par = par_names,
    linear_pred = linear_pred,
    str_val = str_val,
    sigma = sigma,
    lower = lower,
    upper = upper
  ))
}
