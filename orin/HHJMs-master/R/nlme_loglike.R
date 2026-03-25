nlme_loglike_EXP_Delay <- function(nlmeObject){
  
  vars <- fmReverse(nlmeObject$fm)
  resp <- vars$resp
  
  A_names    <- nlmeObject$A_names # c("A1","A2",...,"A7")
  ind_names  <- nlmeObject$ind_names # c("ind_k1",...,"ind_k7")
  lag_names  <- nlmeObject$lag_names # c("lag_k1",...,"lag_k7")
  lambda_nm  <- nlmeObject$lambda_name
  ran_map    <- nlmeObject$ran_map # list which parameter has random effect and the name of random effect, eg. list(lambda="b_lambda", A2="b_A2")  
  # If there is no random effect, write list()
  
  # If any of the parameter has random effect, then return that parameter + random effect
  add_re <- function(par_nm){
    if (!is.null(ran_map[[par_nm]])) {
      paste0("(", par_nm, "+", ran_map[[par_nm]], ")")
    } else {
      par_nm
    } }
  
  lambda_i <- paste0("exp(", add_re(lambda_nm), ")")
  
  A_i <- vapply(A_names, add_re, FUN.VALUE = character(1))
  
  # mu = sum_{k=1}^7 A_k_i * ind_k * exp( -lambda_i * lag_k )
  terms <- character(length(A_names))
  for(k in seq_along(A_names)){
    terms[k] <- paste0(A_i[k], "*", ind_names[k],
                       "*exp(-", lambda_i, "*", lag_names[k], ")")
  }
  mu <- paste(terms, collapse = "+")
  
  # error term ~ normal (NAb | b ~ N(mu, sigma^2))
  sigma <- nlmeObject$sigma
  loglike <- paste0(
    "-0.5*(", resp, "-(", mu, "))^2/", sigma, "^2",
    "-log(", sigma, ")-0.5*log(2*pi)"
  )
  
  fixed <- c(A_names, lambda_nm)
  raneff <- if (length(ran_map) > 0) unlist(ran_map, use.names = FALSE) else NULL
  rvX <- c(ind_names, lag_names)      # NLME: needed data columns
  rvZ <- names(ran_map)               # NLME: which parameters have RE diff with them in glme
  nonlinear_pred <- mu                   # store mu here for downstream reuse
  
  list(
    loglike = loglike,
    fixed = fixed,
    raneff = raneff,
    mu = mu,
    linear_pred = nonlinear_pred,
    rvX = rvX, rvZ = rvZ,
    resp = resp, sigma = sigma
  )
}