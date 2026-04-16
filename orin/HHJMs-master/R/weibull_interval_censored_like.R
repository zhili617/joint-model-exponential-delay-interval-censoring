#' Weibull PH interval-censored log-likelihood (for survival component)
weibull_ph_interval_loglike <- function(survObject, weibPar){
  
  model <- survObject$fm # equation : Surv(L, R) ~ x1 + x2 + ...
  resp  <- fmReverse(model)$resp # Surv(L,R)
  rvX   <- fmReverse(model)$rvX # name of independent variables
  event <- survObject$event
  
  
  # names for L, R
  model_vars <- all.vars(model)
  L_name <- model_vars[1]
  R_name <- model_vars[2]
  
  # beta names + linear_pred 
  if(sum(toupper(rvX)=='NULL')>0 | sum(rvX=='1')>0){ # ensure there are covariance
    linear1 <- "0"
    par <- c()
  } else {
    p <- length(rvX) # number of covariance
    par <- paste(survObject$par, 1:p-1, sep="") # name of parameters eg. lambda0, lambda1, ...
    linear1 <- paste(par,"*", rvX, collapse="+", sep="") # if par = c("lambda0","lambda1") and rvX=c("risk1","risk2")，we have: "lambda0*risk1+lambda1*risk2"
  }
  linear_pred <- linear1
  
  # Weibull PH params
  if(is.null(weibPar)) weibPar <- c(0, 1)  # initial Wlogscale, Wshape, if there is no weibPar
  Wlogscale <- "Wlogscale"
  Wshape    <- "Wshape"
  
  SL <- paste0("exp( -exp(", Wlogscale, "+(", linear_pred, ")) * (", L_name, ")^", Wshape, " )")
  SR <- paste0("exp( -exp(", Wlogscale, "+(", linear_pred, ")) * (", R_name, ")^", Wshape, " )")
  
  
  
  # eps <- 1e-16 # incase S(L)-S(R) = 0
  # loglike <- paste0(
  #  "  log( (", SL, " - ", SR, ") )  "
  #  )
  
  
  eps <- 1e-16
  
  logSL <- paste0(
    "-exp(", Wlogscale, "+(", linear_pred, ")) * (", L_name, ")^", Wshape
  )
  logSR <- paste0(
    "-exp(", Wlogscale, "+(", linear_pred, ")) * (", R_name, ")^", Wshape
  )
  
  
  # The case with only interval censored
  loglike_int <- paste0(
    "(", logSL, ") + log( 1 - exp( (", logSR, ") - (", logSL, ") ) + ", eps, " )"
  )
  
  # The case with only right censored
  loglike_right <- logSL
  
  # The case with both interval-censored and right-censored
    loglike <- paste0(
      "(", event, ")*(", loglike_int, ") + ",
      "(1-(", event, "))*(", loglike_right, ")")
  
  
  # parameters
  str_val <- c(survObject$str_val, weibPar)
  par <- c(par, Wlogscale, Wshape)
  names(str_val) <- par #Notice we need length(survObject$str_val) == length(par) 
  
  return(list(resp=c(L_name, R_name), loglike=loglike, par=par,
              linear_pred=linear_pred, str_val=str_val))
}

