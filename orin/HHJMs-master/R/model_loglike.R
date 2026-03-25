# loglike for longitudinal part 

model_loglike <- function(object, std = TRUE){
  if (!is.null(object$type) && object$type == "nlme") {
    return(nlme_loglike_EXP_Delay(object))
  } else {
    return(glme_loglike(object, std = std))
  }
}