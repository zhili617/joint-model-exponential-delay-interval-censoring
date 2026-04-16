# evaluate a matrix 

evalMat <- function(mat, q, data=NULL, par.val=NULL, raneff.val=NULL){
  D <- matrix(0, q, q)
  
  if(!is.data.frame(par.val)){
    par.val <- data.frame(as.list(par.val))
  }
  
  for(i in 1:q){
    for(j in 1:q){
      kk <- (i-1)*q+j
      Di <- with(data, with(par.val, with(raneff.val,
                                          eval(parse(text=mat[kk])))))
      
      # If the expression is a constant (no data variable), with(data,...) returns
      # a scalar instead of a length-N vector, so we must multiply by nrow(data).
      n_obs <- if (!is.null(data)) nrow(data) else 1
      if (length(Di) == 1 && n_obs > 1) Di <- Di * n_obs
      
      D[i,j] <- sum(Di)
    }
  }
  
  return(D)
}