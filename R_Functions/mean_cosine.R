library(gdata)

mean_cosine <- function(matrix, by.row = FALSE){
  C <- matrix
  if(by.row==TRUE) C <- t(C)
  M <- t(C) %*% C
  T <- diag(M)
  Q <- sqrt(outer(T, T))
  costheta <- as.numeric(upperTriangle(M / Q))
  return(mean(costheta))
}