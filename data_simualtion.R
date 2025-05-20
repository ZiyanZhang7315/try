###data_simualtion.R

################
##load packaghe
################

if (!require(MASS)) {
  install.packages("MASS")
  library(MASS)
}

################
##define function
################

## function tr accepts matrix as input
## output the trace of the matrix
tr <- function(matrix) { 
    sum(diag(matrix)) 
}

### function simulate_genotype accepts scale n, m
### output the genotype matrix X
#simulate_genotype<- function(n, m){
#  X = matrix(rnorm(n*m),nrow = n,ncol =m)
#  list(X=X)
#}


## function data_simulation accpet scale n, m, genotype matrix X and variance component s2_gxg, s2_e as input
## output the exact value of trW, trW2 and yTWy and simulated phenotype.
data_simulation<- function(n, m, s2_gxg, s2_e){
 # generate W and y
 X = matrix(rnorm(n*m),nrow = n,ncol =m) 
 XT = t(X)
 p = m*(m-1)/2
 D = (X^2) %*% t(X^2)
 K = X %*% t(X)
 Q = K^2%*%D
 W = (K*K - D)/(2*p)
 y <- mvrnorm(mu = rep(0, n), Sigma = (1/(2*p))*(K*K - D)*s2_gxg  + s2_e * diag(n))
 Dy <- diag(c(y))
    
 # generate exact value of trW trW2 and yTWy
 trW = tr(W)
 trW2 =(tr(K^2 %*% K^2) - 2*tr(Q) + tr(D%*%D))/(4*p^2)
 yTWy = t(y) %*% W %*%y
    
 c(list(X=X),list(y=y), trW = trW, trW2 = trW2, yTWy = yTWy)
}