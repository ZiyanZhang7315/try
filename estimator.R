###estimator.R

## function estimator accpet scale n, monte carlo simulation number Nmc and variance component s2_gxg, s2_e as input
## output the approximated trW, trW2, yTWy, exact estimator s2_gxg_exact, s2_e_exact approximation estimator s2_gxg_app, s2_gxg_exact
estimator <- function(n, m, Nmc, s2_gxg, s2_e){
    ## load data from data_simulation
    simulation_data <-data_simulation(n, m, s2_gxg, s2_e)
    X     <- simulation_data$X
    y     <- simulation_data$y
    trW   <- simulation_data$trW
    trW2  <- simulation_data$trW2
    yTWy  <- simulation_data$yTWy

    rowsumX2 = rowSums(X^2)
    p = m*(m-1)/2
    ##trw approximation
    res_mat <- replicate(Nmc, {
      u <- rnorm(m)           
      v <- X %*% u           
      (v^2 - rowsumX2)^2     
    })
    est_trW <- sum(res_mat) / (4*p*Nmc)
    ##trw2 approximation
    ## calculate Z_hat
    Zmat <- replicate(Nmc, {
    u1 <- rnorm(m)
    v  <- X %*% u1
    v^2 - rowsumX2
     })
    
    z_hat <- rowSums(Zmat) / (2 * sqrt(Nmc))
    element1 <- sum(replicate(Nmc, {u2 <- rnorm(n); 
                               v1 <- X %*% (t(X) %*% (z_hat * u2)); 
                               v2 <- z_hat * (X %*% (t(X) %*% u2));
                               t(v2) %*% v1 }))/(2*p^2*Nmc)

    est_trW2 = c(element1 - (sum((t((X*X))%*%z_hat)^2))/(2*p^2))

    ##ytwy approximation
    ## esitmate y^T K circ K y 
    est_ytkky <- sum(replicate(Nmc, {u3 <- rnorm(n); 
                             v1 <- X %*% (t(X) %*% (y * u3)); 
                             v2 <- y * (X %*% (t(X) %*% u3));
                             t(v2) %*% v1 }))/Nmc

    est_yTWy<- (est_ytkky - t(y)%*% ((X^2) %*% t(X^2))%*%y)/(2*p) 

    ## out put the exact estimator and approximated estimator
    exact_LHS <- matrix(c(simulation_data$trW2, simulation_data$trW, simulation_data$trW, n), 2, 2)
    exact_RHS <- c(yTWy, t(y)%*% y)
    var_exact <- solve(exact_LHS, exact_RHS) 

    app_LHS <- matrix(c(est_trW2, est_trW, est_trW, n), 2, 2)
    app_RHS <- c(est_yTWy, t(y)%*% y)
    var_app <- solve(app_LHS, app_RHS) 

    ## the output
    c(var_exact[1], var_exact[2], var_app[1], var_app[2],
    trW, est_trW, trW2, est_trW2,
    yTWy, est_yTWy)
}