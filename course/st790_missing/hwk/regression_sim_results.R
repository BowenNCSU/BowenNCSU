#  ST 790, Homework 1

#  Simulate linear regression under various forms
#  of missingness as in EXAMPLE 2 in Chapter 1 of
#  the notes

#  Specify simulation set-up and scenario 

S <- 1000    #  Number of simulated data sets
N <- 200     #  Sample size
set.seed(4)  #  Set the seed

#  The parameter vector psi specifies the probability
#  of a complete case.  See the data generation function below 
#  for the interpretation of elements of psi

 #  (i) Missingness depends only on x (MNAR or MAR)

psi <- c(2,-0.25,0.5,0,0,0)

#  (ii) Missingness depends on y (MNAR or MAR)

#psi <- c(6,0,0,-0.075,-0.003,0.05)

#  (iii) MCAR -- missingness does not depend on y or x

#psi <- c(0.5,0,0,0,0,0) 

#  Data generation function -- outcome y and covariates
#  x and missing data indicator r corresponding to 
#  the variable that may be missing (y or x)

generate <- function(N,psi){

#   Parameter values

    px <- 0.4
    beta <- c(20,5,-5)   #  true value of beta
   
#   Generate covariates and outcome
    
    x1 <- rnorm(N,10,3)       #  continuous 
    x2 <- rbinom(N,1,px)     #  binary 
    x <- cbind(x1,x2)
    y <- beta[1]+beta[2]*x1+beta[3]*x2+rnorm(N,0,8)   #  outcome

#   Generate indicators of complete cases
    
    elinpred <- exp(psi[1]+psi[2]*x1+psi[3]*x2+psi[4]*y+psi[5]*y*x1+psi[6]*y*x2)       
    pifunc <-  elinpred/(elinpred+1)      

    r <- rbinom(N,1,pifunc)

    thedata <- list(y=y,x=x,r=r,pifunc=pifunc)
    return(thedata)
}

#  Simulation loop

results <- NULL
se.results <- NULL

for (s in 1:S){

#  Generate a data set

    thedat <- generate(N,psi)
    y <- thedat$y
    x <- thedat$x
    r <- thedat$r
    pifunc <- thedat$pifunc

#   Ideal full data regression -- this is of 
#   course not possible in practice but serves 
#   as a "gold standard"     

    full.fit <- lm(y ~ x[,1] + x[,2])
    beta.full <- full.fit$coefficients
    beta.full.se <- coef(summary(full.fit))[,2]
    
#   Complete case regression

    Nobs <- sum(r)
    yr <- y[r==1]
    xr <- x[r==1,]
    cc.fit <- lm(yr ~ xr[,1] + xr[,2])
    beta.cc <- cc.fit$coefficients
    beta.cc.se <- coef(summary(cc.fit))[,2]
    
#   Proportion of N that are complete cases 

    pobs <- Nobs/N
    
#   Save results

    results <- rbind(results,c(beta.full,beta.cc,pobs))
    se.results <- rbind(se.results,c(beta.full.se,beta.cc.se))

}


#   End of simulation loop
    
#   Summarize simulation results -- can look at
#   means, sds of estimates to see the extent of
#   bias, efficiency loss, and average proportion
#   of complete cases.  Also look at average of ses
#   to see how well these approximate the true sds

means <- apply(results,2,mean)
sds <- apply(results,2,sd)
ses <- apply(se.results,2,mean)
bias.full <- means[1:3]-c(20,5,-5)
bias.cc <- means[4:6]-c(20,5,-5)

#  Results case (i) MNAR/MAR -- CC estimator seems consistent but
#  somewhat inefficient

## > means
## (Intercept)      x[, 1]      x[, 2] (Intercept)     xr[, 1] 
##   19.970541    5.002845   -4.989779   20.060993    5.002865 
##     xr[, 2]             
##   -5.034305    0.435130 
## > sds
## (Intercept)      x[, 1]      x[, 2] (Intercept)     xr[, 1] 
##  2.07850830  0.19828353  1.18444484  3.13069797  0.33350663 
##     xr[, 2]             
##  1.80801160  0.03491493
## > ses
## (Intercept)      x[, 1]      x[, 2] (Intercept)     xr[, 1] 
##   2.0352825   0.1898984   1.1592384   2.9697481   0.3103363 
##     xr[, 2] 
##   1.7480256 
## > bias.full
##  (Intercept)       x[, 1]       x[, 2] 
## -0.029459219  0.002845078  0.010220794 
## > bias.cc
##  (Intercept)      xr[, 1]      xr[, 2] 
##  0.060992967  0.002865292 -0.034305295 

#  Results case (ii) MNAR/MAR -- CC estimator is inconsistent

## > means
## (Intercept)      x[, 1]      x[, 2] (Intercept)     xr[, 1] 
##   19.970541    5.002845   -4.989779   19.745555    4.691068 
##     xr[, 2]             
##   -2.274049    0.524590 
## > sds
## (Intercept)      x[, 1]      x[, 2] (Intercept)     xr[, 1] 
##  2.07850830  0.19828353  1.18444484  2.63939584  0.31102082 
##     xr[, 2]             
##  1.71566040  0.03457857 
## > ses
## (Intercept)      x[, 1]      x[, 2] (Intercept)     xr[, 1] 
##   2.0352825   0.1898984   1.1592384   2.6040702   0.3015416 
##     xr[, 2] 
##   1.7174501 
## > bias.full
##  (Intercept)       x[, 1]       x[, 2] 
## -0.029459219  0.002845078  0.010220794 
## > bias.cc
## (Intercept)     xr[, 1]     xr[, 2] 
##  -0.2544448  -0.3089315   2.7259508 
                                       
#  Results case (iii) MCAR -- CC estimator is consistent but inefficient

## > means
## (Intercept)      x[, 1]      x[, 2] (Intercept)     xr[, 1] 
##   19.970541    5.002845   -4.989779   20.045149    4.996496 
##     xr[, 2]             
##   -5.007878    0.622360 
## > sds
## (Intercept)      x[, 1]      x[, 2] (Intercept)     xr[, 1] 
##  2.07850830  0.19828353  1.18444484  2.58257942  0.24247188 
##     xr[, 2]             
##  1.46253837  0.03443423 
## > ses
## (Intercept)      x[, 1]      x[, 2] (Intercept)     xr[, 1] 
##   2.0352825   0.1898984   1.1592384   2.5951163   0.2425868 
##     xr[, 2] 
##   1.4775769 
## > bias.full
##  (Intercept)       x[, 1]       x[, 2] 
## -0.029459219  0.002845078  0.010220794 
## > bias.cc
##  (Intercept)      xr[, 1]      xr[, 2] 
##  0.045148986 -0.003504181 -0.007878482 
