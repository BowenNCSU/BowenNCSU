#  Simulation of estimators for the mean outcome 
#  under various forms of missingness --EXAMPLE 1 
#  of Chapter 1 of the notes

#  Specify simulation set-up and scenario 

S <- 1000    #  Number of simulated data sets
N <- 100     #  Sample size
set.seed(4)  #  Set the seed

#  The parameter vector psi specifies the missingness 
#  mechanism.  See the data generation function below 
#  for the interpretation of elements of psi

#psi <- c(1.5,0,0,-0.05)    #  MNAR scenario
psi <- c(0.8,-1.4,0.5,0)  #  MAR scenario
#psi <- c(1.0,0,0,0)       #  MCAR scenario

#  Data generation function -- outcome y and auxiliary 
#  variables v, and missingness indicator c = 1 if 
#  r = (1,1), c = 0 if r = (0,1)

#  In these scenarios, v is always observed, y is 
#  missing

generate <- function(N,psi){

#   Parameter values
    
    pv <- 0.5
    theta <- c(2.3,0.5,1)
   
#   Generate auxilary variables and outcome
    
    v1 <- rnorm(N,0,1)       #  continuous auxiliary 
    v2 <- rbinom(N,1,pv)     #  binary auxiliary
    v <- cbind(v1,v2)
    y <- exp(theta[1]+theta[2]*v1+theta[3]*v2+rnorm(N,0,0.5))    #  outcome

#   Generate missingness indicators
    
    elinpred <- exp(psi[1]+psi[2]*v1+psi[3]*v2+psi[4]*y)       
    pifunc <-  elinpred/(elinpred+1)      

    c <- rbinom(N,1,pifunc)

    thedata <- list(y=y,v=v,c=c,pifunc=pifunc)
    return(thedata)
}

#   Calculate true value of mean outcome by simulation
#   (1 million)

simtrue <- generate(1e6,psi)
truemean <- mean(simtrue$y)

#  Simulation loop

results <- NULL

for (s in 1:S){

#  Generate a data set

    thedat <- generate(N,psi)
    y <- thedat$y
    v <- thedat$v
    c <- thedat$c
    pifunc <- thedat$pifunc

#   Ideal full data sample mean -- this is of 
#   course not possible in practice but serves 
#   as a "gold standard"     

    full <- mean(y)

#   Complete case sample mean

    Nobs <- sum(c)
    ccmean <- sum(c*y)/Nobs

#   Inverse weighted mean with correct missingness 
#   probabilities as weights.  In MNAR case, inverse
#   weighted estimators are not possible in practice    
#   because they depend on the unobserved outcome, 
#   so are of academic interest only

    ipwmean <- mean(c*y/pifunc)

#   Alternative inverse probability weighted estimator 
#   that is more efficient

    ipwmean.alt <- sum(c*y/pifunc)/sum(c/pifunc)

#   Proportion of sample size N observed

    pobs <- Nobs/N
    
#   Save results

    results <- rbind(results,c(full,ccmean,ipwmean,ipwmean.alt,pobs))

}

#   End of simulation loop
    
#   Summarize simulation results -- can look at
#   truemean, means, sds to see the extent of
#   bias and average proportion of observed
#   and compare mses to get an idea of relative
#   efficiency

means <- apply(results,2,mean)
sds <- apply(results,2,sd)

bias <- means[1:4]-truemean
mses <- apply((results[,1:4]-truemean)^2,2,mean)

