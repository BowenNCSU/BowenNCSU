#  EM algorithm for bivariate normal data -- 
#  EXAMPLE 2 of Chapter 3 of the notes

#  Simulate a bivariate normal data set
#  for demonstration (or read in your own data
#  set here)

set.seed(4)   #  Set the seed
N <- 1000     #  Sample size
Emax <- 100      #  max number of EM iterations
tol <- 1e-5    #  convergence tolerance

#  Functions

#  Convergence criterion 

converged <- function(tol,db,iter,imax){
	bmax <- max(abs(db))
	gmax <- iter==imax
	bmax <- bmax<tol
	converged <- gmax | bmax
	converged
}

#  EM algorithm

em.algorithm <- function(tol,Emax,fulldata,robs,pr){

   iter <- 1
   dtheta <- 1
   N <- nrow(fulldata)

#  Change NAs to nonsensical value so can do regular
#  operations

   fulldata[is.na(fulldata)] <- -999

#  Use complete case sample mean and sample covariance
#  matrix as the starting value

#  Get complete cases

   ccdata <- fulldata[complete.cases(fulldata),]
   ccN <- nrow(ccdata)
   
   mu0 <- apply(ccdata,2,mean)
   Sigma0 <- t(sweep(data.matrix(ccdata),2,mu0))%*%sweep(data.matrix(ccdata),2,mu0)/(ccN-1)

#  initial value of theta

   thetat <- c(mu0,Sigma0[1,1],Sigma0[1,2],Sigma0[2,2])

#  Record the results of each iterations   

   itermat <- c(0,thetat)
   
  while (! converged(tol,dtheta,iter,Emax)){

    mut <- thetat[1:2]
    Sigmat <- matrix(c(thetat[3],thetat[4],thetat[4],thetat[5]),2,2)
  
#  E-Step -- calculate conditional expectiations

    ey1 <- (robs[,1]==1)*fulldata[,1] + (robs[,1]==0)*(mut[1]+ Sigmat[1,2]*(fulldata[,2]-mut[2])/Sigmat[2,2])
    ey2 <- (robs[,2]==1)*fulldata[,2] + (robs[,2]==0)*(mut[2]+ Sigmat[1,2]*(fulldata[,1]-mut[1])/Sigmat[1,1])
    ey1.2 <- (robs[,1]==1)*fulldata[,1]^2 + (robs[,1]==0)*(Sigmat[1,1]-Sigmat[1,2]^2/Sigmat[2,2]+ ey1^2)
    ey2.2 <- (robs[,2]==1)*fulldata[,2]^2 + (robs[,2]==0)*(Sigmat[2,2]-Sigmat[1,2]^2/Sigmat[1,1]+ ey2^2)
    ey1y2 <- (robs[,1]*robs[,2]*fulldata[,1]*fulldata[,2] + robs[,1]*(1-robs[,2])*fulldata[,1]*ey2
         + (1-robs[,1])*robs[,2]*ey1*fulldata[,2])

#  M-step -- calculate update

   munew <- c(mean(ey1),mean(ey2))
   s11 <- mean(ey1.2)-munew[1]^2
   s22 <- mean(ey2.2)-munew[2]^2
   s12 <- mean(ey1y2)-munew[1]*munew[2]

   thetanew <- c(munew,s11,s12,s22)

#  If requested, print the results of each iteration
    
   if (pr==1) print(c(iter,thetanew))

#  Record this iteration result

   itermat <- rbind(itermat,c(iter,thetanew))
    
# Compute relative change

    dtheta <- (thetanew-thetat)/thetat
    thetat <- thetanew
    iter <- iter+1

  }

#  final estimates

  results <- list(theta.em=thetat,itermat=itermat)
  return(results)

 }

#  Nonparametric bootstrap

np.bootstrap <- function(B,fulldata,robs){

  theta.boot <- NULL

  for (b in 1:B){

#  Sample with replacement from the original data set

    brow <- sample(N,N,replace=TRUE)

    fulldata.b <- fulldata[brow,]
    robs.b <- robs[brow,]

#  Call the EM algorithm

    em.b <- em.algorithm(tol,Emax,fulldata.b,robs.b,0)

    theta.b <- em.b$theta.em

    theta.boot <- rbind(theta.boot,theta.b)

  }

    return(list(theta.boot=theta.boot))
}

#############################################

#  Start of program                     

#  True values of theta = (mu1,mu2,sigma1^2,sigma12,sigma2^2)

theta <- c(5,8,1,0.5,1)

mu <- theta[1:2]
Sigma <- matrix(c(theta[3],theta[4],theta[4],theta[5]),2,2)

#  Missingness mechanism -- 3 possibilities
#  r = (1,1)  y1 and y2 both observed
#  r = (1,0)  y1 observed, y2 missing
#  r = (0,1)  y1 missing, y2 observed

#  To generate nonomontone MAR data is not
#  trivial -- we use the "randomized monotone
#  missingness process" approach in Robins
#  and Gill (1997), Statistics in Medicine

probobs <- 0.6     #  prob observe y1

psi1 <- c(-1,0.3)
psi2 <- c(-1,0.1)

#  Simulation loop -- set S=1 for single data set
#  If S > 1 set pr = 0 so it won't print internal
#  EM algorithm iterations

S <- 1
pr <- 1

theta.sim <- NULL
theta.boot.sim <- NULL
se.boot.sim <- NULL

for (s in 1:S){

fulldata <- NULL
robs <- NULL

#  Generate a data set

for (i in 1:N){

#  Generate observation indicator  

  obsind <- rbinom(1,1,probobs)
  
#  Generate full data

  y <- mu + t(chol(Sigma))%*%rnorm(2,0,1)

#  Create r = (r1,r2)

  ep1 <- exp(psi1[1]+psi1[2]*y[1])
  p1 <- ep1/(1+ep1)
  ep2 <- exp(psi2[1]+psi2[2]*y[2])
  p2 <- ep2/(1+ep2)
  
  r1 <- obsind + (1-obsind)*rbinom(1,1,p2)
  r2 <- (1-obsind) + obsind*rbinom(1,1,p1)  

  fulldata <- rbind(fulldata,t(y))
  robs <- rbind(robs,c(r1,r2))

}

#   To make like real data, round to 3 decimal
#   places and replace positions where y1 or y2 
#   is missing with NA 

fulldata <- round(fulldata,3)

fulldata[robs==0] <- NA

#  Write the data set to a text file for use in SAS, etc

#write.table(fulldata,"bvnormal.dat",row.names=FALSE,col.names=FALSE)

#  Call the EM algorithm

em.estimate <- em.algorithm(tol,Emax,fulldata,robs,pr)

theta.em <- em.estimate$theta.em

#  Final EM estimates

 muhat <- theta.em[1:2]
 Sigmahat <- matrix(c(theta.em[3],theta.em[4],theta.em[4],theta.em[5]),2,2)

 #  Now use a nonparametric bootstrap to get standard errors

B <- 250    #  Number of bootstrap samples

em.boot <- np.bootstrap(B,fulldata,robs)
theta.boot <- em.boot$theta.boot
boot.mean <- apply(theta.boot,2,mean)
boot.se <- apply(theta.boot,2,sd)

#  Append the results for this simulated data set

theta.sim <- rbind(theta.sim,theta.em)
theta.boot.sim <- rbind(theta.boot.sim,boot.mean)
se.boot.sim <- rbind(se.boot.sim,boot.se)

#  If this is not a simulation (S=1), print 
#  final results in R session and print results
#  to a file

if (pr==1){
    print("final estimates")
    print(muhat)
    print(Sigmahat)
    print("Bootstrap standard errors")
    print(boot.se)

    outfile="bvnormal_em.Rout"
    
   cat("EM Algorithm for Bivariate Normal Data", file=outfile,"\n","\n",append=FALSE)
   cat("History of iterations",file=outfile,"\n","\n",append=TRUE)
   write.table(round(em.estimate$itermat,6),file=outfile,col.names=FALSE,row.names=FALSE,append=TRUE)
   cat(file=outfile,"\n","\n",append=TRUE)    
   cat("Final estimate of mu",file=outfile,"\n","\n",append=TRUE)
   cat(muhat,file=outfile,"\n","\n",append=TRUE)
   cat("Final estimate of Sigma",file=outfile,"\n","\n",append=TRUE)
   write.table(round(Sigmahat,6),file=outfile,col.names=FALSE,row.names=FALSE,append=TRUE)
   cat(file=outfile,"\n","\n",append=TRUE)
   cat("Bootstrap standard errors for mu and distinct elements of Sigma",file=outfile,"\n","\n",append=TRUE)
   cat(boot.se,file=outfile,"\n",append=TRUE) 
}

}
