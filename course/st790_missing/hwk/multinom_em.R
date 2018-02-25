#  EM algorithm for multinomial data -- ST 790 Homework 2

#  Functions

#  Convergence criterion 

converged <- function(tol,db,iter,imax){
	bmax <- max(abs(db))
	gmax <- iter==imax
	bmax <- bmax<tol
	converged <- gmax | bmax
	converged
}

#  EM algorithm -- this is specific to G=3 and H=2 

em.algorithm <- function(tol,Emax,thedata,robs,pr){

   iter <- 1
   dtheta <- 1
   N <- nrow(thedata)

#  Use complete case sample mean and sample covariance
#  matrix as the starting value

#  complete cases

   ccdata <- thedata[complete.cases(thedata),]

#  initial value of theta

   thetat <- c(sum(ccdata[,1]==1 & ccdata[,2]==1),
               sum(ccdata[,1]==2 & ccdata[,2]==1),
               sum(ccdata[,1]==3 & ccdata[,2]==1),
               sum(ccdata[,1]==1 & ccdata[,2]==2),
               sum(ccdata[,1]==2 & ccdata[,2]==2),
               sum(ccdata[,1]==3 & ccdata[,2]==2))/nrow(ccdata)
               
#  Record the results of each iteration   

   itermat <- c(0,thetat)
   
  while (!converged(tol,dtheta,iter,Emax)){

    thetatmat <- matrix(thetat,3,2)
    thetacol <- thetatmat/apply(thetatmat,1,sum)
    thetarow <- t(t(thetatmat)/apply(thetatmat,2,sum))

    thetatnew <- matrix(0,3,2)

    for (i in 1:N){
        
#  E-Step -- calculate conditional expectations
         
       indmat <- matrix(0,3,2)
       indmatcol <- matrix(0,3,2)
       indmatrow <- matrix(0,3,2)

       zi <- unlist(thedata[i,])
       if (sum(is.na(zi))==0) {indmat[zi[1],zi[2]] <- 1 }
       if (is.na(zi[1])==0) { indmatcol[zi[1],] <- c(1,1) }
       if (is.na(zi[2])==0) { indmatrow[,zi[2]] <- c(1,1,1) }
        
       thetatmati <- (robs[i,1]*robs[i,2]*indmat +
                       robs[i,1]*(1-robs[i,2])*thetacol*indmatcol +
                           (1-robs[i,1])*robs[i,2]*thetarow*indmatrow)

#  M-step calculate update sequentially       

       thetatnew <- thetatnew + thetatmati
   }
       
   thetanew <- matrix(thetatnew,1,6)/N

#  If requested, print the results of each iteration
    
   if (pr==1) print(c(iter,thetanew))

#  Record this iteration result

   itermat <- rbind(itermat,c(iter,thetanew))
    
# Compute change

    thetat <- (abs(thetat)<0.001)*0.001 + (abs(thetat) > 0.001)*thetat
    dtheta <- (abs(thetat)>0.001)*(thetanew-thetat)/thetat+(abs(thetat)<=0.001)*(thetanew-thetat)
    thetat <- thetanew
    iter <- iter+1

  }

#  final estimates

  results <- list(theta.em=thetat,itermat=itermat)
  return(results)

 }

#  Nonparametric bootstrap

np.bootstrap <- function(B,thedata,robs){

  theta.boot <- NULL

  for (b in 1:B){

#  Sample with replacement from the original data set

    brow <- sample(N,N,replace=TRUE)

    thedata.b <- thedata[brow,]
    robs.b <- robs[brow,]

#  Call the EM algorithm

    em.b <- em.algorithm(tol,Emax,thedata.b,robs.b,0)

    theta.b <- em.b$theta.em

    theta.boot <- rbind(theta.boot,theta.b)

  }

    return(list(theta.boot=theta.boot))
}

#############################################

#  Start of program                     

#  pr = 1 to print each EM iteration

  pr <- 1

#  Set EM parameters

Emax <- 200      #  max number of EM iterations
tol <- 1e-4    #  convergence tolerance

#  Read in the data

thedata <- read.table("multinom.dat")
colnames(thedata) <- c("y1","y2")

#  Create matrix of missingness indicators

robs <- matrix(as.numeric(!is.na(thedata)),ncol=2)

#  Call the EM algorithm

em.estimate <- em.algorithm(tol,Emax,thedata,robs,pr)

#  Final EM estimates

theta.em <- em.estimate$theta.em

#  Now use a nonparametric bootstrap to get standard errors

B <- 250    #  Number of bootstrap samples

em.boot <- np.bootstrap(B,thedata,robs)
theta.boot <- em.boot$theta.boot
boot.mean <- apply(theta.boot,2,mean)
boot.se <- apply(theta.boot,2,sd)

#  final results in R session and print results
#  to a file

if (pr==1){
    print("final estimates")
    print(round(theta.em,6))

##          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
## [1,] 0.249733 0.334381 0.085944 0.044937 0.182569 0.102436
    
    print("Bootstrap standard errors")
    print(boot.se)

## [1] 0.03326904 0.04019921 0.02517883 0.01986783 0.03852786
## [6] 0.02733862

    outfile="multinom_em.Rout"
    
    cat("EM Algorithm for Multinomial Data", file=outfile,"\n","\n",append=FALSE)
    cat("History of iterations",file=outfile,"\n","\n",append=TRUE)
    write.table(round(em.estimate$itermat,6),file=outfile,col.names=FALSE,row.names=FALSE,append=TRUE)
    cat(file=outfile,"\n","\n",append=TRUE)    
    cat("Bootstrap standard errors for theta",file=outfile,"\n","\n",append=TRUE)
    cat(boot.se,file=outfile,"\n",append=TRUE)
    
}


