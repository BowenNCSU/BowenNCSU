#  Use the "norm" package downloaded from CRAN to
#  carry out proper multiple imputation using a 
#  multivariate normal model to impute on the 
#  bivariate normal data 
#
#  install.packages("norm") and choose a mirror site
#
#  Visit the CRAN website for the documentation

#  Load the libary for the norm package

library(norm)
library(nlme)

#  number of imputed data sets

M <- 10

#  Read in the data set -- the missing values are denoted
#  by the usual R convention NA
#  The norm wants the data in a matrix, so we read it in
#  directly to a matrix

datamat <- matrix(scan("bvnormal.R.dat"),ncol=2,byrow=TRUE)

N <- nrow(datamat)
nvar <- ncol(datamat)

#  Run prelim.norm() to set up the data for the algorithm

prelim.mi <- prelim.norm(datamat)

#  Look at the numbers of missing values and missing data patterns

prelim.mi$nmis
prelim.mi$r

#  This produces

#  > prelim.mi$nmis
#    [1] 214 238
#
#> prelim.mi$r
#    [,1] [,2]
#  548    1    1
#  214    0    1
#  238    1    0

#  Run the EM algorithm initially with the default settings to 
#  get the form of the starting values

theta.init <- em.norm(prelim.mi,showits=FALSE)

#  Extract the parameters so you can see the format

theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix

#  Convert these to the form used by the function em.norm

theta.init <- makeparam.norm(prelim.mi,theta.init)

#  Create M imputed data sets using proper imputation
#  and carry out the desired analysis on each one; save
#  them in a big matrix in case you want to use them 
#  for something else

imputed.data <- NULL

#  the mi.inference function wants the parameters and their
#  standard errors saved in lists.  We'll also save the full
#  asymptotic covariance matrices for use with our own multivariate
#  mi inference function below

mu.list <- vector("list",M)
ses.list <- vector("list",M)
covs.list <- vector("list",M)
covparms.list <- vector("list",M)

#   Set the seed for the MCMC

rngseed(82)

for (m in 1:M){

#  Impute the missing values to create a single data set using MCMC.
#  Use the da.norm (for "data augmentation") function in the norm package
#  to impute the missing values. Here, steps is apparently the number
#  of iterations in the chain between draws - the documentation does not give
#  details, and it is not clear if there is a burn-in period
    
   imp.dat <- da.norm(prelim.mi,theta.init,steps=200,return.ymis=TRUE)   

#  misobs is a vector containing the imputed missing values in the order
#  they are missing in the original data set, so insert them into a copy
#  of the original data set.
   
   misobs <- imp.dat$ymis
   this.imp <- datamat
   this.imp[is.na(this.imp)] <- misobs
   
#  Save this imputed data set in case you want to do other stuff with it

   imputed.data <- rbind(imputed.data,this.imp)

#  Reconfigure the imputed data set to 1 record per observation 
#  for use with gls()
   
   this.imp.alt <- NULL
   for (i in 1:N){
       this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
       this.imp.alt <- rbind(this.imp.alt,this)
   }
     
#  Call gls() to fit the multivariate normal model to the reconfigured
#  imputed data set

   id <- factor(this.imp.alt[,1])
   ind <- factor(this.imp.alt[,2])
   y <- this.imp.alt[,3]
   time <- rep(seq(1,2,1),N)
   
   gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                  weights=varIdent(form= ~1 | ind),method="ML")

#  Get the estimates of the mean parameters and their standard errors
#  and also the full asymptotic covariance matrices
   
   this.mu <- gls.fit$coef
   this.ses <- sqrt(diag(gls.fit$varBeta))
   this.cov <- gls.fit$varBeta

#  Save these in the lists, as the mi.inference function
#  requires that the parameter estimates from each imputed data set and
#  their standard errors be in lists
   
   mu.list[[m]] <- this.mu
   ses.list[[m]] <- this.ses
   covs.list[[m]] <- this.cov
 
#  Save the estimates of the variance and covariance parameters   

   this.covmat <- getVarCov(gls.fit)
   this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])

   covparms.list[[m]] <- this.covparms
   
#  Unfortunately, getting their standard errors is hard -- they are computed
#  but are available in a different (and not clear) parameterization via 

#  gls.fit$apVar

#  We won't bother trying to extract them; usually these aren't of
#  interest anyway 

}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){

#  Get the mean of estimates and mean of (within) covariance matrices 
#  (presumably these lists are of the same length; don't bother
#  checking this)

    m <- length(est)
    qmat <- simplify2array(est)
    qbar <- apply(qmat,1,mean)
    wcov <- Reduce('+',covmat)/m

#  Get among-imputation covariance matrix    
    
    bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)

#   Rubin covariance matrix and diagonal elements of each component
    
    qcovmat <- wcov + (1+1/m)*bcov
    bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
    ubar <- diag(wcov)  

#   This code is from mi.inference - CIs, DFs, etc for each component

    tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
    rem <- (1 + (1/m)) * bm/ubar
    nu <- (m - 1) * (1 + (1/rem))^2
    alpha <- 1 - (1 - confidence)/2
    low <- qbar - qt(alpha, nu) * sqrt(tm)
    up <- qbar + qt(alpha, nu) * sqrt(tm)
    pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
    fminf <- (rem + 2/(nu + 3))/(rem + 1)

#   First 3 elements are the mean of estimates, their SEs,   
#   entire covariance matrix using Rubin's formula; last 2
#   are the within and among components 
    
    result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                   df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
    result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
      mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

#> within.cov
#              ind1         ind2
# ind1 0.0010229366 0.0005132697
# ind2 0.0005132697 0.0011069902
#> between.cov
#              ind1         ind2
# ind1 3.222650e-04 6.834349e-06
# ind2 6.834349e-06 2.179478e-04
#> Rubin.cov
#              ind1         ind2
# ind1 0.0013774280 0.0005207875
# ind2 0.0005207875 0.0013467327



