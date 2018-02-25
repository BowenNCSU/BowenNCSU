##################################################################
#
#  Use the mice package to carry out multivariate imputation
#  via chained equations on data from a hypertension study.  The
#  columns are (note that there is no id variable):
#
#   1-3     hypertensive status (0 = no, 1 = yes) at 1, 2, and 4
#           months after starting treatment
#   4       age range (1 = < 40, 2 = 40 - 65, 3 = > 65 years)
#   5       cholesterol (mg/DL) at baseline
#
##################################################################

#  You will need to install.packages("mice").  Documentation is
#  available on CRAN -- also see the excellent Journal of Statistical
#  Software (JSS) article (link on the course website).  Also need gee
#  package to carry out longitudinal analysis of each imputed data
#  set and the MASS package to get the function polr() used to 
#  implement equations for ordered categorical variables 

library(mice)
library(gee)
library(MASS)
library(reshape)
library(norm)

#  Number of imputed data sets

M <- 10

#  Read in the data 

hyperdat <- read.table("hyper.R.dat")
colnames(hyperdat) <- c("h1","h2","h4","age","chol")

#  create an id variable for model fits later

N <- nrow(hyperdat)
id <- seq(1,N,1)

#  convert columns 1-4 to factors for later use by mice() and check
#  on the result

hyperdat[c("h1","h2","h4","age")] <- lapply(hyperdat[c("h1","h2","h4","age")],as.factor)
str(hyperdat)

#'data.frame':	300 obs. of  5 variables:'
# $ h1  : Factor w/ 2 levels "0","1": 1 2 1 2 2 2 2 1 2 1 ...
# $ h2  : Factor w/ 2 levels "0","1": 1 2 2 2 2 2 1 1 2 1 ...
# $ h4  : Factor w/ 2 levels "0","1": 1 2 NA 1 1 1 2 1 2 1 ...
# $ age : Factor w/ 3 levels "1","2","3": 1 1 3 1 3 2 NA 2 1 1 ...
# $ chol: num  NA 185 264 227 161 ...

#  View missing data patterns and numbers of missing values

md.pattern(hyperdat)

#  Call mice to create m=M imputed data sets.  Here, we set the argument
#  set.seed() using the seed= argument so that we can reproduce the results.
#  The method= argument allows the user to set the type of equation used
#  for each column (there are defaults; see the documentation).  Here
#  we use logistic regression for the 3 binary hypertension indicators,
#  an ordered (cumulative) logit model for the categorical age variable,
#  and Bayesian normal imputation for the continuous variable
#  chol. 

hyper.imp <- mice(hyperdat,m=M,seed=3485,method=c("logreg","logreg","logreg","polr","norm"),maxit=20)
        
#  In the following, the PredictorMatrix shows that each variable was imputed
#  by regressing it (using the method in the method= argument) on all others.
#  This is the default and can be changed by specifying the
#  predictorMatrix= argument -- see the documentation and JSS paper.

# > hyper.imp
# Multiply imputed data set
# Call:
# mice(data = hyperdat, m = M, method = c("logreg", "logreg", 
#    "logreg", "polr", "norm"), seed = 3485)
# Number of multiple imputations:  10
# Missing cells per column:
#    h1   h2   h4  age chol 
#    28   24   74   57  120 
# Imputation methods:
#      id       h1       h2       h4      age     chol 
#      "" "logreg" "logreg" "logreg"   "polr"   "norm" 
# VisitSequence:
#  h1   h2   h4  age chol 
#   1    2    3    4    5 
# PredictorMatrix:
#      h1 h2 h4 age chol
# h1    0  1  1   1    1
# h2    1  0  1   1    1
# h4    1  1  0   1    1
# age   1  1  1   0    1
# chol  1  1  1   1    0

# Random generator seed value:  3485 

#  Can assess convergence by plotting the m parallel
#  imputation streams

#  plot(hyper.imp,c("h1","h2","h4"))
#  plot(hyper.imp,c("age","chol"))

#  Can look at the imputed data sets via the complete()
#  function -- here we look at the first few rows of the
#  imputed data set 1 next to the original data set

cbind(head(complete(hyper.imp,1)),head(hyperdat))

#   h1 h2 h4 age     chol h1 h2   h4 age  chol
# 1  0  0  0   1 163.4881  0  0    0   1    NA
# 2  1  1  1   1 185.1000  1  1    1   1 185.1
# 3  0  1  0   3 263.9000  0  1 <NA>   3 263.9
# 4  1  1  0   1 227.4000  1  1    0   1 227.4
# 5  1  1  0   3 160.6000  1  1    0   3 160.6
# 6  1  1  0   2 206.0000  1  1    0   2 206.0

#  all.imp contains all M imputed data sets, stacked; the 
#  first column is the imputation number; extracting the
#  data sets this way allows them to be exported for other
#  uses

all.imp <- complete(hyper.imp,action='long')

head(all.imp)

#   .imp .id h1 h2 h4 age     chol
# 1    1   1  0  0  0   1 163.4881
# 2    1   2  1  1  1   1 185.1000
# 3    1   3  0  1  0   3 263.9000
# 4    1   4  1  1  0   1 227.4000
# 5    1   5  1  1  0   3 160.6000
# 6    1   6  1  1  0   2 206.0000

#  If we were just interested in a simple analysis, say the 
#  regression of chol on age, carrying this out on each imputed
#  data set and combining the results is accomplished 
#  trivially using the with() and pool() functions applied
#  to the mice object hyper.imp

simple.anal <- with(hyper.imp,lm(chol ~ age))

#  summary(simple.anal) will show the summary for each imputed
#  data set

simple.anal.combined <- pool(simple.anal)
summary(simple.anal.combined)

#                est       se         t       df     Pr(>|t|)     lo 95
# (Intercept) 180.99017 3.493947 51.801057 39.70541 0.000000e+00 173.92701
# age2         33.30364 7.362016  4.523712 79.54964 2.097131e-05  18.65149
# age3         48.75391 6.924915  7.040363 41.21562 1.405866e-08  34.77099
#                 hi 95 nmis       fmi    lambda
# (Intercept) 188.05334   NA 0.4447712 0.4174909
# age2         47.95579   NA 0.2851375 0.2673879
# age3         62.73684   NA 0.4351675 0.4084081

#  However, we are interested in a longitudinal analysis of the 
#  binary hypertension status using gee().  Accordingly, we need
#  to reconfigure the data to 1 record per observation.

#  One way to do this is convert the imputed data to a list of M 
#  imputed data sets, reconfigure, analyze each, and combine yourself.
#  Add the id variabe back in.

hyper.imp.data <- as.list(1:M)
for (i in 1:M){
    hyper.imp.data[[i]] <- cbind(id,complete(hyper.imp,action=i))
}

#  Reconfigure each data set and add a month variable

hyper.imp.data2 <- lapply(hyper.imp.data,reshape,varying=2:4,v.names="hyper",idvar="id",timevar="month",
       times=c(1,2,4),direction="long")
hyper.imp.data2 <- lapply(hyper.imp.data2,FUN=function(u){ u[order(u$id),]})

#  Here is the top of the first imputed data set
#  id age     chol month hyper
# 1.1  1   1 163.4881     1     0
# 1.2  1   1 163.4881     2     0
# 1.4  1   1 163.4881     4     0
# 2.1  2   1 185.1000     1     1
# 2.2  2   1 185.1000     2     1
# 2.4  2   1 185.1000     4     1

#  Do the gee analysis on each data set

hyper.imp.gee <- lapply(hyper.imp.data2, FUN=function(u){
gee(hyper ~ age + month+ chol,id=id,family=binomial,corstr="unstructured",data=u)
})

#  We can now extract whatever we want from the gee fits; e.g., 
#  here are the parameter estimates from the M analyses

hyper.imp.parmmat <- sapply(hyper.imp.gee,coefficients)
hyper.imp.parmmat

# > hyper.imp.parmmat
#                    [,1]        [,2]        [,3]          [,4]        [,5]
# (Intercept)  4.29401379  4.07411914  3.89640372  5.41557424  3.71879073
# age2        -0.14392260 -0.11918271  0.09596662  0.23988734 -0.16512978
# age3         0.39409536  0.13003971  0.28099951  0.60825267  0.02166555
# month       -0.68150710 -0.66786870 -0.65818956 -0.73021995 -0.64075761
# chol        -0.01384963 -0.01233472 -0.01179351 -0.01886716 -0.01059747
#                    [,6]       [,7]        [,8]        [,9]       [,10]
# (Intercept)  5.12406843  4.6248814  4.61118855  4.56800446  4.55831333
# age2         0.11836326  0.2261034 -0.26339066 -0.15884231 -0.14745528
# age3         0.45115187  0.3190503  0.28755927  0.02574516  0.37925325
# month       -0.71647125 -0.7621449 -0.70871770 -0.69341044 -0.71965907
# chol        -0.01766432 -0.0147185 -0.01464137 -0.01447368 -0.01445712

minparm <- apply(hyper.imp.parmmat,1,min)
maxparm <- apply(hyper.imp.parmmat,1,max)

#  Put into lists and combine them 

hyper.imp.parms <- lapply(hyper.imp.gee,coefficients)
hyper.imp.ses <- lapply(hyper.imp.gee,FUN=function(u){
    sqrt(diag(u$naive.variance))
})

#  Note that one can actually extract the entire asymptotic
#  variance covariance matrix from eac fit (naive.variance),
#  one could do multivariate inference if one chose; we just
#  get the standard errors

#  We can't use the mice pool() function on these, but we
#  CAN use the mi.inference function in the norm package

mi.parms <- mi.inference(hyper.imp.parms,hyper.imp.ses,confidence=0.95)
mi.results <- cbind(mi.parms$est,mi.parms$std.err,mi.parms$lower,mi.parms$upper,
      mi.parms$df,minparm,maxparm)
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

# > mi.results
#                     Est      StdErr       Lower        Upper        DF
# (Intercept)  4.48853577 0.750291363  2.95934488  6.017726665 31.52822
# age2        -0.03176027 0.329851751 -0.68830491  0.624784371 79.06851
# age3         0.28978127 0.315262766 -0.34079363  0.920356171 60.20668
# month       -0.69789462 0.067928955 -0.83290997 -0.562879281 87.03704
# chol        -0.01433975 0.003632308 -0.02173164 -0.006947858 32.77800
#                     Min         Max
# (Intercept)  3.71879073  5.41557424
# age2        -0.26339066  0.23988734
# age3         0.02166555  0.60825267
# month       -0.76214485 -0.64075761
# chol        -0.01886716 -0.01059747

#  If we are willing to do a little of our own programming, though, we
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

#  Get the covariance matrices in a list

hyper.imp.covs <- lapply(hyper.imp.gee,FUN=function(u){
    (u$naive.variance)})

#  Call the multivariate function

mi.mv.parms <-
    mi.mv.inference(hyper.imp.parms,hyper.imp.covs,confidence=0.95)

#  This should be identical to the results from mi.inference

mi.mv.results <- cbind(mi.mv.parms$est,mi.mv.parms$std.err,mi.mv.parms$lower,mi.mv.parms$upper,
      mi.mv.parms$df,minparm,maxparm)
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#                     Est      StdErr       Lower        Upper        DF
# (Intercept)  4.48853577 0.750291363  2.95934488  6.017726665 31.52822
# age2        -0.03176027 0.329851751 -0.68830491  0.624784371 79.06851
# age3         0.28978127 0.315262766 -0.34079363  0.920356171 60.20668
# month       -0.69789462 0.067928955 -0.83290997 -0.562879281 87.03704
# chol        -0.01433975 0.003632308 -0.02173164 -0.006947858 32.77800
#                     Min         Max
# (Intercept)  3.71879073  5.41557424
# age2        -0.26339066  0.23988734
# age3         0.02166555  0.60825267
# month       -0.76214485 -0.64075761
# chol        -0.01886716 -0.01059747

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.parms$within
between.cov <- mi.mv.parms$between
Rubin.cov <- mi.mv.parms$cov.mat

# > within.cov
#             (Intercept)          age2          age3         month      
# age2         0.022538258  0.0720944879  0.0209816660  4.357840e-05
# age3         0.042259772  0.0209816660  0.0609629211 -5.457212e-04
# month       -0.012318580  0.0000435784 -0.0005457212  3.130530e-03
# chol        -0.001208678 -0.0001898994 -0.0002903326  2.617331e-05
#                     chol
# (Intercept) -1.208678e-03
# age2        -1.898994e-04
# age3        -2.903326e-04
# month        2.617331e-05
# chol         6.280203e-06
# > between.cov
#              (Intercept)         age2          age3         month          chol
# (Intercept)  0.273425315  0.046765487  0.0715843459 -1.512727e-02 -1.303554e-03
# age2         0.046765487  0.033370627  0.0198782331 -3.410499e-03 -2.367000e-04
# age3         0.071584346  0.019878233  0.0349342640 -4.115017e-03 -3.620864e-04
# month       -0.015127273 -0.003410499 -0.0041150172  1.348921e-03  6.897981e-05
# chol        -0.001303554 -0.000236700 -0.0003620864  6.897981e-05  6.284963e-06
# > Rubin.cov
#              (Intercept)          age2          age3         month
# (Intercept)  0.562937129  0.0739802935  0.1210025522 -0.0289585804
# age2         0.073980293  0.1088021777  0.0428477224 -0.0037079702
# age3         0.121002552  0.0428477224  0.0993906116 -0.0050722401
# month       -0.028958580 -0.0037079702 -0.0050722401  0.0046143429
# chol        -0.002642588 -0.0004502694 -0.0006886276  0.0001020511
#                      chol
# (Intercept) -2.642588e-03
# age2        -4.502694e-04
# age3        -6.886276e-04
# month        1.020511e-04
# chol         1.319366e-05



