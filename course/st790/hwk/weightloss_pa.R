#########################################################################
#
#   Homework 2, Weightloss data
#
#########################################################################

library(nlme)
library(Matrix)
library(magic)

#  u <- model object from fit, m <- total # individuals

robust.cov <- function(u,m){
  form <-  formula(u)
  mf <- model.frame(form,getData(u))
  Xmat <- model.matrix(form,mf)
  Vlist <-  as.list(1:m)
  for (i in 1:m){
    Vlist[[i]] <- getVarCov(u,individual=i,type="marginal")
  }
  V <- Reduce(adiag,Vlist)
  Vinv <- solve(V)
  Sig.model <- solve(t(Xmat)%*%Vinv%*%Xmat)
  resid <- diag(residuals(u,type="response"))
  ones.list <- lapply(Vlist,FUN=function(u){matrix(1,nrow(u),ncol(u))})
  Ones <- Reduce(adiag,ones.list)
  meat <- t(Xmat)%*%Vinv%*%resid%*%Ones%*%resid%*%Vinv%*%Xmat
  Sig.robust <- Sig.model%*%meat%*%Sig.model
  se.robust <- sqrt(diag(Sig.robust))
  se.model <- sqrt(diag(Sig.model))
  return(list(Sig.model=Sig.model,se.model=se.model,Sig.robust=Sig.robust,se.robust=se.robust))
}
#  Read in the data -- they are in the "wide" form of one record per observation

thedat.wide <- read.table("weightloss.dat",row.names=NULL)
colnames(thedat.wide) <-  c("id","month0","month1","month2","month3","month4","program")

#  Reconfigure as one observation per record ("long" format) and sort
#  in id order

thedat <- reshape(thedat.wide,varying=c("month0","month1","month2","month3"
,"month4"),v.names="weight",idvar="id",times=c(0,3,6,9,12),timevar="month",direction="long")
thedat <- thedat[order(thedat$id),]

m <- max(thedat$id)

month2 <- thedat$month^2
diet <- as.factor(thedat$program)

# time variable for specifying covariance matrix - data are balanced
# so not needed

thedat$time <- as.numeric(factor(thedat$month))

#  Fit different covariance structures and a quadratic model to allow for
#  possibly nonconstant rate of change.  From the diagnostics we did in 
#  Homework 1, there seemed to be little evidence that variance changes 
#  over time and only a hint that variance might be a tad bit smaller
#  for the control group.  A common compound symmetric structure seemed
#  plausible for groups 2 and 3, with a somewhat different structure for 
#  group 1.  As we know, we cannot fit separate covariance structures in 
#  gls().  The only thing we can allow to vary with group is the variance.  
#  However, with a separate quadratic model for each program
#  group, we could fit the model separately by group and then add up
#  the objective functions and compute AIC and BIC ourselves to
#  compare to the AIC and BIC for the common covariance fits.  But rather
#  than go to all that trouble, we can just use SAS...  

#  As in the SAS program, we maintain separate intercepts for each
#  program; you might have decided based on randomization and the
#  results to collapse to a single intercept for all groups. We fit a
#  series of models involving unstructured and compound symmetric
#  overall covariance structures and then compare via AIC and BIC.

#  increase number of iterations

weight.cntl <- glsControl(maxIter=200,msMaxIter=200,tolerance=1e-5,msTol=1e-6)

# As a first shot, fit the most general model possible with gls() with
# a different quadratic mean for each group, unstructured correlation
# with different variances for each group different at each time point
# Of course, in R we are stuck with having the correlation matrix being the
# same for each group, so we can't compare directly to SAS

weight.un <- gls(weight ~ -1 + diet + month:diet + month2:diet,data=thedat,
              correlation=corSymm(form = ~ 1 | as.factor(id)),
              weights = varIdent(form= ~ 1 | time*diet),method="ML",control=weight.cntl)
beta.un <- coef(weight.un)
sebeta.un <- robust.cov(weight.un,m)$se.model
V.un.program1 <- getVarCov(weight.un, individual=1)  #  or corMatrix(weight.un$modelStruct$corStruct)[[1]]
V.un.program2 <- getVarCov(weight.un, individual=35)
V.un.program3 <- getVarCov(weight.un, individual=63)
Gamma.un <- cov2cor(V.un)

## > cbind(beta.un,sebeta.un)
##                  beta.un  sebeta.un
## diet1        244.3937825 3.11264903
## diet2        244.5584215 3.45099378
## diet3        240.4860006 3.29073278
## diet1:month   -2.3485761 0.54010285
## diet2:month  -12.0225079 0.69513382
## diet3:month   -6.3248930 0.57232362
## diet1:month2   0.2223571 0.04606908
## diet2:month2   0.3641271 0.06204013
## diet3:month2   0.3009843 0.04928280

## > V.un.program1
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 367.26 285.20 303.72 386.67 307.20
## [2,] 285.20 452.19 298.66 314.73 290.65
## [3,] 303.72 298.66 377.98 351.62 304.28
## [4,] 386.67 314.73 351.62 587.42 341.70
## [5,] 307.20 290.65 304.28 341.70 389.46
##   Standard Deviations: 19.164 21.265 19.442 24.237 19.735 
## > V.un.program2
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 429.34 291.35 457.78 367.16 431.30
## [2,] 291.35 403.66 393.36 261.15 356.58
## [3,] 457.78 393.36 734.50 430.46 550.78
## [4,] 367.16 261.15 430.46 453.05 389.65
## [5,] 431.30 356.58 550.78 389.65 656.64
##   Standard Deviations: 20.721 20.091 27.102 21.285 25.625 
## > V.un.program3
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 458.93 322.00 401.27 428.35 392.92
## [2,] 322.00 461.27 356.51 315.01 335.88
## [3,] 401.27 356.51 527.97 411.83 411.47
## [4,] 428.35 315.01 411.83 576.90 387.44
## [5,] 392.92 335.88 411.47 387.44 509.85
##   Standard Deviations: 21.423 21.477 22.978 24.019 22.58 

## > Gamma.un
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]
## [1,] 1.00000 0.68876 0.81253 0.81813 0.80406
## [2,] 0.68876 1.00000 0.71849 0.61933 0.69167
## [3,] 0.81253 0.71849 1.00000 0.72429 0.78893
## [4,] 0.81813 0.61933 0.72429 1.00000 0.69538
## [5,] 0.80406 0.69167 0.78893 0.69538 1.00000

#  quadratic means, unstructured correlation, variances same for all groups
#  different at each time point; this is directly comparable to SAS with common
#  unstructured covariance matrix

weight.un.same <- gls(weight ~ -1 + diet + month:diet + month2:diet,data=thedat,
              correlation=corSymm(form = ~ 1 | as.factor(id)),
              weights = varIdent(form= ~ 1 | time),method="ML",control=weight.cntl)
beta.un.same <- coef(weight.un.same)
sebeta.un.same <- robust.cov(weight.un.same,m)$se.model
V.un.same <- getVarCov(weight.un.same, individual=1)  
Gamma.un.same <- cov2cor(V.un.same)

##  These results are identical to PROC MIXED with common unstructured
##  covariance structure

## > cbind(beta.un.same,sebeta.un.same)
##              beta.un.same sebeta.un.same
## diet1         244.0298249     3.31009425
## diet2         249.1228489     3.64754621
## diet3         240.9359801     3.13103570
## diet1:month    -2.4085748     0.58874947
## diet2:month   -11.3990354     0.64877032
## diet3:month    -6.4356272     0.55690125
## diet1:month2    0.2335804     0.05122961
## diet2:month2    0.3220838     0.05645228
## diet3:month2    0.3120334     0.04845836

## > V.un.same
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 413.42 292.76 373.70 390.76 366.47
## [2,] 292.76 434.40 329.07 296.95 318.88
## [3,] 373.70 329.07 506.48 382.51 394.23
## [4,] 390.76 296.95 382.51 545.60 361.04
## [5,] 366.47 318.88 394.23 361.04 498.68
##   Standard Deviations: 20.333 20.842 22.505 23.358 22.331 

## > Gamma.un.same
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]
## [1,] 1.00000 0.69083 0.81667 0.82276 0.80712
## [2,] 0.69083 1.00000 0.70156 0.60995 0.68512
## [3,] 0.81667 0.70156 1.00000 0.72766 0.78444
## [4,] 0.82276 0.60995 0.72766 1.00000 0.69216
## [5,] 0.80712 0.68512 0.78444 0.69216 1.00000

#  With gls(), we can also specify and unstructured correlation matrix
#  with the same variance at all time points, which can be either
#  different by group or the same for all groups.  We can't fit this in SAS.
#  Given that the diagnostics we did in Homework 1 seem to suggest 
#  that variance might be similar at all time points (perhaps with a different
#  variance for each group), fit these models

weight.un.var <- gls(weight ~ -1 + diet + month:diet + month2:diet,data=thedat,
              correlation=corSymm(form = ~ 1 | as.factor(id)),
              weights = varIdent(form= ~ 1 | diet),method="ML",control=weight.cntl)
beta.un.var <- coef(weight.un.var)
sebeta.un.var <- robust.cov(weight.un.var,m)$se.model
V.un.var.program1 <- getVarCov(weight.un.var,individual=1)
V.un.var.program2 <- getVarCov(weight.un.var,individual=35)
V.un.var.program3 <- getVarCov(weight.un.var,individual=63)
Gamma.un.var <- cov2cor(V.un.var)

## > cbind(beta.un.var,sebeta.un.var)
##              beta.un.var sebeta.un.var
## diet1        244.4266383    3.41567161
## diet2        248.8239701    4.20793061
## diet3        241.4134930    3.47397902
## diet1:month   -2.0802938    0.56877419
## diet2:month  -11.3383272    0.70070036
## diet3:month   -5.9960117    0.57848348
## diet1:month2   0.2099395    0.04857415
## diet2:month2   0.3172741    0.05984084
## diet3:month2   0.2803115    0.04940334

## > V.un.var.program1
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 429.20 295.61 348.74 351.14 345.10
## [2,] 295.61 429.20 308.37 265.81 296.86
## [3,] 348.74 308.37 429.20 310.86 338.61
## [4,] 351.14 265.81 310.86 429.20 298.45
## [5,] 345.10 296.86 338.61 298.45 429.20
##   Standard Deviations: 20.717 20.717 20.717 20.717 20.717 
## > V.un.var.program2
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 536.44 369.47 435.87 438.87 431.33
## [2,] 369.47 536.44 385.43 332.23 371.04
## [3,] 435.87 385.43 536.44 388.54 423.21
## [4,] 438.87 332.23 388.54 536.44 373.03
## [5,] 431.33 371.04 423.21 373.03 536.44
##   Standard Deviations: 23.161 23.161 23.161 23.161 23.161 
## > V.un.var.program3
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 496.21 341.76 403.18 405.96 398.98
## [2,] 341.76 496.21 356.52 307.32 343.21
## [3,] 403.18 356.52 496.21 359.40 391.47
## [4,] 405.96 307.32 359.40 496.21 345.05
## [5,] 398.98 343.21 391.47 345.05 496.21
## Standard Deviations: 22.276 22.276 22.276 22.276 22.276

## > Gamma.un.var
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]
## [1,] 1.00000 0.68876 0.81253 0.81813 0.80406
## [2,] 0.68876 1.00000 0.71849 0.61933 0.69167
## [3,] 0.81253 0.71849 1.00000 0.72429 0.78893
## [4,] 0.81813 0.61933 0.72429 1.00000 0.69538
## [5,] 0.80406 0.69167 0.78893 0.69538 1.00000

weight.un.var.same <- gls(weight ~ -1 + diet + month:diet + month2:diet,data=thedat,
              correlation=corSymm(form = ~ 1 | as.factor(id)),
              method="ML",control=weight.cntl)
beta.un.var.same <- coef(weight.un.var.same)
sebeta.un.var.same <- robust.cov(weight.un.var.same,m)$se.model
V.un.var.same <- getVarCov(weight.un.var.same)
Gamma.un.var.same <- cov2cor(V.un.var.same)

## > cbind(beta.un.var.same,sebeta.un.var.same)
##              beta.un.var.same sebeta.un.var.same
## diet1             244.4472924         3.64540345
## diet2             248.9026647         4.01703895
## diet3             241.4519884         3.44820645
## diet1:month        -2.0896665         0.60051727
## diet2:month       -11.3738827         0.66173780
## diet3:month        -6.0134582         0.56803247
## diet1:month2        0.2105807         0.05141391
## diet2:month2        0.3206737         0.05665537
## diet3:month2        0.2816450         0.04863269

## > V.un.var.same
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 487.99 337.85 397.50 401.00 393.94
## [2,] 337.85 487.99 349.93 303.03 339.34
## [3,] 397.50 349.93 487.99 353.10 384.96
## [4,] 401.00 303.03 353.10 487.99 341.24
## [5,] 393.94 339.34 384.96 341.24 487.99
##   Standard Deviations: 22.091 22.091 22.091 22.091 22.091 
## > Gamma.un.var.same
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]
## [1,] 1.00000 0.69232 0.81457 0.82172 0.80726
## [2,] 0.69232 1.00000 0.71708 0.62098 0.69537
## [3,] 0.81457 0.71708 1.00000 0.72356 0.78887
## [4,] 0.82172 0.62098 0.72356 1.00000 0.69928
## [5,] 0.80726 0.69537 0.78887 0.69928 1.00000

#  quadratic means, compound symmetric correlation same for all
#  groups, different variances for each group changing over time;
#  the correlation matrix is the same for all groups because gls()
#  does not support it being different, so this is not comparable 
#  the fit with different heterogeneous CS for each group in SAS

weight.csh <- gls(weight ~ -1 + diet + month:diet + month2:diet,data=thedat,
              correlation=corCompSymm(form = ~ 1 | as.factor(id)),
              weights = varIdent(form= ~ 1 | time*diet),method="ML")
beta.csh <- coef(weight.csh)
sebeta.csh <- robust.cov(weight.csh,m)$se.model
V.csh.program1 <- getVarCov(weight.csh,individual=1)  
V.csh.program2 <- getVarCov(weight.csh,individual=35)
V.csh.program3 <- getVarCov(weight.csh,individual=63)
Gamma.csh <- cov2cor(V.csh.program1)

## > cbind(beta.csh,sebeta.csh)
##                 beta.csh sebeta.csh
## diet1        247.3314054 2.93620679
## diet2        246.1927783 3.65215576
## diet3        249.0163366 3.23059808
## diet1:month   -2.4632420 0.61337306
## diet2:month  -12.2084991 0.80283201
## diet3:month   -6.1919514 0.65699183
## diet1:month2   0.2390529 0.04975678
## diet2:month2   0.3880185 0.06627260
## diet3:month2   0.3041219 0.05445448
## > V.csh.program1
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 311.87 302.94 247.96 305.73 263.50
## [2,] 302.94 523.52 321.27 396.11 341.40
## [3,] 247.96 321.27 350.74 324.22 279.44
## [4,] 305.73 396.11 324.22 533.19 344.54
## [5,] 263.50 341.40 279.44 344.54 396.07
##   Standard Deviations: 17.66 22.881 18.728 23.091 19.902 
## > V.csh.program2
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 421.22 316.85 416.22 336.93 390.06
## [2,] 316.85 424.02 417.60 338.05 391.35
## [3,] 416.22 417.60 731.69 444.06 514.08
## [4,] 336.93 338.05 444.06 479.45 416.15
## [5,] 390.06 391.35 514.08 416.15 642.58
##   Standard Deviations: 20.524 20.592 27.05 21.896 25.349 
## > V.csh.program3
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 401.27 392.41 337.71 333.16 358.51
## [2,] 392.41 682.68 440.48 434.56 467.62
## [3,] 337.71 440.48 505.62 373.98 402.44
## [4,] 333.16 434.56 373.98 492.11 397.02
## [5,] 358.51 467.62 402.44 397.02 569.84
##   Standard Deviations: 20.032 26.128 22.486 22.184 23.871 
## > Gamma.csh
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]
## [1,] 1.00000 0.74973 0.74973 0.74973 0.74973
## [2,] 0.74973 1.00000 0.74973 0.74973 0.74973
## [3,] 0.74973 0.74973 1.00000 0.74973 0.74973
## [4,] 0.74973 0.74973 0.74973 1.00000 0.74973
## [5,] 0.74973 0.74973 0.74973 0.74973 1.00000

#  quadratic means, compound symmetric correlation same for all
#  groups, same variances for each group changing over time - this
#  is directly comparable to the analysis in SAS

weight.csh.same <- gls(weight ~ -1 + diet + month:diet + month2:diet,data=thedat,
              correlation=corCompSymm(form = ~ 1 | as.factor(id)),
              weights = varIdent(form= ~ 1 | time),method="ML")
beta.csh.same <- coef(weight.csh.same)
sebeta.csh.same <- robust.cov(weight.csh.same,m)$se.model
V.csh.same <- getVarCov(weight.csh.same)
Gamma.csh.same <- cov2cor(V.csh.same)

## > cbind(beta.csh.same,sebeta.csh.same)
##              beta.csh.same sebeta.csh.same
## diet1          246.5712713      3.29897532
## diet2          251.3401037      3.63529374
## diet3          244.5922266      3.12051825
## diet1:month     -2.5149483      0.69456414
## diet2:month    -11.4433928      0.76537240
## diet3:month     -6.5816495      0.65699190
## diet1:month2     0.2503531      0.05660195
## diet2:month2     0.3357998      0.06237232
## diet3:month2     0.3360307      0.05354009

## > V.csh.same
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 374.56 317.73 310.61 328.50 318.91
## [2,] 317.73 499.34 358.63 379.29 368.22
## [3,] 310.61 358.63 477.20 370.78 359.96
## [4,] 328.50 379.29 370.78 533.76 380.70
## [5,] 318.91 368.22 359.96 380.70 503.06
##   Standard Deviations: 19.354 22.346 21.845 23.103 22.429 
## > Gamma.csh.same
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]
## [1,] 1.00000 0.73468 0.73468 0.73468 0.73468
## [2,] 0.73468 1.00000 0.73468 0.73468 0.73468
## [3,] 0.73468 0.73468 1.00000 0.73468 0.73468
## [4,] 0.73468 0.73468 0.73468 1.00000 0.73468
## [5,] 0.73468 0.73468 0.73468 0.73468 1.00000

#  quadratic means, compound symmetric correlation same for all
#  groups, different variances for each group constant over time

weight.cs <- gls(weight ~ -1 + diet + month:diet + month2:diet,data=thedat,
              correlation=corCompSymm(form = ~ 1 | as.factor(id)),
              weights = varIdent(form= ~ 1 | diet),method="ML")
beta.cs <- coef(weight.cs)
sebeta.cs <- robust.cov(weight.cs,m)$se.model
V.cs.program1 <- getVarCov(weight.cs,individual=1)  
V.cs.program2 <- getVarCov(weight.cs,individual=35)
V.cs.program3 <- getVarCov(weight.cs,individual=63)
Gamma.cs <- cov2cor(V.cs.program1)

## > cbind(beta.cs,sebeta.cs)
##                  beta.cs  sebeta.cs
## diet1        245.0836134 3.38855267
## diet2        250.4356122 4.23197502
## diet3        242.5089474 3.55848529
## diet1:month   -2.2769188 0.66904180
## diet2:month  -11.5963605 0.83556859
## diet3:month   -6.2914035 0.70259360
## diet1:month2   0.2355275 0.05346350
## diet2:month2   0.3465420 0.06677075
## diet3:month2   0.3181287 0.05614464

## > V.cs.program1
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 402.99 292.79 292.79 292.79 292.79
## [2,] 292.79 402.99 292.79 292.79 292.79
## [3,] 292.79 292.79 402.99 292.79 292.79
## [4,] 292.79 292.79 292.79 402.99 292.79
## [5,] 292.79 292.79 292.79 292.79 402.99
##   Standard Deviations: 20.075 20.075 20.075 20.075 20.075 
## > V.cs.program2
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 517.65 376.09 376.09 376.09 376.09
## [2,] 376.09 517.65 376.09 376.09 376.09
## [3,] 376.09 376.09 517.65 376.09 376.09
## [4,] 376.09 376.09 376.09 517.65 376.09
## [5,] 376.09 376.09 376.09 376.09 517.65
##   Standard Deviations: 22.752 22.752 22.752 22.752 22.752 
## > V.cs.program3
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 496.71 360.88 360.88 360.88 360.88
## [2,] 360.88 496.71 360.88 360.88 360.88
## [3,] 360.88 360.88 496.71 360.88 360.88
## [4,] 360.88 360.88 360.88 496.71 360.88
## [5,] 360.88 360.88 360.88 360.88 496.71
##   Standard Deviations: 22.287 22.287 22.287 22.287 22.287 

## > Gamma.cs
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]
## [1,] 1.00000 0.72653 0.72653 0.72653 0.72653
## [2,] 0.72653 1.00000 0.72653 0.72653 0.72653
## [3,] 0.72653 0.72653 1.00000 0.72653 0.72653
## [4,] 0.72653 0.72653 0.72653 1.00000 0.72653
## [5,] 0.72653 0.72653 0.72653 0.72653 1.00000

#  quadratic means, compound symmetric correlation same for all
#  groups, same variance for all groups, constant over time

weight.cs.same <- gls(weight ~ -1 + diet + month:diet + month2:diet,data=thedat,
              correlation=corCompSymm(form = ~ 1 | as.factor(id)),
              method="ML")
beta.cs.same <- coef(weight.cs.same)
sebeta.cs.same <- robust.cov(weight.cs.same,m)$se.model
V.cs.same <- getVarCov(weight.cs.same)  
Gamma.cs.same <- cov2cor(V.cs.same)

## > cbind(beta.cs.same,sebeta.cs.same)
##              beta.cs.same sebeta.cs.same
## diet1         245.0836134     3.67703725
## diet2         250.4356122     4.05189770
## diet3         242.5089474     3.47812903
## diet1:month    -2.2769188     0.72210204
## diet2:month   -11.5963605     0.79571769
## diet3:month    -6.2914035     0.68304015
## diet1:month2    0.2355275     0.05770357
## diet2:month2    0.3465420     0.06358624
## diet3:month2    0.3181287     0.05458211

## > V.cs.same
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 474.37 345.99 345.99 345.99 345.99
## [2,] 345.99 474.37 345.99 345.99 345.99
## [3,] 345.99 345.99 474.37 345.99 345.99
## [4,] 345.99 345.99 345.99 474.37 345.99
## [5,] 345.99 345.99 345.99 345.99 474.37
##   Standard Deviations: 21.78 21.78 21.78 21.78 21.78 

## > Gamma.cs.same
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]
## [1,] 1.00000 0.72937 0.72937 0.72937 0.72937
## [2,] 0.72937 1.00000 0.72937 0.72937 0.72937
## [3,] 0.72937 0.72937 1.00000 0.72937 0.72937
## [4,] 0.72937 0.72937 0.72937 1.00000 0.72937
## [5,] 0.72937 0.72937 0.72937 0.72937 1.00000

#  Get AIC/BIC values

## > anova(weight.un,weight.un.same,weight.un.var,weight.un.var.same,
## weight.csh,weight.csh.same,weight.cs,weight.cs.same)
##                    Model df      AIC      BIC    logLik   Test   L.Ratio
## weight.un              1 34 4133.842 4277.139 -2032.921                 
## weight.un.same         2 24 4124.350 4225.501 -2038.175 1 vs 2 10.508163
## weight.un.var          3 22 4125.886 4218.607 -2040.943 2 vs 3  5.535175
## weight.un.var.same     4 20 4123.841 4208.134 -2041.921 3 vs 4  1.955797
## weight.csh             5 25 4143.834 4249.199 -2046.917 4 vs 5  9.992756
## weight.csh.same        6 15 4136.056 4199.275 -2053.028 5 vs 6 12.221444
## weight.cs              7 13 4136.806 4191.596 -2055.403 6 vs 7  4.750886
## weight.cs.same         8 11 4135.680 4182.041 -2056.840 7 vs 8  2.873928
##                    p-value
## weight.un                 
## weight.un.same      0.3971
## weight.un.var       0.0628
## weight.un.var.same  0.3761
## weight.csh          0.0754
## weight.csh.same     0.2705
## weight.cs           0.0930
## weight.cs.same      0.2376

#  Based on the information criteria, the model with common
#  unstructured correlation matrix with common variance that is the SAME
#  over time.  This model is NOT possible in proc mixed.  In the SAS program
#  given we cannot fit this model, we use the common unstructured model (whose
#  variances do change over time, which is the model weight.un.same above).

u <- weight.un.var.same

#  We use model-based standard errors here, given the unstructured covariance;
#  you may have used robust standard errors

Sig.model <- robust.cov(u,m)$Sig.model
t.model <- beta.un.var.same/sebeta.un.var.same
df <- nrow(thedat)-length(beta.un.var.same)    #  this is the residual df used
                                               #  by gls; SAS uses a different approximation
p.value <- round(2*pt(-abs(t.model),df),4)

#  Results with using residual degreees of freedom and t-test -
#  slightly different from those in SAS but qualitatively same

## > cbind(beta.un.var.same,sebeta.un.var.same,t.model,p.value)
##              beta.un.var.same sebeta.un.var.same    t.model p.value
## diet1             244.4472924         3.64540345  67.056307   0e+00
## diet2             248.9026647         4.01703895  61.961726   0e+00
## diet3             241.4519884         3.44820645  70.022486   0e+00
## diet1:month        -2.0896665         0.60051727  -3.479778   5e-04
## diet2:month       -11.3738827         0.66173780 -17.187899   0e+00
## diet3:month        -6.0134582         0.56803247 -10.586469   0e+00
## diet1:month2        0.2105807         0.05141391   4.095792   0e+00
## diet2:month2        0.3206737         0.05665537   5.660077   0e+00
## diet3:month2        0.2816450         0.04863269   5.791269   0e+00

#  Function to contstruct F test statistics

f.test <- function(L,beta,Sig,dendf){
    numdf <- nrow(L)
    F <- t(L%*%beta)%*%solve(L%*%Sig%*%t(L))%*%(L%*%beta)/numdf
    p.value <- pf(F,numdf,dendf,lower.tail=FALSE)
    return(list(F=F,p.value=p.value,F.pvalue=c(F,round(p.value,4))))
}
    
#  Test of whether or not all quadratic terms = 0

L <- matrix(c(0,0,0,0,0,0,1,0,0,
              0,0,0,0,0,0,0,1,0,
              0,0,0,0,0,0,0,0,1),3,9,byrow=TRUE)

## > f.test(L,beta.un.var.same,Sig.model,df)$F.pvalue
## [1] 27.45026  0.00000

#  Test of whether quadratic terms differ

L  <- matrix(c(0,0,0,0,0,0,1,-1,0,
               0,0,0,0,0,0,1,0,-1),2,9,byrow=TRUE)

## > f.test(L,beta.un.var.same,Sig.model,df)$F.pvalue
## [1] 1.095084 0.335300

#  Function to get estimate of form L beta and standard error
#  L has one row in this case; alpha is conf level for conf interval

est.L <- function(L,beta,Sig,alpha,df){
    est <- L%*%beta
    se <- sqrt(L%*%Sig%*%t(L))
    crit.val <- abs(qt(alpha/2,df))
    CI <- c(est-crit.val*se,est+crit.val*se)
    result <- matrix(c(est,se,CI),1,4,byrow=TRUE)
    colnames(result) <- c("Estimate","Std Error","Lower CI", "Upper CI")
    return(result)
}

#  Estimates of mean @ month 12 for each diet and rate of change
#  at 8 months

#  mean @ month 12, diet 1

L <- matrix(c(1,0,0,12,0,0,144,0,0),1,9,byrow=TRUE)

## > est.L(L,beta.un.var.same,Sig.model,0.05,df)
##     Estimate Std Error Lower CI Upper CI
## [1,] 249.6949  3.674325 242.4756 256.9143

#  mean @ month 12, diet 2

L <- matrix(c(0,1,0,0,12,0,0,144,0),1,9,byrow=TRUE)

## > est.L(L,beta.un.var.same,Sig.model,0.05,df)
##      Estimate Std Error Lower CI Upper CI
## [1,] 158.5931  4.048909 150.6378 166.5484

#  mean @ month 12, diet 3

L <- matrix(c(0,0,1,0,0,12,0,0,144),1,9,byrow=TRUE)

## >  est.L(L,beta.un.var.same,Sig.model,0.05,df)
##      Estimate Std Error Lower CI Upper CI
## [1,] 209.8474  3.475563 203.0186 216.6762

#  difference in means @ month 12

L  <- matrix(c(1,-1,0,12,-12,0,144,-144,0,
               1,0,-1,12,0,-12,144,0,-144),2,9,byrow=TRUE)

## > f.test(L,beta.un.var.same,Sig.model,df)$F.pvalue
## [1] 138.8577   0.0000

#  rate of change @ month 8, diet 1

L <- matrix(c(0,0,0,1,0,0,16,0,0),1,9,byrow=TRUE)

## > est.L(L,beta.un.var.same,Sig.model,0.05,df)
##      Estimate Std Error  Lower CI Upper CI
## [1,] 1.279625 0.2927269 0.7044728 1.854777

#  rate of change @ month 8, diet 2

L <- matrix(c(0,0,0,0,1,0,0,16,0),1,9,byrow=TRUE)

## > est.L(L,beta.un.var.same,Sig.model,0.05,df)
##       Estimate Std Error  Lower CI  Upper CI
## [1,] -6.243103 0.3225693 -6.876889 -5.609316

#  rate of change @ month 8, diet 3

L <- matrix(c(0,0,0,0,0,1,0,0,16),1,9,byrow=TRUE)

## > est.L(L,beta.un.var.same,Sig.model,0.05,df)
##       Estimate Std Error  Lower CI   Upper CI
## [1,] -1.507138 0.2768919 -2.051177 -0.9630988

#  Fit the model in alternative parametrization -- define program 3
#  (diet 3) to be the reference level to match with SAS default

diet3 <- relevel(diet,ref="3")

weight.un.alt <- gls(weight ~ diet3 + month + month2+ month:diet3 + month2:diet3,data=thedat,
                    correlation=corSymm(form = ~ 1 | as.factor(id)),
                    method="ML",control=weight.cntl)
beta.un.alt <- coef(weight.un.alt)
u <- weight.un.alt
Sig.model.alt <- robust.cov(u,m)$Sig.model
sebeta.un.alt <- robust.cov(u,m)$se.model
t.alt <- beta.un.alt/sebeta.un.alt
p.value.alt <- round(2*pt(-abs(t.alt),df),4)

## > cbind(beta.un.alt,sebeta.un.alt,t.alt,p.value.alt)
##                beta.un.alt sebeta.un.alt       t.alt p.value.alt
## (Intercept)   241.45198843    3.44820645  70.0224862      0.0000
## diet31          2.99530393    5.01787744   0.5969265      0.5508
## diet32          7.45067630    5.29402773   1.4073739      0.1599
## month          -6.01345817    0.56803247 -10.5864691      0.0000
## month2          0.28164501    0.04863269   5.7912690      0.0000
## diet31:month    3.92379165    0.82660866   4.7468552      0.0000
## diet32:month   -5.36042450    0.87209965  -6.1465734      0.0000
## diet31:month2  -0.07106431    0.07077096  -1.0041450      0.3158
## diet32:month2   0.03902874    0.07466572   0.5227129      0.6014

#  test of whether or not quadratic effects are the same 

L <- matrix(c(0,0,0,0,0,0,0,1,0,
              0,0,0,0,0,0,0,0,1),2,9,byrow=TRUE)

## > f.test(L,beta.un.alt,Sig.model.alt,df)$F.pvalue
## [1] 1.095084 0.335300

#  There does not seem to be strong evidence that the quadratic effect is different 
#  across the three groups, so it would be possible to move to a reduced model with 
#  common quadratic term and do all of the above analyses with this model

weight.un.var.same.quad <- gls(weight ~ -1 + diet + month:diet + month2,data=thedat,
              correlation=corSymm(form = ~ 1 | as.factor(id)),
              method="ML",control=weight.cntl)
beta.un.var.same.quad <- coef(weight.un.var.same.quad)
sebeta.un.var.same.quad <- robust.cov(weight.un.var.same.quad,m)$se.model
V.un.var.same.quad <- getVarCov(weight.un.var.same.quad)
Gamma.un.var.same.quad <- cov2cor(V.un.var.same.quad)

u <- weight.un.var.same.quad

Sig.model <- robust.cov(u,m)$Sig.model
t.model <- beta.un.var.same.quad/sebeta.un.var.same.quad
df <- nrow(thedat)-length(beta.un.var.same.quad)    #  this is the residual df used
                                               #  by gls; SAS uses a different approximation
p.value <- round(2*pt(-abs(t.model),df),4)

## > cbind(beta.un.var.same.quad,sebeta.un.var.same.quad,t.model,p.value)
##             beta.un.var.same.quad sebeta.un.var.same.quad    t.model p.value
## diet1                 245.5174759              3.58971584  68.394683       0
## diet2                 248.0277718              3.94871730  62.812238       0
## diet3                 241.3234202              3.39950660  70.987778       0
## month2                  0.2730455              0.03014172   9.058725       0
## diet1:month            -2.7969784              0.37585125  -7.441716       0
## diet2:month           -10.8321209              0.38323995 -28.264592       0
## diet3:month            -5.9214369              0.37216799 -15.910656       0






