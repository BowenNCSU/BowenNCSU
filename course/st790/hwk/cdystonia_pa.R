#########################################################################
#
#   Cdystonia data - three weightloss programs
#
#   0 = Placebo
#   5000 = 5000 units botulinum toxin type B (BotB)
#   10000 = 10000 units BotB
#
#   Homework 2 population averaged analyses with spline model
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

#  Read in the data -- they are in the "long" form of one record per observation

thedat <- read.table("cdystonia.dat",row.names=NULL)
colnames(thedat) <- c("id","week","grp","age","gender","twstrs")

#  Total number of individuals

m <- length(unique(thedat$id))

week2 <- thedat$week^2
weekplus <- (thedat$week-4)*(thedat$week-4 > 0)
dose <- as.factor(thedat$grp)

# time variable for specifying covariance matrix

thedat$time <- as.numeric(factor(thedat$week))

#  Fit different covariance structures -- first fit a common
#  unstructured correlation with variances that change
#  over time to get a sense of the structure; as we know, we cannot  
#  fit separate separate covariance structures in gls().  We could 
#  fit three separate models by group, but I don't do that here; we'd
#  have to add up the objective functions and compute AIC and BIC ourselves.
#  But I already fit models with separate covariance structures in SAS.

#  As in the SAS program, we maintain separate intercepts for each
#  program; you might have decided based on randomization and the
#  results to collapse to a single intercept for all groups. We fit a
#  series of models involving unstructured and compound symmetric
#  overall covariance structures and then compare via AIC and BIC.

#  increase number of iterations

twstrs.cntl <- glsControl(maxIter=200,msMaxIter=200,tolerance=1e-5,msTol=1e-6)

# As a first shot, fit the most general model possible with gls() with
# a intercept, 1st stage and 2nd stage slope for each group, unstructured correlation
# with different variances for each group different at each time point
# Of course, in R we are stuck with having the correlation matrix being the
# same for each group, so we can't compare directly to SAS

#  For some reason getVarCov does not work with this.  I will
#  calculate the covariance matrix 

twstrs.un <- gls(twstrs ~ -1 + dose + week:dose + weekplus:dose,data=thedat,
              correlation=corSymm(form = ~ time | as.factor(id)),
              weights = varIdent(form= ~ 1 | week*dose),method="ML",control=twstrs.cntl)
beta.un <- coef(twstrs.un)
sebeta.un <- robust.cov(twstrs.un,m)$se.model
V.un.dose1 <- getVarCov(twstrs.un, individual=1)  #  or corMatrix(twstrs.un$modelStruct$corStruct)[[1]]
V.un.dose2 <- getVarCov(twstrs.un, individual=6)
V.un.dose3 <- getVarCov(twstrs.un, individual=8)
Gamma.un <- cov2cor(V.un.dose1)

## > cbind(beta.un,sebeta.un)
##                      beta.un sebeta.un
## dose0              43.147647 1.5212205
## dose5000           46.419266 1.4817168
## dose10000          47.267741 1.6515984
## dose0:week         -1.097173 0.3154363
## dose5000:week      -2.135302 0.3106792
## dose10000:week     -2.808894 0.3370663
## dose0:weekplus      1.433937 0.3879561
## dose5000:weekplus   2.777506 0.4018853
## dose10000:weekplus  3.943620 0.4233584
## > V.un.dose1
## Marginal variance covariance matrix
##        [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
## [1,] 80.229  82.009  90.169  81.749  69.498  67.115
## [2,] 82.009 166.950 131.570 108.250  92.049  93.057
## [3,] 90.169 131.570 163.730 129.940 103.180 106.670
## [4,] 81.749 108.250 129.940 135.870 102.790  93.060
## [5,] 69.498  92.049 103.180 102.790 105.020  96.418
## [6,] 67.115  93.057 106.670  93.060  96.418 124.430
##   Standard Deviations: 8.9571 12.921 12.796 11.656 10.248 11.155 
## > V.un.dose2
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]   [,5]    [,6]
## [1,] 112.780 110.930 113.480 121.640 102.67  67.779
## [2,] 110.930 217.290 159.340 154.990 130.84  90.428
## [3,] 113.480 159.340 184.500 173.120 136.46  96.447
## [4,] 121.640 154.990 173.120 214.020 160.74  99.484
## [5,] 102.670 130.840 136.460 160.740 163.03 102.330
## [6,]  67.779  90.428  96.447  99.484 102.33  90.276
##   Standard Deviations: 10.62 14.741 13.583 14.629 12.768 9.5014 
## > V.un.dose3
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]   [,4]    [,5]    [,6]
## [1,]  88.217  81.797  92.903 107.32  97.572  95.489
## [2,]  81.797 151.050 122.970 128.90 117.230 120.100
## [3,]  92.903 122.970 158.070 159.84 135.730 142.200
## [4,] 107.320 128.900 159.840 212.95 172.290 158.080
## [5,]  97.572 117.230 135.730 172.29 188.260 175.160
## [6,]  95.489 120.100 142.200 158.08 175.160 229.070
##   Standard Deviations: 9.3924 12.29 12.573 14.593 13.721 15.135 
## > Gamma.un
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
## [1,] 1.00000 0.70859 0.78673 0.78298 0.75714 0.67174
## [2,] 0.70859 1.00000 0.79578 0.71870 0.69517 0.64564
## [3,] 0.78673 0.79578 1.00000 0.87119 0.78683 0.74731
## [4,] 0.78298 0.71870 0.87119 1.00000 0.86050 0.71572
## [5,] 0.75714 0.69517 0.78683 0.86050 1.00000 0.84347
## [6,] 0.67174 0.64564 0.74731 0.71572 0.84347 1.00000
##   Standard Deviations: 1 1 1 1 1 1 

#  unstructured correlation, variances same for all groups
#  different at each time point; this is directly comparable to SAS with common
#  unstructured covariance matrix

twstrs.un.same <- gls(twstrs ~ -1 + dose + week:dose + weekplus:dose,data=thedat,
              correlation=corSymm(form = ~ time | as.factor(id)),
              weights = varIdent(form= ~ 1 | week),method="ML",control=twstrs.cntl)
beta.un.same <- coef(twstrs.un.same)
sebeta.un.same <- robust.cov(twstrs.un.same,m)$se.model
V.un.same <- getVarCov(twstrs.un.same, individual=1)  
Gamma.un.same <- cov2cor(V.un.same)

##  These results are identical to PROC MIXED with common unstructured
##  covariance structure

## > cbind(beta.un.same,sebeta.un.same)
##                    beta.un.same sebeta.un.same
## dose0                 43.261296      1.5699751
## dose5000              46.072302      1.5715273
## dose10000             46.667938      1.5491940
## dose0:week            -1.059819      0.3307144
## dose5000:week         -2.147837      0.3327841
## dose10000:week        -2.865269      0.3237450
## dose0:weekplus         1.379673      0.4195902
## dose5000:weekplus      2.787521      0.4258288
## dose10000:weekplus     4.027760      0.4111935
## > V.un.same
## Marginal variance covariance matrix
##        [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
## [1,] 91.416  87.037  95.937  98.425  84.218  73.781
## [2,] 87.037 172.910 134.010 123.960 106.320  95.479
## [3,] 95.937 134.010 167.030 148.810 118.440 109.110
## [4,] 98.425 123.960 148.810 177.400 133.440 108.880
## [5,] 84.218 106.320 118.440 133.440 142.680 117.990
## [6,] 73.781  95.479 109.110 108.880 117.990 141.260
##   Standard Deviations: 9.5612 13.15 12.924 13.319 11.945 11.885 
## > Gamma.un.same
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
## [1,] 1.00000 0.69227 0.77639 0.77288 0.73743 0.64928
## [2,] 0.69227 1.00000 0.78854 0.70774 0.67690 0.61093
## [3,] 0.77639 0.78854 1.00000 0.86448 0.76722 0.71036
## [4,] 0.77288 0.70774 0.86448 1.00000 0.83872 0.68779
## [5,] 0.73743 0.67690 0.76722 0.83872 1.00000 0.83116
## [6,] 0.64928 0.61093 0.71036 0.68779 0.83116 1.00000

#  With gls(), we can also specify an unstructured correlation matrix
#  with the same variance at all time points, which can be either
#  different by group or the same for all groups.  We can't fit this in SAS.

twstrs.un.var <- gls(twstrs ~ -1 + dose + week:dose + weekplus:dose,data=thedat,
              correlation=corSymm(form = ~ 1 | as.factor(id)),
              weights = varIdent(form= ~ 1 | dose),method="ML",control=twstrs.cntl)
beta.un.var <- coef(twstrs.un.var)
sebeta.un.var <- robust.cov(twstrs.un.var,m)$se.model
V.un.var.dose1 <- getVarCov(twstrs.un.var,individual=1)
V.un.var.dose2 <- getVarCov(twstrs.un.var,individual=6)
V.un.var.dose3 <- getVarCov(twstrs.un.var,individual=8)
Gamma.un.var <- cov2cor(V.un.var.dose1)

## > cbind(beta.un.var,sebeta.un.var)
##                   beta.un.var sebeta.un.var
## dose0                43.045346     1.9366445
## dose5000             43.909837     1.7459441
## dose10000            44.792524     1.8846401
## dose0:week           -1.043694     0.3602305
## dose5000:week        -2.098770     0.3283415
## dose10000:week       -2.780761     0.3490987
## dose0:weekplus        1.423848     0.4505271
## dose5000:weekplus     2.820712     0.4097429
## dose10000:weekplus    3.927834     0.4344264

## ## > V.un.var.dose1
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
## [1,] 127.020  82.410  92.864  96.733  96.725  80.203
## [2,]  82.410 127.020  95.061  86.544  88.680  72.362
## [3,]  92.864  95.061 127.020 109.900  98.180  86.481
## [4,]  96.733  86.544 109.900 127.020 104.590  84.539
## [5,]  96.725  88.680  98.180 104.590 127.020 108.390
## [6,]  80.203  72.362  86.481  84.539 108.390 127.020
##   Standard Deviations: 11.27 11.27 11.27 11.27 11.27 11.27 
## > V.un.var.dose2
## Marginal variance covariance matrix
##         [,1]    [,2]   [,3]   [,4]   [,5]    [,6]
## [1,] 153.570  99.639 112.28 116.96 116.95  96.971
## [2,]  99.639 153.570 114.94 104.64 107.22  87.491
## [3,] 112.280 114.940 153.57 132.87 118.71 104.560
## [4,] 116.960 104.640 132.87 153.57 126.45 102.210
## [5,] 116.950 107.220 118.71 126.45 153.57 131.050
## [6,]  96.971  87.491 104.56 102.21 131.05 153.570
##   Standard Deviations: 12.392 12.392 12.392 12.392 12.392 12.392 
## > V.un.var.dose3
## Marginal variance covariance matrix
##         [,1]    [,2]   [,3]   [,4]   [,5]    [,6]
## [1,] 157.250 102.030 114.97 119.76 119.75  99.293
## [2,] 102.030 157.250 117.69 107.14 109.79  89.586
## [3,] 114.970 117.690 157.25 136.06 121.55 107.070
## [4,] 119.760 107.140 136.06 157.25 129.48 104.660
## [5,] 119.750 109.790 121.55 129.48 157.25 134.190
## [6,]  99.293  89.586 107.07 104.66 134.19 157.250
##   Standard Deviations: 12.54 12.54 12.54 12.54 12.54 12.54 
## > Gamma.un.var
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
## [1,] 1.00000 0.64882 0.73112 0.76158 0.76152 0.63144
## [2,] 0.64882 1.00000 0.74842 0.68136 0.69818 0.56971
## [3,] 0.73112 0.74842 1.00000 0.86523 0.77297 0.68087
## [4,] 0.76158 0.68136 0.86523 1.00000 0.82342 0.66558
## [5,] 0.76152 0.69818 0.77297 0.82342 1.00000 0.85337
## [6,] 0.63144 0.56971 0.68087 0.66558 0.85337 1.00000

#  variance common across time same for all groups

twstrs.un.var.same <- gls(twstrs ~ -1 + dose + week:dose + weekplus:dose,data=thedat,
              correlation=corSymm(form = ~ time | as.factor(id)),
              method="ML",control=twstrs.cntl)
beta.un.var.same <- coef(twstrs.un.var.same)
sebeta.un.var.same <- robust.cov(twstrs.un.var.same,m)$se.model
V.un.var.same <- getVarCov(twstrs.un.var.same)
Gamma.un.var.same <- cov2cor(V.un.var.same)

## > cbind(beta.un.var.same,sebeta.un.var.same)
##                    beta.un.var.same sebeta.un.var.same
## dose0                    43.0329338          1.8732062
## dose5000                 43.9758733          1.8716612
## dose10000                44.5557565          1.8419694
## dose0:week               -0.9999867          0.3289368
## dose5000:week            -2.0065192          0.3316030
## dose10000:week           -2.7240744          0.3221556
## dose0:weekplus            1.3204992          0.4091524
## dose5000:weekplus         2.6567987          0.4159986
## dose10000:weekplus        3.8968798          0.4009856
## > V.un.var.same
## Marginal variance covariance matrix
##         [,1]    [,2]   [,3]    [,4]    [,5]    [,6]
## [1,] 143.740  92.171 108.55 106.000 106.740  94.922
## [2,]  92.171 143.740 110.24  97.712  94.153  84.950
## [3,] 108.550 110.240 143.74 121.410 108.080 101.080
## [4,] 106.000  97.712 121.41 143.740 118.340  96.842
## [5,] 106.740  94.153 108.08 118.340 143.740 120.020
## [6,]  94.922  84.950 101.08  96.842 120.020 143.740
##   Standard Deviations: 11.989 11.989 11.989 11.989 11.989 11.989 
## > Gamma.un.var.same
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
## [1,] 1.00000 0.64123 0.75519 0.73742 0.74257 0.66037
## [2,] 0.64123 1.00000 0.76692 0.67978 0.65502 0.59099
## [3,] 0.75519 0.76692 1.00000 0.84461 0.75192 0.70320
## [4,] 0.73742 0.67978 0.84461 1.00000 0.82331 0.67372
## [5,] 0.74257 0.65502 0.75192 0.82331 1.00000 0.83494
## [6,] 0.66037 0.59099 0.70320 0.67372 0.83494 1.00000

#  compound symmetric correlation same for all
#  groups, different variances for each group changing over time;
#  the correlation matrix is the same for all groups because gls()
#  does not support it being different, so this is not comparable 
#  the fit with different heterogeneous CS for each group in SAS

twstrs.csh <- gls(twstrs ~ -1 + dose + week:dose + weekplus:dose,data=thedat,
              correlation=corCompSymm(form = ~ 1 | as.factor(id)),
              weights = varIdent(form= ~ 1 | week*dose),method="ML")
beta.csh <- coef(twstrs.csh)
sebeta.csh <- robust.cov(twstrs.csh,m)$se.model
V.csh.dose1 <- getVarCov(twstrs.csh,individual=1)  
V.csh.dose2 <- getVarCov(twstrs.csh,individual=6)
V.csh.dose3 <- getVarCov(twstrs.csh,individual=8)
Gamma.csh <- cov2cor(V.csh.dose1)

## > cbind(beta.csh,sebeta.csh)
##                     beta.csh sebeta.csh
## dose0              43.520997  1.5859796
## dose5000           47.217008  1.4723333
## dose10000          47.810593  1.7374654
## dose0:week         -1.035993  0.3085061
## dose5000:week      -2.267725  0.3106958
## dose10000:week     -2.966964  0.3293355
## dose0:weekplus      1.362944  0.3888334
## dose5000:weekplus   3.012887  0.4011724
## dose10000:weekplus  4.210435  0.4160822

## > V.csh.dose1
## Marginal variance covariance matrix
##        [,1]   [,2]    [,3]    [,4]    [,5]    [,6]
## [1,] 80.860  91.23  88.118  83.656  70.726  76.798
## [2,] 91.230 172.96 128.880 122.350 103.440 112.320
## [3,] 88.118 128.88 161.360 118.180  99.911 108.490
## [4,] 83.656 122.35 118.180 145.430  94.852 102.990
## [5,] 70.726 103.44  99.911  94.852 103.950  87.076
## [6,] 76.798 112.32 108.490 102.990  87.076 122.570
##   Standard Deviations: 8.9922 13.152 12.703 12.06 10.196 11.071 
## > V.csh.dose2
## Marginal variance covariance matrix
##         [,1]   [,2]   [,3]   [,4]   [,5]    [,6]
## [1,] 118.020 130.48 111.06 118.07 108.37  86.122
## [2,] 130.480 242.39 159.17 169.21 155.31 123.430
## [3,] 111.060 159.17 175.64 144.04 132.21 105.060
## [4,] 118.070 169.21 144.04 198.49 140.55 111.690
## [5,] 108.370 155.31 132.21 140.55 167.23 102.520
## [6,]  86.122 123.43 105.06 111.69 102.52 105.610
##   Standard Deviations: 10.864 15.569 13.253 14.089 12.932 10.277 
## > V.csh.dose3
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]   [,4]    [,5]   [,6]
## [1,]  93.097  94.314  90.981 104.86  97.987 118.83
## [2,]  94.314 160.560 119.480 137.71 128.680 156.05
## [3,]  90.981 119.480 149.410 132.85 124.130 150.54
## [4,] 104.860 137.710 132.850 198.49 143.080 173.51
## [5,]  97.987 128.680 124.130 143.08 173.310 162.13
## [6,] 118.830 156.050 150.540 173.51 162.130 254.88
##   Standard Deviations: 9.6487 12.671 12.223 14.089 13.165 15.965 
## > Gamma.csh
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
## [1,] 1.00000 0.77143 0.77143 0.77143 0.77143 0.77143
## [2,] 0.77143 1.00000 0.77143 0.77143 0.77143 0.77143
## [3,] 0.77143 0.77143 1.00000 0.77143 0.77143 0.77143
## [4,] 0.77143 0.77143 0.77143 1.00000 0.77143 0.77143
## [5,] 0.77143 0.77143 0.77143 0.77143 1.00000 0.77143
## [6,] 0.77143 0.77143 0.77143 0.77143 0.77143 1.00000

#  compound symmetric correlation same for all groups, 
#  same variances for each group changing over time - this
#  is directly comparable to the analysis in SAS and is the one we chose
#  from SAS results 

twstrs.csh.same <- gls(twstrs ~ -1 + dose + week:dose + weekplus:dose,data=thedat,
              correlation=corCompSymm(form = ~ 1 | as.factor(id)),
              weights = varIdent(form= ~ 1 | time),method="ML")
beta.csh.same <- coef(twstrs.csh.same)
sebeta.csh.same <- robust.cov(twstrs.csh.same,m)$se.model
V.csh.same <- getVarCov(twstrs.csh.same)
Gamma.csh.same <- cov2cor(V.csh.same)

## > cbind(beta.csh.same,sebeta.csh.same)
##                    beta.csh.same sebeta.csh.same
## dose0                  43.627249       1.6037983
## dose5000               46.884982       1.6041651
## dose10000              47.309093       1.5820023
## dose0:week             -1.002853       0.3251695
## dose5000:week          -2.300083       0.3293273
## dose10000:week         -3.040340       0.3196097
## dose0:weekplus          1.292958       0.4118029
## dose5000:weekplus       3.067702       0.4198651
## dose10000:weekplus      4.352583       0.4045895

## > V.csh.same
## Marginal variance covariance matrix
##        [,1]    [,2]    [,3]   [,4]    [,5]    [,6]
## [1,] 93.816  98.892  92.562  96.29  86.455  90.318
## [2,] 98.892 183.610 129.490 134.71 120.950 126.350
## [3,] 92.562 129.490 160.850 126.08 113.210 118.260
## [4,] 96.290 134.710 126.080 174.07 117.770 123.030
## [5,] 86.455 120.950 113.210 117.77 140.330 110.460
## [6,] 90.318 126.350 118.260 123.03 110.460 153.150
##   Standard Deviations: 9.6859 13.55 12.683 13.194 11.846 12.375

## > Gamma.csh.same
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
## [1,] 1.00000 0.75349 0.75349 0.75349 0.75349 0.75349
## [2,] 0.75349 1.00000 0.75349 0.75349 0.75349 0.75349
## [3,] 0.75349 0.75349 1.00000 0.75349 0.75349 0.75349
## [4,] 0.75349 0.75349 0.75349 1.00000 0.75349 0.75349
## [5,] 0.75349 0.75349 0.75349 0.75349 1.00000 0.75349
## [6,] 0.75349 0.75349 0.75349 0.75349 0.75349 1.00000

#  compound symmetric correlation same for all
#  groups, different variances for each group constant over time

twstrs.cs <- gls(twstrs ~ -1 + dose + week:dose + weekplus:dose,data=thedat,
              correlation=corCompSymm(form = ~ 1 | as.factor(id)),
              weights = varIdent(form= ~ 1 | dose),method="ML")
beta.cs <- coef(twstrs.cs)
sebeta.cs <- robust.cov(twstrs.cs,m)$se.model
V.cs.dose1 <- getVarCov(twstrs.cs,individual=1)  
V.cs.dose2 <- getVarCov(twstrs.cs,individual=6)
V.cs.dose3 <- getVarCov(twstrs.cs,individual=8)
Gamma.cs <- cov2cor(V.cs.dose1)

## > cbind(beta.cs,sebeta.cs)
##                       beta.cs sebeta.cs
## dose0              43.0226873 2.0784508
## dose5000           44.8140382 1.8688352
## dose10000          45.2910081 2.0437160
## dose0:week         -0.9556101 0.3517708
## dose5000:week      -2.2501421 0.3204879
## dose10000:week     -2.9780574 0.3455527
## dose0:weekplus      1.2519602 0.4345738
## dose5000:weekplus   3.0354231 0.3985165
## dose10000:weekplus  4.3073372 0.4266637

## > V.cs.dose1
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
## [1,] 131.160  98.087  98.087  98.087  98.087  98.087
## [2,]  98.087 131.160  98.087  98.087  98.087  98.087
## [3,]  98.087  98.087 131.160  98.087  98.087  98.087
## [4,]  98.087  98.087  98.087 131.160  98.087  98.087
## [5,]  98.087  98.087  98.087  98.087 131.160  98.087
## [6,]  98.087  98.087  98.087  98.087  98.087 131.160
##   Standard Deviations: 11.452 11.452 11.452 11.452 11.452 11.452 
## > V.cs.dose2
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]   [,6]
## [1,] 161.44 120.74 120.74 120.74 120.74 120.74
## [2,] 120.74 161.44 120.74 120.74 120.74 120.74
## [3,] 120.74 120.74 161.44 120.74 120.74 120.74
## [4,] 120.74 120.74 120.74 161.44 120.74 120.74
## [5,] 120.74 120.74 120.74 120.74 161.44 120.74
## [6,] 120.74 120.74 120.74 120.74 120.74 161.44
##   Standard Deviations: 12.706 12.706 12.706 12.706 12.706 12.706 
## > V.cs.dose3
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]   [,6]
## [1,] 162.14 121.26 121.26 121.26 121.26 121.26
## [2,] 121.26 162.14 121.26 121.26 121.26 121.26
## [3,] 121.26 121.26 162.14 121.26 121.26 121.26
## [4,] 121.26 121.26 121.26 162.14 121.26 121.26
## [5,] 121.26 121.26 121.26 121.26 162.14 121.26
## [6,] 121.26 121.26 121.26 121.26 121.26 162.14
##   Standard Deviations: 12.733 12.733 12.733 12.733 12.733 12.733 
## > Gamma.cs
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
## [1,] 1.00000 0.74787 0.74787 0.74787 0.74787 0.74787
## [2,] 0.74787 1.00000 0.74787 0.74787 0.74787 0.74787
## [3,] 0.74787 0.74787 1.00000 0.74787 0.74787 0.74787
## [4,] 0.74787 0.74787 0.74787 1.00000 0.74787 0.74787
## [5,] 0.74787 0.74787 0.74787 0.74787 1.00000 0.74787
## [6,] 0.74787 0.74787 0.74787 0.74787 0.74787 1.00000

#  compound symmetric correlation same for all
#  groups, same variance for all groups, constant over time

twstrs.cs.same <- gls(twstrs ~ -1 + dose + week:dose + weekplus:dose,data=thedat,
              correlation=corCompSymm(form = ~ 1 | as.factor(id)),
              method="ML")
beta.cs.same <- coef(twstrs.cs.same)
sebeta.cs.same <- robust.cov(twstrs.cs.same,m)$se.model
V.cs.same <- getVarCov(twstrs.cs.same)  
Gamma.cs.same <- cov2cor(V.cs.same)

## > cbind(beta.cs.same,sebeta.cs.same)
##                    beta.cs.same sebeta.cs.same
## dose0                43.0228599      1.9877334
## dose5000             44.8138716      1.9871713
## dose10000            45.2910438      1.9586539
## dose0:week           -0.9552474      0.3415959
## dose5000:week        -2.2499120      0.3460461
## dose10000:week       -2.9780883      0.3362959
## dose0:weekplus        1.2513966      0.4220074
## dose5000:weekplus     3.0350951      0.4302921
## dose10000:weekplus    4.3073696      0.4152312

## > V.cs.same
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]   [,6]
## [1,] 148.48 109.93 109.93 109.93 109.93 109.93
## [2,] 109.93 148.48 109.93 109.93 109.93 109.93
## [3,] 109.93 109.93 148.48 109.93 109.93 109.93
## [4,] 109.93 109.93 109.93 148.48 109.93 109.93
## [5,] 109.93 109.93 109.93 109.93 148.48 109.93
## [6,] 109.93 109.93 109.93 109.93 109.93 148.48
##   Standard Deviations: 12.185 12.185 12.185 12.185 12.185 12.185 

## > Gamma.cs.same
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
## [1,] 1.00000 0.74035 0.74035 0.74035 0.74035 0.74035
## [2,] 0.74035 1.00000 0.74035 0.74035 0.74035 0.74035
## [3,] 0.74035 0.74035 1.00000 0.74035 0.74035 0.74035
## [4,] 0.74035 0.74035 0.74035 1.00000 0.74035 0.74035
## [5,] 0.74035 0.74035 0.74035 0.74035 1.00000 0.74035
## [6,] 0.74035 0.74035 0.74035 0.74035 0.74035 1.00000

#  Get AIC/BIC values

## > anova(twstrs.un,twstrs.un.same,twstrs.un.var,twstrs.un.var.same,
## twstrs.csh,twstrs.csh.same,twstrs.cs,twstrs.cs.same)
##                    Model df      AIC      BIC    logLik   Test   L.Ratio p-value
## twstrs.un              1 42 4109.012 4293.048 -2012.506                         
## twstrs.un.same         2 30 4109.029 4240.484 -2024.515 1 vs 2 24.017361  0.0202
## twstrs.un.var          3 27 4118.240 4236.549 -2032.120 2 vs 3 15.211133  0.0016
## twstrs.un.var.same     4 25 4128.334 4237.879 -2039.167 3 vs 4 14.093257  0.0009
## twstrs.csh             5 28 4140.936 4263.627 -2042.468 4 vs 5  6.602835  0.0857
## twstrs.csh.same        6 16 4140.680 4210.789 -2054.340 5 vs 6 23.743244  0.0220
## twstrs.cs              7 13 4162.248 4219.212 -2068.124 6 vs 7 27.568656  <.0001
## twstrs.cs.same         8 11 4160.962 4209.162 -2069.481 7 vs 8  2.714076  0.2574

#  Based on BIC, the model with common compound symmetric
#  correlation matrix with common variance that is the same over time
#  over time is preferred.  Thus, we can get the same analyses as in SAS

u <- twstrs.cs.same

#  We use model-based standard errors here, given the missing data

Sig.model <- robust.cov(u,m)$Sig.model
t.model <- beta.cs.same/sebeta.cs.same
df <- nrow(thedat)-length(beta.cs.same)    #  this is the residual df used
                                               #  by gls; SAS uses a different approximation
p.value <- round(2*pt(-abs(t.model),df),4)

#  Results with using residual degreees of freedom and t-test -
#  slightly different from those in SAS but qualitatively same

## > cbind(beta.cs.same,sebeta.cs.same,t.model,p.value)
##                    beta.cs.same sebeta.cs.same   t.model
## dose0                43.0228599      1.9877334 21.644181
## dose5000             44.8138716      1.9871713 22.551589
## dose10000            45.2910438      1.9586539 23.123556
## dose0:week           -0.9552474      0.3415959 -2.796426
## dose5000:week        -2.2499120      0.3460461 -6.501769
## dose10000:week       -2.9780883      0.3362959 -8.855560
## dose0:weekplus        1.2513966      0.4220074  2.965343
## dose5000:weekplus     3.0350951      0.4302921  7.053570
## dose10000:weekplus    4.3073696      0.4152312 10.373425
##                    p.value
## dose0               0.0000
## dose5000            0.0000
## dose10000           0.0000
## dose0:week          0.0053
## dose5000:week       0.0000
## dose10000:week      0.0000
## dose0:weekplus      0.0031
## dose5000:weekplus   0.0000
## dose10000:weekplus  0.0000

#  Function to contstruct F test statistics

f.test <- function(L,beta,Sig,dendf){
    numdf <- nrow(L)
    F <- t(L%*%beta)%*%solve(L%*%Sig%*%t(L))%*%(L%*%beta)/numdf
    p.value <- pf(F,numdf,dendf,lower.tail=FALSE)
    return(list(F=F,p.value=p.value,F.pvalue=c(F,round(p.value,4))))
}
    
#  Test of common intercept

L <- matrix(c(1,-1,0,0,0,0,0,0,0,
              1,0,-1,0,0,0,0,0,0),2,9,byrow=TRUE)

## > f.test(L,beta.cs.same,Sig.model,df)$F.pvalue
## [1] 0.3650811 0.6943000

#  Test of common 1st stage slope

L  <- matrix(c(0,0,0,1,-1,0,0,0,0,
               0,0,0,1,0,-1,0,0,0),2,9,byrow=TRUE)

## > f.test(L,beta.cs.same,Sig.model,df)$F.pvalue
## [1] 9.105642 0.000100

#  Test of common 2nd stage slope

L  <- matrix(c(0,0,0,1,-1,0,1,-1,0,
               0,0,0,1,0,-1,1,0,-1),2,9,byrow=TRUE)

## > f.test(L,beta.cs.same,Sig.model,df)$F.pvalue
## [1] 18.06999  0.00000

#  Test of difference in means at end of study (week 16)

L <- matrix(c(1,-1,0,16,-16,0,12,-12,0,
              1,0,-1,16,0,-16,12,0,-12),2,9,byrow=TRUE)

## >  f.test(L,beta.cs.same,Sig.model,df)$F.pvalue
## [1] 2.757605 0.064300

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

#  Estimates of mean @ week 16 for each dose and rate of change
#  at 8 months

#  mean @ week 16, Placebo

L <- matrix(c(1,0,0,16,0,0,12,0,0),1,9,byrow=TRUE)

## > est.L(L,beta.cs.same,Sig.model,0.05,df)
##      Estimate Std Error Lower CI Upper CI
## [1,] 42.75566  2.019271 38.78971 46.72161

#  mean @ week 16, 5000 U

L <- matrix(c(0,1,0,0,16,0,0,12,0),1,9,byrow=TRUE)

## > est.L(L,beta.cs.same,Sig.model,0.05,df)
##      Estimate Std Error Lower CI Upper CI
## [1,] 45.23642  2.029655 41.25008 49.22276

#  mean @ week 16, 10000 U

L <- matrix(c(0,0,1,0,0,16,0,0,12),1,9,byrow=TRUE)

## >  est.L(L,beta.cs.same,Sig.model,0.05,df)
##      Estimate Std Error Lower CI Upper CI
## [1,] 49.33007  1.984136 45.43313 53.22701

#  difference in means @ month 12

L  <- matrix(c(1,-1,0,12,-12,0,144,-144,0,
               1,0,-1,12,0,-12,144,0,-144),2,9,byrow=TRUE)

## > f.test(L,beta.un.var.same,Sig.model,df)$F.pvalue
## [1] 138.8577   0.0000

#  2nd stage slope, Placebo

L <- matrix(c(0,0,0,1,0,0,1,0,0),1,9,byrow=TRUE)

## > est.L(L,beta.cs.same,Sig.model,0.05,df)
##      Estimate Std Error   Lower CI  Upper CI
## [1,] 0.2961491  0.122987 0.05459664 0.5377017

#  2nd stage slope, 5000 U

L <- matrix(c(0,0,0,0,1,0,0,1,0),1,9,byrow=TRUE)

## > est.L(L,beta.cs.same,Sig.model,0.05,df)
##      Estimate Std Error  Lower CI Upper CI
## [1,] 0.7851832 0.1267187 0.5363014 1.034065

#  2nd stage slope, 10000 U

L <- matrix(c(0,0,0,0,0,1,0,0,1),1,9,byrow=TRUE)

## > est.L(L,beta.cs.same,Sig.model,0.05,df)
##      Estimate Std Error Lower CI Upper CI
## [1,] 1.329281  0.120257 1.093091 1.565472

#  Include age and gender effects

twstrs.cs.ag <- gls(twstrs ~ -1 + dose + dose:age + dose:gender +
              week:dose + week:dose:age + week:dose:gender+
              weekplus:dose,data=thedat,
              correlation=corCompSymm(form = ~ 1 | as.factor(id)),
              method="ML")
beta.cs.ag <- coef(twstrs.cs.ag)
sebeta.cs.ag <- robust.cov(twstrs.cs.ag,m)$se.model
V.cs.ag <- getVarCov(twstrs.cs.ag)  
Gamma.cs.ag <- cov2cor(V.cs.ag)

u <- twstrs.cs.ag

#  We use model-based standard errors here, given the missing data

Sig.model <- robust.cov(u,m)$Sig.model
t.model <- beta.cs.ag/sebeta.cs.ag
df <- nrow(thedat)-length(beta.cs.ag)    #  this is the residual df used
                                               #  by gls; SAS uses a different approximation
p.value <- round(2*pt(-abs(t.model),df),4)

## > cbind(beta.cs.ag,sebeta.cs.ag,t.model,p.value)
##                         beta.cs.ag sebeta.cs.ag    t.model
## dose0                 36.339338418  8.860697344  4.1011827
## dose5000              58.710995499  8.797073742  6.6739233
## dose10000             51.499755648  8.740198205  5.8922869
## dose0:age              0.153028864  0.152251263  1.0051074
## dose5000:age          -0.264541651  0.149303011 -1.7718440
## dose10000:age         -0.081429997  0.154273461 -0.5278289
## dose0:gender          -3.697405469  3.758963136 -0.9836238
## dose5000:gender        2.410453688  3.631687890  0.6637282
## dose10000:gender      -6.897265911  4.190851180 -1.6457912
## dose0:week            -1.890064689  0.530032995 -3.5659378
## dose5000:week         -2.886916411  0.561136739 -5.1447646
## dose10000:week        -3.499016776  0.514741851 -6.7976147
## dose0:weekplus         1.322832814  0.418351405  3.1620136
## dose5000:weekplus      3.091967538  0.428374684  7.2179044
## dose10000:weekplus     4.326857020  0.411166159 10.5233783
## dose0:age:week         0.016183955  0.007034380  2.3006940
## dose5000:age:week      0.010303640  0.007445843  1.3838111
## dose10000:age:week     0.008552257  0.007086039  1.2069163
## dose0:gender:week      0.062297225  0.169502053  0.3675308
## dose5000:gender:week   0.054284794  0.167013721  0.3250319
## dose10000:gender:week  0.175613383  0.198500440  0.8847002
##                       p.value
## dose0                  0.0000
## dose5000               0.0000
## dose10000              0.0000
## dose0:age              0.3153
## dose5000:age           0.0770
## dose10000:age          0.5978
## dose0:gender           0.3257
## dose5000:gender        0.5071
## dose10000:gender       0.1004
## dose0:week             0.0004
## dose5000:week          0.0000
## dose10000:week         0.0000
## dose0:weekplus         0.0017
## dose5000:weekplus      0.0000
## dose10000:weekplus     0.0000
## dose0:age:week         0.0218
## dose5000:age:week      0.1670
## dose10000:age:week     0.2280
## dose0:gender:week      0.7134
## dose5000:gender:week   0.7453
## dose10000:gender:week  0.3767
