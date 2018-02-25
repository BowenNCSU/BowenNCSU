########################################################################
#
#   Age-related macular degeneration data - ML analyses under MAR
#
#  There are 8 individuals with nonmonotone patterns; all the rest
#  exhibit monotone dropout.  The baseline measure is available on
#  everyone.
#
#  Several analyses are performed here on the visual acuity measures
#  (visual0 - visual52 in the data set).
#
#  (a) Fit a multivariate normal model to the available data, LOCF
#      data, and complete cases with distinct means at each of 0,
#      4, 12, 24, 52 weeks for each treatment and common unstructured
#      covariance matrix for both treatments by ML using lme
#
#  (b) Fit a multivariate normal model to the available data separately
#      by treatment group by ML using lme
#    
#  (c) Fit the same models to all available data using EM algorithm
#      in the norm package
#
########################################################################

library(norm)   #  to do EM
library(nlme)   #  to use gls()
library(zoo)    #  to get LOCF data
library(mice)   #  to use md.pattern function

armd.data <- read.table("armd.R.dat")

#  Extract just the columns of interest

armd.data <- armd.data[,c(1,7:11,13)]

#  Recode treatment

armd.data[,7] <- 0*(armd.data[,7]==1) + 1*(armd.data[,7]==4)
colnames(armd.data) <- c("id","y1","y2","y3","y4","y5","trt")
N <- nrow(armd.data)

#  Data set with complete cases only

armd.data.cc <- armd.data[complete.cases(armd.data),]
N.cc <- nrow(armd.data.cc)

#  LOCF data set - use the zoo package and function na.locf

visual <- as.matrix(armd.data[,2:6])
for (r in 1:nrow(visual)){
    v <- visual[r,]
    vz <- na.locf(zoo(v))
    visual[r,] <- vz
}
armd.data.locf <- cbind(armd.data[,1],visual,armd.data[,7])
armd.data.locf <- data.frame(armd.data.locf)
colnames(armd.data.locf) <- c("id","y1","y2","y3","y4","y5","trt")

# (a)  missingness patterns

## > md.pattern(armd.data[2:6])
##     y1 y2 y3 y4 y5   
## 188  1  1  1  1  1  0
##   2  1  0  1  1  1  1
##   4  1  1  1  0  1  1
##  24  1  1  1  1  0  1
##   1  1  1  0  0  1  2
##   8  1  1  1  0  0  2
##   1  1  0  1  0  0  3
##   6  1  1  0  0  0  3
##   6  1  0  0  0  0  4
##      0  9 13 26 45 93

#  these data sets are in the "wide" format; reconfigure in the "long" format

armd.data.alt <- reshape(armd.data,varying=list(names(armd.data)[2:6]),v.names="visual",times=c(0,4,12,24,52),timevar="weeks",idvar="id",direction="long")
armd.data.alt <- armd.data.alt[order(armd.data.alt$id),]

armd.data.cc.alt <- reshape(armd.data.cc,varying=list(names(armd.data.cc)[2:6]),v.names="visual",times=c(0,4,12,24,52),timevar="weeks",idvar="id",direction="long")
armd.data.cc.alt <- armd.data.cc.alt[order(armd.data.cc.alt$id),]

armd.data.locf.alt <- reshape(armd.data.locf,varying=list(names(armd.data.locf)[2:6]),v.names="visual",times=c(0,4,12,24,52),timevar="weeks",idvar="id",direction="long")
armd.data.locf.alt <- armd.data.locf.alt[order(armd.data.locf.alt$id),]

#  time variables for unstructured correlation in gls()

armd.data.alt$time <- as.numeric(factor(armd.data.alt$weeks))
armd.data.cc.alt$time <- as.numeric(factor(armd.data.cc.alt$weeks))
armd.data.locf.alt$time <- as.numeric(factor(armd.data.locf.alt$weeks))

#  (b)  Fit multivariate normal model to available data from both
#  groups using ML

pid <- as.factor(armd.data.alt$id)
week <- as.factor(armd.data.alt$weeks)
trt <- as.factor(armd.data.alt$trt)  

avail.mle <- gls(visual ~ -1 + week + week:trt,correlation=corSymm(form = ~ time | pid),
                 weights = varIdent(form = ~ 1 | week),method="ML",na.action=na.omit,data=armd.data.alt)
summary(avail.mle)

## Coefficients:
##               Value Std.Error  t-value p-value
## week0      55.33613  1.367406 40.46795  0.0000
## week4      54.05487  1.460954 36.99971  0.0000
## week12     52.98447  1.588584 33.35327  0.0000
## week24     49.31615  1.721048 28.65473  0.0000
## week52     44.02513  1.766815 24.91780  0.0000
## week0:trt  -0.75762  1.925797 -0.39341  0.6941
## week4:trt  -2.96186  2.063044 -1.43567  0.1514
## week12:trt -4.26556  2.253360 -1.89298  0.0586
## week24:trt -3.82741  2.453133 -1.56021  0.1190
## week52:trt -5.62378  2.545241 -2.20953  0.0273
## > getVarCov(avail.mle,individual=2)
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 220.50 201.61 188.79 182.45 143.79
## [2,] 201.61 250.77 228.67 219.40 175.34
## [3,] 188.79 228.67 295.71 262.48 224.76
## [4,] 182.45 219.40 262.48 342.44 290.86
## [5,] 143.79 175.34 224.76 290.86 350.16

#  (c) Fit same model to CC data using ML

pid <- as.factor(armd.data.cc.alt$id)
week <- as.factor(armd.data.cc.alt$weeks)
trt <- as.factor(armd.data.cc.alt$trt)  
cc.mle <- gls(visual ~  -1 + week + week:trt,correlation=corSymm(form = ~ time | pid),
                 weights = varIdent(form = ~ 1 | week),method="ML",data=armd.data.cc.alt)
summary(cc.mle)

## Coefficients:
##               Value Std.Error  t-value p-value
## week0      55.39216  1.473343 37.59624  0.0000
## week4      54.47059  1.551355 35.11163  0.0000
## week12     53.07843  1.673451 31.71795  0.0000
## week24     49.79412  1.808359 27.53553  0.0000
## week52     44.43137  1.835730 24.20365  0.0000
## week0:trt  -0.54332  2.178380 -0.24941  0.8031
## week4:trt  -2.86594  2.293722 -1.24947  0.2118
## week12:trt -2.89238  2.474245 -1.16900  0.2427
## week24:trt -3.27086  2.673710 -1.22334  0.2215
## week52:trt -4.71044  2.714179 -1.73549  0.0830

## > getVarCov(cc.mle,individual=2)
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 219.06 205.22 194.81 191.84 151.26
## [2,] 205.22 242.87 221.54 217.76 174.92
## [3,] 194.81 221.54 282.61 254.05 218.11
## [4,] 191.84 217.76 254.05 330.01 279.89
## [5,] 151.26 174.92 218.11 279.89 340.07

#  (d) Fit same model to LOCF data using ML

pid <- as.factor(armd.data.locf.alt$id)
week <- as.factor(armd.data.locf.alt$weeks)
trt <- as.factor(armd.data.locf.alt$trt)  

locf.mle <- gls(visual ~ -1 + week + week:trt,correlation=corSymm(form = ~ time | pid),
                weights = varIdent(form = ~ 1 |week),method="ML",data=armd.data.locf.alt)
summary(locf.mle)

## Coefficients:
##               Value Std.Error  t-value p-value
## week0      55.33613  1.366928 40.48211  0.0000
## week4      54.05882  1.457555 37.08869  0.0000
## week12     53.01681  1.579395 33.56780  0.0000
## week24     49.46218  1.700184 29.09226  0.0000
## week52     44.70588  1.729726 25.84564  0.0000
## week0:trt  -0.75762  1.925124 -0.39354  0.6940
## week4:trt  -2.78610  2.052759 -1.35724  0.1750
## week12:trt -3.91763  2.224353 -1.76125  0.0785
## week24:trt -2.87541  2.394467 -1.20086  0.2300
## week52:trt -3.69762  2.436073 -1.51786  0.1293

## > getVarCov(locf.mle)
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 220.50 202.55 190.01 184.35 152.43
## [2,] 202.55 250.70 230.54 224.94 186.46
## [3,] 190.01 230.54 294.37 267.19 232.11
## [4,] 184.35 224.94 267.19 341.12 296.17
## [5,] 152.43 186.46 232.11 296.17 353.08

## (e)  Comparison of estimates of beta.5

## Available  week52:trt -5.62378  2.545241 -2.20953  0.0273
## Comp Case  week52:trt -4.71044  2.714179 -1.73549  0.0830
## LOCF       week52:trt -3.69762  2.436073 -1.51786  0.1293

#  (f)  Now fit multivariate normal model to avaiable data by treatment

armd.data.alt0 <- armd.data.alt[armd.data.alt$trt==0,]
armd.data.alt1 <-armd.data.alt[armd.data.alt$trt==1,]

#   Placebo (trt == 0)

pid <- as.factor(armd.data.alt0$id)
week <- as.factor(armd.data.alt0$weeks)

trt0.mle <- gls(visual ~  -1 + week,correlation=corSymm(form = ~ time | pid),
                 weights = varIdent(form = ~ 1 |week),method="ML",na.action=na.omit,data=armd.data.alt0)
summary(trt0.mle)

## Coefficients:
##           Value Std.Error  t-value p-value
## week0  55.33613  1.375417 40.23226       0
## week4  54.05293  1.452678 37.20918       0
## week12 52.98078  1.576440 33.60786       0
## week24 49.34091  1.748577 28.21775       0
## week52 44.03456  1.764363 24.95777       0

## > getVarCov(trt0.mle,individual=4)
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 223.15 206.17 188.80 188.81 148.25
## [2,] 206.17 248.10 223.43 219.36 176.35
## [3,] 188.80 223.43 291.25 251.30 217.76
## [4,] 188.81 219.36 251.30 353.02 296.47
## [5,] 148.25 176.35 217.76 296.47 349.42

#  Active Treatment (trt == 1)

pid <- as.factor(armd.data.alt1$id)
week <- as.factor(armd.data.alt1$weeks)

trt1.mle <- gls(visual ~  -1 + week,correlation=corSymm(form = ~ time | pid),
                 weights = varIdent(form = ~ 1 |week),method="ML",na.action=na.omit,data=armd.data.alt1)
summary(trt1.mle)

## Coefficients:
##           Value Std.Error  t-value p-value
## week0  54.57851  1.348227 40.48171       0
## week4  51.09662  1.465549 34.86517       0
## week12 48.71815  1.609522 30.26870       0
## week24 45.52128  1.702002 26.74573       0
## week52 38.45591  1.836415 20.94076       0

## >  getVarCov(trt1.mle,individual=2)
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]   [,5]
## [1,] 217.90 196.93 189.01 174.60 139.30
## [2,] 196.93 253.38 233.97 217.40 173.76
## [3,] 189.01 233.97 300.14 272.19 230.92
## [4,] 174.60 217.40 272.19 326.02 282.19
## [5,] 139.30 173.76 230.92 282.19 351.13

#  (g) Fit the same models to each treatment using the EM algorithm

#   Placebo (trt == 0)

#  Get the visual acuity measures for trt == 0

visual <- as.matrix(armd.data[armd.data[,7]==0,2:6])

#  Run prelim.norm to set up the data for the algorithm

prelim.em0 <- prelim.norm(visual)

#  Run the EM algorithm initially with default settings to get the form
#  of the starting values

theta.init0 <- em.norm(prelim.em0,showits=FALSE)

 #  Extract the parameters to see the format; if the user wants to change
#  these values to something else /she can; we won't bother

theta.init0 <- getparam.norm(prelim.em0,theta.init0,corr=TRUE)

#  Convert these values to the form used by the function em.norm

theta.init0 <- makeparam.norm(prelim.em0,theta.init0)

trt0.em <- em.norm(prelim.em0,start=theta.init0,showits=TRUE,maxits=200,criterion=1e-5)
theta.em.trt0 <- getparam.norm(prelim.em0,trt0.em,corr=TRUE)

#  Final parameter values for trt == 0 by EM

## > theta.em.trt0$mu       
## [1] 55.33613 54.05293 52.98078 49.34093 44.03457

##  covariance matrix

theta.em.trt0.cov <- diag(theta.em.trt0$sdv)%*%theta.em.trt0$r%*%diag(theta.em.trt0$sdv)
> theta.em.trt0.cov
##          [,1]     [,2]     [,3]     [,4]     [,5]
## [1,] 223.1475 206.1692 188.7986 188.8160 148.2489
## [2,] 206.1692 248.1052 223.4291 219.3650 176.3470
## [3,] 188.7986 223.4291 291.2531 251.3030 217.7606
## [4,] 188.8160 219.3650 251.3030 353.0215 296.4667
## [5,] 148.2489 176.3470 217.7606 296.4667 349.4235

#   Activex Treatment (trt == 1) -- use the same steps

#  Get the visual acuity measures for trt == 0

visual <- as.matrix(armd.data[armd.data[,7]==1,2:6])

prelim.em1 <- prelim.norm(visual)
theta.init1 <- em.norm(prelim.em1,showits=FALSE)
theta.init1 <- getparam.norm(prelim.em1,theta.init1,corr=TRUE)
theta.init1 <- makeparam.norm(prelim.em1,theta.init1)
trt1.em <- em.norm(prelim.em1,start=theta.init1,showits=TRUE,maxits=200,criterion=1e-5)
theta.em.trt1 <- getparam.norm(prelim.em1,trt1.em,corr=TRUE)

#  Final parameter values for trt == 1 by EM

## > theta.em.trt1$mu     
## [1] 54.57851 51.09662 48.71815 45.52130 38.45595

##  covariance matrix

## > theta.em.trt1.cov
##          [,1]     [,2]     [,3]     [,4]     [,5]
## [1,] 217.8967 196.9365 189.0164 174.5998 139.3036
## [2,] 196.9365 253.3823 233.9794 217.4068 173.7648
## [3,] 189.0164 233.9794 300.1460 272.1957 230.9290
## [4,] 174.5998 217.4068 272.1957 326.0305 282.1967
## [5,] 139.3036 173.7648 230.9290 282.1967 351.1272



