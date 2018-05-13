##  Homework 5 PK analysis

## Read in data

thedat1 <- read.table("pk.dat")
colnames(thedat1) <- c("id","time","conc","weight","age","metab")

##  plot the PK data on the same graph

pdf("pkdata_plot.pdf",width=8)
yname <- "Concentration (mg/L)"
xname <- "Time (hours)"
i <- 1
ti <- thedat1[thedat1[,1]==i,2]
yi <- thedat1[thedat1[,1]==i,3]
plot(ti,yi,type="b",xlab=xname,ylab=yname,ylim=c(0,max(thedat1[,3])),lty=2,pch=18)
for (i in 2:m1){
    ti <- thedat1[thedat1[,1]==i,2]
    yi <- thedat1[thedat1[,1]==i,3]
    points(ti,yi,type="b",lty=2,pch=18)
}
dev.off()

## covariates for each individual

wt <- thedat1[thedat1[,2]==0.25,4]
metab <- thedat1[thedat1[,2]==0.25,6]
ag <- thedat1[thedat1[,2]==0.25,5]

library(nlme)

pkmodel <- function(x,b1,b2,b3,dose){
  ka <- exp(b1)
  cl <- exp(b2)
  v <- exp(b3)
  ke <- cl/v
  kk <- ka-ke
  f1 <- dose*ka*(exp(-ke*x)-exp(-ka*x))/(v*kk)
  tt <- dose*x*ka/(v*kk)
  return(f1)
  
#  compute analytical dervivatives -- create the gradient matrix X
        
   meangrad <- array(0,c(length(x),3),list(NULL,c("b1","b2","b3")))
   meangrad[,"b1"] <- -f1*(ke/kk)+tt*ka*exp(-ka*x)
   meangrad[,"b2"] <- f1*(ke/kk)-tt*ke*exp(-ke*x)
   meangrad[,"b3"] <- -f1*(ka/kk)+tt*ke*exp(-ke*x)
   attr(f1,"gradient") <- meangrad
  
}

logwt <- log(thedat1$weight)  #  could also fit with weight on log scale
meta <- as.factor(thedat1$metab)

## Model with no covariates

pkfit.nocov <- nlme(conc ~ pkmodel(time,b1,b2,b3,1000),
    fixed = list(b1 ~ 1,b2 ~ 1,b3 ~ 1),
    random = list(b1 ~ 1,b2 ~ 1,b3 ~ 1),
    groups = ~ id,
    weights = varPower(1),
    start = list(fixed = c(0.3,-0.2,1.6)),
    data = thedat1,                  
    method="ML",verbose=TRUE)

## > summary(pkfit.nocov)
## Nonlinear mixed-effects model fit by maximum likelihood
##   Model: conc ~ pkmodel(time, b1, b2, b3, 1000) 
##  Data: thedat1 
##        AIC      BIC    logLik
##   4440.089 4489.096 -2209.044

## Random effects:
##  Formula: list(b1 ~ 1, b2 ~ 1, b3 ~ 1)
##  Level: id
##  Structure: General positive-definite, Log-Cholesky parametrization
##          StdDev     Corr       
## b1       0.11093363 b1    b2   
## b2       0.55249548 0.339      
## b3       0.14395331 0.575 0.337
## Residual 0.08558814            

## Variance function:
##  Structure: Power of variance covariate
##  Formula: ~fitted(.) 
##  Parameter estimates:
##    power 
## 1.064212 
## Fixed effects: list(b1 ~ 1, b2 ~ 1, b3 ~ 1) 
##         Value  Std.Error  DF  t-value p-value
## b1  0.2180528 0.02250306 581  9.68992  0.0000
## b2 -0.1218259 0.07630353 581 -1.59660  0.1109
## b3  1.5955573 0.02173870 581 73.39710  0.0000
##  Correlation: 
##    b1    b2   
## b2 0.258      
## b3 0.585 0.330

D.T <- diag(VarCorr(pkfit.nocov)[1:3,2])
D.corrs <- as.numeric(VarCorr(pkfit.nocov)[2:3,3:4])
D.Gam <- diag(3)
D.Gam[1,2]<-D.corrs[1]; D.Gam[1,3]<-D.corrs[2]; D.Gam[2,3]<-D.corrs[4]
 D.Gam[2,1]<-D.Gam[1,2]; D.Gam[3,1]<-D.Gam[1,3]; D.Gam[3,2]<-D.Gam[2,3]
D <- D.T%*%D.Gam%*%D.T

## > D
##             [,1]       [,2]        [,3]
## [1,] 0.012306270 0.02077742 0.009182326
## [2,] 0.020777422 0.30525126 0.026802807
## [3,] 0.009182326 0.02680281 0.020722555
## > D.Gam
##       [,1]  [,2]  [,3]
## [1,] 1.000 0.339 0.575
## [2,] 0.339 1.000 0.337
## [3,] 0.575 0.337 1.000
## > sqrt(diag(D))
## [1] 0.1109336 0.5524955 0.1439533
## > pkfit.nocov$sigma^2
## [1] 0.007325329

## get the fixed and random effects estimates

pk.beta <- pkfit.nocov$coefficients$fixed
pk.randef <- random.effects(pkfit.nocov)

b1.nocov <- pk.randef$b1
b2.nocov <- pk.randef$b2
b3.nocov <- pk.randef$b3

## Plot them against covariates to look for evidence of
## associations

plotre <- b3.nocov    ##  change this to plot different components
                      ##  against covariates
boxplot(plotre ~ metab)
plot(wt,plotre)
plot(ag,plotre)

## Model with all covariates associated with all 3 PK parameters

## The model fit is VERY sensitive to starting values with respect to
## D.  I gave it the same starting values are for the nlinmix macro
## got an almost singular D matrix.  Here, I am giving it starting
## values based on the final estimates from nlinmix and it is better
## behaved but still yields a singular D matrix

pkfit.allcov <- nlme(conc ~ pkmodel(time,b1,b2,b3,1000),
    fixed = list(b1 ~ 1 + meta + age + weight,b2 ~ 1+ meta + age +weight,b3 ~ 1+meta+age+weight),
    random = list(b1 ~ 1,b2 ~ 1,b3 ~ 1),
    groups = ~ id,
    weights = varPower(1),
##  start = list(fixed = c(-0.3,0,0,0,0,-0.2,0,0,0,0,1.6,0,0,0,0)),    
    start = list(fixed =  c(0.2,-0.03,0.07,0.001,-0.0007,-0.8,0.9,1.5,-0.009,0.001,1.8,0.0005,
    0.08,0.0009,-0.004)),
    data = thedat1, method="ML",verbose=TRUE)

## > summary(pkfit.allcov)
## Nonlinear mixed-effects model fit by maximum likelihood
##   Model: conc ~ pkmodel(time, b1, b2, b3, 1000) 
##  Data: thedat1 
##        AIC      BIC    logLik
##   4235.977 4338.447 -2094.989

## Random effects:
##  Formula: list(b1 ~ 1, b2 ~ 1, b3 ~ 1)
##  Level: id
##  Structure: General positive-definite, Log-Cholesky parametrization
##                StdDev     Corr         
## b1.(Intercept) 0.02468919 b1.(I) b2.(I)
## b2.(Intercept) 0.10460455 0.710        
## b3.(Intercept) 0.10202909 1.000  0.701 
## Residual       0.08898988              

## Variance function:
##  Structure: Power of variance covariate
##  Formula: ~fitted(.) 
##  Parameter estimates:
##    power 
## 1.061075 
## Fixed effects: list(b1 ~ 1 + meta + age + weight, b2 ~ 1 + meta + age + weight,      b3 ~ 1 + meta + age + weight) 
##                     Value  Std.Error  DF   t-value p-value
## b1.(Intercept)  0.1798051 0.09039888 569   1.98902  0.0472
## b1.meta2       -0.0338067 0.05946406 569  -0.56852  0.5699
## b1.meta3        0.0706506 0.05240790 569   1.34809  0.1782
## b1.age          0.0010714 0.00106511 569   1.00587  0.3149
## b1.weight      -0.0005198 0.00085492 569  -0.60803  0.5434
## b2.(Intercept) -0.8050304 0.07988085 569 -10.07789  0.0000
## b2.meta2        0.9107101 0.04910041 569  18.54791  0.0000
## b2.meta3        1.4634147 0.04365928 569  33.51898  0.0000
## b2.age         -0.0093625 0.00098242 569  -9.53002  0.0000
## b2.weight       0.0011824 0.00079373 569   1.48967  0.1369
## b3.(Intercept)  1.8336977 0.08572887 569  21.38950  0.0000
## b3.meta2        0.0000446 0.05425866 569   0.00082  0.9993
## b3.meta3        0.0796903 0.04772087 569   1.66993  0.0955
## b3.age          0.0008564 0.00104396 569   0.82034  0.4124
## b3.weight      -0.0040926 0.00084287 569  -4.85562  0.0000

D.T <- diag(VarCorr(pkfit.allcov)[1:3,2])
D.corrs <- as.numeric(VarCorr(pkfit.allcov)[2:3,3:4])
D.Gam <- diag(3)
D.Gam[1,2]<-D.corrs[1]; D.Gam[1,3]<-D.corrs[2]; D.Gam[2,3]<-D.corrs[4]
 D.Gam[2,1]<-D.Gam[1,2]; D.Gam[3,1]<-D.Gam[1,3]; D.Gam[3,2]<-D.Gam[2,3]
D <- D.T%*%D.Gam%*%D.T

## It looks like putting all these covariates in is causing problems;
## the correlation between log ka and log Cl is 1. This did not happen
## with the nlinmix macro in SAS nor with proc nlmixed and Laplace.

## > D
##              [,1]        [,2]        [,3]
## [1,] 0.0006095561 0.001833647 0.002519016
## [2,] 0.0018336471 0.010942112 0.007481568
## [3,] 0.0025190156 0.007481568 0.010409935
## > D.Gam
##      [,1]  [,2]  [,3]
## [1,] 1.00 0.710 1.000
## [2,] 0.71 1.000 0.701
## [3,] 1.00 0.701 1.000

## Reduced model with only the associations suggested by fit with all
## covariates and the plots of empirical Bayes estimates
## vs. covariates -- base starting values on fits from nlinmix and
## proc nlmixed in SAS

pkfit.cov1 <- nlme(conc ~ pkmodel(time,b1,b2,b3,1000),
    fixed = list(b1 ~ 1,b2 ~ 1 + meta + age,b3 ~ 1 + weight),
    random = list(b1 ~ 1,b2 ~ 1,b3 ~ 1),
    groups = ~ id,
    weights = varPower(1),
    start = list(fixed = c(0.2,-0.7,0.9,1.4,-0.01,1.97,-0.004)),
    data = thedat1, method="ML",verbose=TRUE)

## > summary(pkfit.cov1)
## Nonlinear mixed-effects model fit by maximum likelihood
##   Model: conc ~ pkmodel(time, b1, b2, b3, 1000) 
##  Data: thedat1 
##        AIC      BIC    logLik
##   4233.173 4300.001 -2101.586

## Random effects:
##  Formula: list(b1 ~ 1, b2 ~ 1, b3 ~ 1)
##  Level: id
##  Structure: General positive-definite, Log-Cholesky parametrization
##                StdDev     Corr        
## b1             0.03403102 b1    b2.(I)
## b2.(Intercept) 0.10982711 0.624       
## b3.(Intercept) 0.10957916 0.993 0.714 
## Residual       0.08965207             

## Variance function:
##  Structure: Power of variance covariate
##  Formula: ~fitted(.) 
##  Parameter estimates:
##    power 
## 1.060835 
## Fixed effects: list(b1 ~ 1, b2 ~ 1 + meta + age, b3 ~ 1 + weight) 
##                     Value  Std.Error  DF   t-value p-value
## b1              0.2284121 0.01764076 577  12.94797       0
## b2.(Intercept) -0.6701238 0.04436608 577 -15.10442       0
## b2.meta2        0.9170534 0.03831430 577  23.93502       0
## b2.meta3        1.4187276 0.03504851 577  40.47898       0
## b2.age         -0.0097412 0.00067894 577 -14.34755       0
## b3.(Intercept)  1.9656331 0.04451984 577  44.15184       0
## b3.weight      -0.0046201 0.00051571 577  -8.95873       0

D.T <- diag(VarCorr(pkfit.cov1)[1:3,2])
D.corrs <- as.numeric(VarCorr(pkfit.cov1)[2:3,3:4])
D.Gam <- diag(3)
D.Gam[1,2]<-D.corrs[1]; D.Gam[1,3]<-D.corrs[2]; D.Gam[2,3]<-D.corrs[4]
D.Gam[2,1]<-D.Gam[1,2]; D.Gam[3,1]<-D.Gam[1,3]; D.Gam[3,2]<-D.Gam[2,3]
D <- D.T%*%D.Gam%*%D.T

## Note the almost perfect correlation between b1 and b3; something
## weird is going on here.  The fits in SAS did not result in this, so
## I am skeptical this is reflecting a real phenomenon.

## > D
##             [,1]        [,2]        [,3]
## [1,] 0.001158110 0.002332218 0.003702987
## [2,] 0.002332218 0.012061994 0.008592820
## [3,] 0.003702987 0.008592820 0.012007592
## > D.Gam
##       [,1]  [,2]  [,3]
## [1,] 1.000 0.624 0.993
## [2,] 0.624 1.000 0.714
## [3,] 0.993 0.714 1.000

## Despite the weirdness, compare the model fits.  Both AIC and BIC
## choose the model with metab and age in clearance and log weight in
## volume

## > anova(pkfit.nocov,pkfit.allcov,pkfit.cov1)
##              Model df      AIC      BIC    logLik   Test
## pkfit.nocov      1 11 4440.089 4489.096 -2209.044       
## pkfit.allcov     2 23 4235.977 4338.447 -2094.989 1 vs 2
## pkfit.cov1       3 15 4233.173 4300.001 -2101.586 2 vs 3
##               L.Ratio p-value
## pkfit.nocov                  
## pkfit.allcov 228.1116  <.0001
## pkfit.cov1    13.1956  0.1053

## Use the preferred fit to estimate the PK parameters for everyone

pk.beta <- pkfit.cov1$coefficients$fixed
pk.randef <- random.effects(pkfit.cov1)

b1.cov <- pk.randef$b1
b2.cov <- pk.randef$b2
b3.cov <- pk.randef$b3

#  Here is an example of getting the estimate of mean clearance on the
#  original scale and a delta method standard error for it based on
#  the linear approximation in the solution.  We evaluate at the mean
#  value of age for each level of CYP2D6

avgage <- mean(ag)

cov.beta <- pkfit.cov1$varFix  #  covariance matrix of betahat

#  Estimate for poor metabolizers

a.Cl <- matrix(c(0,1,0,0,avgage,0,0),nrow=7)
Cl <- exp(t(a.Cl)%*%pk.beta)
se.Cl <- Cl*sqrt(t(a.Cl)%*%cov.beta%*%a.Cl)

## > cbind(Cl,se.Cl)
##           [,1]       [,2]
## [1,] 0.3159924 0.01032256

#  Estimate for extensive metabolizers

a.Cl <- matrix(c(0,1,0,1,avgage,0,0),nrow=7)
Cl <- exp(t(a.Cl)%*%pk.beta)
se.Cl <- Cl*sqrt(t(a.Cl)%*%cov.beta%*%a.Cl)

## > cbind(Cl,se.Cl)
##          [,1]       [,2]
## [1,] 1.305636 0.02480321

## Estimates of ka, Cl, and V for each individual

kai <- exp(pk.beta[1]+b1.cov)
Cli <- exp(pk.beta[2]+pk.beta[3]*(metab==2)+pk.beta[4]*(metab==3)+pk.beta[5]*ag
        + b2.cov)
Vi <- exp(pk.beta[6]+pk.beta[7]*wt+b3.cov)

## > cbind(kai,Cli,Vi)
##            kai       Cli       Vi
##  [1,] 1.313021 1.4977937 5.603882
##  [2,] 1.301786 0.7812308 5.307079
##  [3,] 1.383658 1.7414557 6.333534
##  [4,] 1.270810 1.4527079 5.826962
##  [5,] 1.318259 1.1425257 6.398643
##  [6,] 1.241582 1.3055414 5.474610
##  [7,] 1.208106 0.8682660 4.632686
##  [8,] 1.242064 1.3043519 5.174902
##  [9,] 1.276723 1.5611515 4.394933
## [10,] 1.201741 0.9450485 4.411419
## [11,] 1.232237 0.3034487 3.908335
## [12,] 1.270218 0.7668584 4.877640
## [13,] 1.256464 1.2003958 5.327050
## [14,] 1.201956 0.2790574 4.606822
## [15,] 1.232878 0.3309084 4.135249
## [16,] 1.202903 0.6345198 3.706770
## [17,] 1.283489 1.1532461 5.971166
## [18,] 1.224269 1.1511824 4.150786
## [19,] 1.265220 0.3083931 4.412071
## [20,] 1.292847 0.6322565 4.486967
## [21,] 1.209744 0.2657581 4.927266
## [22,] 1.201219 0.2573191 4.508633
## [23,] 1.228871 1.3425713 4.783875
## [24,] 1.231438 1.4423512 5.245055
## [25,] 1.226671 0.9460508 5.193241
## [26,] 1.212151 1.1592987 4.573093
## [27,] 1.258531 1.2133557 5.249810
## [28,] 1.271964 1.1894630 4.824491
## [29,] 1.255325 0.8649583 5.086885
## [30,] 1.309922 1.6362522 5.150922
## [31,] 1.259235 0.6440832 5.061668
## [32,] 1.245812 1.0542232 4.990044
## [33,] 1.262552 1.1365936 4.911798
## [34,] 1.255476 0.7827338 4.755812
## [35,] 1.247378 1.3018653 4.015200
## [36,] 1.339664 0.4267755 6.933870
## [37,] 1.224313 0.3718100 4.575348
## [38,] 1.236286 1.5696203 4.785379
## [39,] 1.182029 1.1526863 4.493734
## [40,] 1.310685 1.2232475 5.774319
## [41,] 1.250794 0.8394140 4.534189
## [42,] 1.220172 1.3587757 4.181253
## [43,] 1.251192 0.7822137 5.499755
## [44,] 1.302557 1.1165398 5.162102
## [45,] 1.296707 0.3985255 5.770890
## [46,] 1.309633 0.9489612 5.648099
## [47,] 1.297724 1.4411309 5.619642
## [48,] 1.255770 0.8678926 4.599822
## [49,] 1.270136 1.5798142 5.363696
## [50,] 1.261064 1.7245072 4.738264
## [51,] 1.242732 0.3281806 4.884020
## [52,] 1.246216 1.5297260 4.545024
## [53,] 1.238009 0.8455456 5.310134

## Compare these to those from nlinmix 

## First 2 individuals

beta.sas <- c(0.2234,-0.6663,1.9688,0.9145,1.4157,-0.00978,-0.00469)

b.id1 <- c(0.1206,0.08102,0.1637)
b.id2 <- c(-0.06255,-0.1117,0.05250)

## compare these to those from nlme fit; very different because D
## matrices are very different

##>  cbind(b1.cov[1:2],b2.cov[1:2],b3.cov[1:2])
##            [,1]        [,2]       [,3]
## [1,] 0.04391888  0.06451941 0.13667875
## [2,] 0.03532502 -0.09442540 0.08734298

ka1 <- exp(beta.sas[1]+b.id1[1])
Cl1 <- exp(beta.sas[2]+beta.sas[4]*(metab[1]==2)+beta.sas[5]*(metab[1]==3)+beta.sas[6]*ag[1]+b.id1[2])
V1 <- exp(beta.sas[3]+beta.sas[7]*wt[1]+b.id1[3])

ka2 <- exp(beta.sas[1]+b.id2[1])
Cl2 <- exp(beta.sas[2]+beta.sas[4]*(metab[2]==2)+beta.sas[5]*(metab[2]==3)+beta.sas[6]*ag[2]+b.id2[2])
V2 <- exp(beta.sas[3]+beta.sas[7]*wt[2]+b.id2[3])

## > rbind(c(1,ka1,Cl1,V1),c(2,ka2,Cl2,V2))
##      [,1]     [,2]     [,3]     [,4]
## [1,]    1 1.410579 1.521444 5.742646
## [2,]    2 1.174509 0.767605 5.111848


