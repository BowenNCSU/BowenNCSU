##  Homework 4, Problem 2.  Fit the nonlinear model to the HCV data
##  from "Subject 3" using gnls() with efficacy parameterized by
##  logistic function

library(nlme)

##  define the HCV dynamic model function 

meanfunc <- function(x,b1,b2,b3){

    v0 <- exp(b1)
    cc <- exp(b2)
    e <- exp(b3)/(1+exp(b3))
    t0 <- 0.20      #  t0 is fixed
    tt <- (x>t0)*(x-t0)
    fpl <- v0*(1-e+e*exp(-cc*tt))

    return(fpl)

##  compute analytical dervivatives

   meangrad <- array(0,c(length(x),3),list(NULL,c("b1","b2","b3")))
   meangrad[,"b1"] <- fpl
   meangrad[,"b2"] <- -v0*e*exp(-cc*tt)*cc*tt
   meangrad[,"b3"] <- v0*e*(1-e)*(exp(-cc*tt)-1)
   attr(fpl,"gradient") <- meangrad

}

##  read in the data

thedat <- read.table("hcv3.dat")
colnames(thedat) <-  c("days","vl106")

##  rescaled viral load -- alternatively we could keep it on the
##  original scale and redefine the parameter v0 in the model to
##  scaled viral load

vl <- thedat$vl106/10^6

##  Call the gnls() function to implement the fit.  A starting value
##  for beta must be provided for a nonlinear model like this.
##  Reasonable starting values can be deduced by considerint the form
##  of the model

hcv3.gls <- gnls(vl ~ meanfunc(days,b1,b2,b3),
   start = c(b1=log(4),b2=1.5,b3=1.9),
   weights = varPower(1,form = ~ fitted(.)),
   data=thedat)

## >  summary (hcv3.gls)
## Generalized nonlinear least squares fit
##   Model: vl ~ meanfunc(days, b1, b2, b3) 
##   Data: thedat 
##         AIC       BIC   logLik
##   -5.790593 -3.801117 7.895297

## Variance function:
##  Structure: Power of variance covariate
##  Formula: ~fitted(.) 
##  Parameter estimates:
##     power 
## 0.8241415 

## Coefficients:
##       Value  Std.Error  t-value p-value
## b1 1.438664 0.03891406 36.97028       0
## b2 1.312434 0.08478266 15.47998       0
## b3 1.816623 0.08980596 20.22831       0

##  Correlation: 
##    b1     b2    
## b2  0.479       
## b3  0.342 -0.384

## Standardized residuals:
##        Min         Q1        Med         Q3        Max 
## -1.6069731 -0.2661226  0.3801448  0.4432829  1.2127223 

## Residual standard error: 0.09515951 
## Degrees of freedom: 11 total; 8 residual

# Extract the parameters

beta <- coef(hcv3.gls)
delta <- hcv3.gls$modelStruct$varStruct[1]
sigma <- hcv3.gls$sigma

## > beta
##       b1       b2       b3 
## 1.438663 1.312434 1.816623

## > c(delta,sigma)
## [1] 0.82414151 0.09515951

## Get delta method standard errors

SEs <- sqrt(diag(hcv3.gls$varBeta))

v0 <- exp(beta[1])*10^6
SEv0 <- exp(beta[1])*SEs[1]*10^6

cc <- exp(beta[2])
SEcc <- exp(beta[2])*SEs[2]

eff <- exp(beta[3])/(1+exp(beta[3]))
SEeff <- eff*(1-eff)*SEs[3]

parms <- c(v0,cc,eff)
SEparms <- c(SEv0,SEcc,SEeff)

## > round(cbind(parms,SEparms),5)
##           parms     SEparms
## b1 4.215059e+06 1.64025e+05
## b2 3.715200e+00 3.14980e-01
## b3 8.601600e-01 1.08000e-02

## Plot of the data with fit superimposed

tgrid <- seq(0,2,0.01) 
fgrid <- 10^6*meanfunc(tgrid,beta[1],beta[2],beta[3])

pdf("hcv3_plot.pdf",width=8)
plot(thedat$days,thedat$vl106,,xlab="Days",ylab="Viral Load (copies/ml)",
     cex=1.2,pch=18,cex.lab=1.2)
lines(tgrid,fgrid,lty=2)
dev.off()
