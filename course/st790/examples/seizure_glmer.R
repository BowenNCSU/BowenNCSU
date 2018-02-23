#########################################################################
#
#   CHAPTER 9, EXAMPLE 5, Epileptic Seizure Study 
# 
#   Fitting subject-specific generalized linear mixed effects model
#   using the glmer() function in the lme4 package; you can install.packages("lme4")
#
#   The method is to evaluate the likelihood using adaptive Gaussian quadrature
#   to evaluate the integrals.  The option nAGQ specifies the number
#   of abscissae in each dimension of the random effects.  The default
#   is nAGQ = 1, which corresponds the Laplace approximation (roughly PQL).
#   Apparently, adaptive Gaussian quadrature is implemented only for 
#   a single scalar random effect in the currently available version
#   of glmer(), and the number of quadrature points nAGQ should be <= 25.
#
#   The rest of the syntax is the same as for the glm() function for fitting
#   ordinary generalized linear models and for lmer() for specifying fixed and 
#   random effects.  See the documentation on CRAN for more.
#
#   na.action specifies what should happen for missing
#   observations. We really don't need it here, as the data are balanced
#
#   Note that adaptive Gaussian quadrature and the Laplace approximation
#   yield virtually IDENTICAL results in this case.  See if you can
#   figure out why.
#
#########################################################################

library(lme4)

#  Read in the data

thedat <- read.table("seize.dat")
colnames(thedat) <- c("subj","seize","visit","trt","base","age")

#  Create other covariates

logage <- log(thedat$age)
o <- 8*(thedat$visit==0)+2*(thedat$visit>0)
logo <- log(o)
v <- as.numeric(thedat$visit>0)
v4 <- as.numeric(thedat$visit==4)

#  Fit the model with 2 random effects for intercept and time.  The
#  only option is to fit using the Laplace approximation

seize.glmm.laplace <- glmer(seize ~ offset(logo) + v + trt + v:trt + (1 + v | subj),
     family=poisson(link="log"), nAGQ=1, data=thedat, na.action=na.omit)

## > summary(seize.glmm.laplace)
## Generalized linear mixed model fit by maximum likelihood (Laplace
##   Approximation) [glmerMod]
##  Family: poisson  ( log )
## Formula: seize ~ offset(logo) + v + trt + v:trt + (1 + v | subj)
##    Data: thedat

##      AIC      BIC   logLik deviance df.resid 
##   1863.3   1889.1   -924.7   1849.3      288 

## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.1388 -0.7118 -0.0607  0.5189  6.9652 

## Random effects:
##  Groups Name        Variance Std.Dev. Corr
##  subj   (Intercept) 0.4999   0.7070       
##         v           0.2319   0.4815   0.17
## Number of obs: 295, groups:  subj, 59

## Fixed effects:
##              Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  1.071299   0.140267   7.638 2.21e-14 ***
## v           -0.002394   0.109092  -0.022   0.9825    
## trt          0.049481   0.192716   0.257   0.7974    
## v:trt       -0.307159   0.150452  -2.042   0.0412 *  
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Correlation of Fixed Effects:
##       (Intr) v      trt   
## v      0.016              
## trt   -0.725 -0.017       
## v:trt -0.018 -0.709  0.030

#  Fit using adaptive Gaussian quadrature to implement ML -- this is
#  only implemented in the case of a single scalar random effect for the intercept

seize.glmm.AGQ.1 <- glmer(seize ~ offset(logo) + v + trt + v:trt + (1 | subj),
     family=poisson, nAGQ=25, data=thedat, na.action=na.omit)

## > summary(seize.glmm.AGQ.1)
## Generalized linear mixed model fit by maximum likelihood (Adaptive
##   Gauss-Hermite Quadrature, nAGQ = 25) [glmerMod]
##  Family: poisson  ( log )
## Formula: seize ~ offset(logo) + v + trt + v:trt + (1 | subj)
##    Data: thedat

##      AIC      BIC   logLik deviance df.resid 
##    970.1    988.5   -480.0    960.1      290 

## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.8660 -0.9376 -0.1677  0.5852  9.8779 

## Random effects:
##  Groups Name        Variance Std.Dev.
##  subj   (Intercept) 0.609    0.7804  
## Number of obs: 295, groups:  subj, 59

## Fixed effects:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  1.03259    0.15269   6.763 1.35e-11 ***
## v            0.11080    0.04689   2.363   0.0181 *  
## trt         -0.02387    0.21067  -0.113   0.9098    
## v:trt       -0.10368    0.06505  -1.594   0.1110    
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Correlation of Fixed Effects:
##       (Intr) v      trt   
## v     -0.162              
## trt   -0.724  0.117       
## v:trt  0.117 -0.721 -0.159

#  For comparision fit the same model using the Laplace approximation

seize.glmm.laplace.1 <- glmer(seize ~ offset(logo) + v + trt + v:trt + (1 | subj),
     family=poisson(link="log"), nAGQ=1, data=thedat, na.action=na.omit)

## > summary(seize.glmm.laplace.1)
## Generalized linear mixed model fit by maximum likelihood (Laplace
##   Approximation) [glmerMod]
##  Family: poisson  ( log )
## Formula: seize ~ offset(logo) + v + trt + v:trt + (1 | subj)
##    Data: thedat

##      AIC      BIC   logLik deviance df.resid 
##   2031.4   2049.8  -1010.7   2021.4      290 

## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.8659 -0.9377 -0.1678  0.5853  9.8781 

## Random effects:
##  Groups Name        Variance Std.Dev.
##  subj   (Intercept) 0.6083   0.7799  
## Number of obs: 295, groups:  subj, 59

## Fixed effects:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  1.03265    0.15254   6.769 1.29e-11 ***
## v            0.11080    0.04672   2.371   0.0177 *  
## trt         -0.02385    0.21047  -0.113   0.9098    
## v:trt       -0.10368    0.06482  -1.599   0.1097    
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Correlation of Fixed Effects:
##       (Intr) v      trt   
## v     -0.162              
## trt   -0.724  0.117       
## v:trt  0.116 -0.721 -0.159

