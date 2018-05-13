#########################################################################
##
##   Homework 5, arthritis trial, using glmer
##
##   glmer() can do quadrature only for models with a single scalar
##   random effect
##   
#########################################################################

library(lme4)

##  Read in data in long format
                      
thedat <- read.table("arthritis.dat")
thedat <- thedat[,c(1:4,6)]
colnames(thedat) <- c("id","trt","age","month","dscale")

##  Random effect for intercept, Laplace approximation - this agrees
##  with the fit using proc glimmix in SAS

arth.glmer.1.laplace <- glmer(dscale ~ month + trt:month + (1 | id),
     family=binomial(link="logit"), nAGQ=1, data=thedat, na.action=na.omit)

## > summary(arth.glmer.1.laplace)
## Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) [
## glmerMod]
##  Family: binomial  ( logit )
## Formula: dscale ~ month + trt:month + (1 | id)
##    Data: thedat

##      AIC      BIC   logLik deviance df.resid 
##   1151.0   1171.2   -571.5   1143.0     1156 

## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.4019 -0.3436 -0.2684 -0.1637  2.9887 

## Random effects:
##  Groups Name        Variance Std.Dev.
##  id     (Intercept) 3.376    1.837   
## Number of obs: 1160, groups:  id, 290

## Fixed effects:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -1.32304    0.19477  -6.793  1.1e-11 ***
## month       -0.08605    0.04716  -1.825   0.0681 .  
## month:trt   -0.16112    0.06301  -2.557   0.0106 *  
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Correlation of Fixed Effects:
##           (Intr) month 
## month     -0.376       
## month:trt -0.004 -0.565

##  Random effect for intercept, adaptive quadrature approximation, 10
##  quadrature points - this appears to agree with that from proc glimmix

arth.glmer.1.quad10 <- glmer(dscale ~ month + trt:month + (1 | id),
     family=binomial(link="logit"), nAGQ=10, data=thedat, na.action=na.omit)

## > summary(arth.glmer.1.quad10)
## Generalized linear mixed model fit by maximum likelihood (Adaptive Gauss-Hermite
##   Quadrature, nAGQ = 10) [glmerMod]
##  Family: binomial  ( logit )
## Formula: dscale ~ month + trt:month + (1 | id)
##    Data: thedat

##      AIC      BIC   logLik deviance df.resid 
##   1142.7   1162.9   -567.3   1134.7     1156 

## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.4454 -0.3376 -0.2627 -0.1591  2.9736 

## Random effects:
##  Groups Name        Variance Std.Dev.
##  id     (Intercept) 3.744    1.935   
## Number of obs: 1160, groups:  id, 290

## Fixed effects:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -1.30435    0.19480  -6.696 2.14e-11 ***
## month       -0.08974    0.04776  -1.879   0.0603 .  
## month:trt   -0.16099    0.06437  -2.501   0.0124 *  
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Correlation of Fixed Effects:
##           (Intr) month 
## month     -0.379       
## month:trt -0.008 -0.567

##  Random effect for intercept, adaptive quadrature approximation, 20
##  quadrature points - this appears to agree with that from proc glimmix

arth.glmer.2.quad20 <- glmer(dscale ~ month + trt:month + (1 | id),
     family=binomial(link="logit"), nAGQ=20, data=thedat, na.action=na.omit)

## > summary(arth.glmer.2.quad20)

## Generalized linear mixed model fit by maximum likelihood (Adaptive Gauss-Hermite
##   Quadrature, nAGQ = 20) [glmerMod]
##  Family: binomial  ( logit )
## Formula: dscale ~ month + trt:month + (1 | id)
##    Data: thedat

##      AIC      BIC   logLik deviance df.resid 
##   1142.6   1162.9   -567.3   1134.6     1156 

## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.4461 -0.3373 -0.2625 -0.1590  2.9735 

## Random effects:
##  Groups Name        Variance Std.Dev.
##  id     (Intercept) 3.753    1.937   
## Number of obs: 1160, groups:  id, 290

## Fixed effects:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -1.30553    0.19516  -6.690 2.24e-11 ***
## month       -0.08978    0.04777  -1.879   0.0602 .  
## month:trt   -0.16096    0.06440  -2.499   0.0124 *  
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Correlation of Fixed Effects:
##           (Intr) month 
## month     -0.378       
## month:trt -0.008 -0.567

##  Random effects for intercept and time, Laplace approximation.  This
##  seems to be close to but not the same as the results from proc
##  glimmix in SAS.  Note that the D matrix has correlation of 1.

arth.glmer.2.laplace <- glmer(dscale ~ month + trt:month + (1 + month | id),
     family=binomial(link="logit"), nAGQ=1, data=thedat, na.action=na.omit)

## > summary(arth.glmer.2.laplace)
## Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) [
## glmerMod]
##  Family: binomial  ( logit )
## Formula: dscale ~ month + trt:month + (1 + month | id)
##    Data: thedat

##      AIC      BIC   logLik deviance df.resid 
##   1143.4   1173.7   -565.7   1131.4     1154 

## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.2996 -0.4784 -0.1920 -0.0741  2.0173 

## Random effects:
##  Groups Name        Variance Std.Dev. Corr
##  id     (Intercept) 1.1273   1.062        
##         month       0.1325   0.364    1.00
## Number of obs: 1160, groups:  id, 290

## Fixed effects:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -0.95927    0.16304  -5.884 4.01e-09 ***
## month       -0.27993    0.08547  -3.275  0.00106 ** 
## month:trt   -0.21283    0.09140  -2.329  0.01988 *  
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Correlation of Fixed Effects:
##           (Intr) month 
## month     -0.342       
## month:trt -0.087 -0.3261

