#########################################################################
#
#   CHAPTER 7, Fitting a logistic regression model using glm()
# 
#   We consider data with binary response (woman has suffered a 
#   myocardial infarction (MI, 0=no, 1=yes), with covariates given 
#   below.
#
#   We demonstrate use of the glm() function for fitting stanard 
#   univariate genearlized linear models
#
###########################################################################

thedat <- read.table("infarc.dat")
colnames(thedat) <- c("id","oral","age","smoke","mi")

infarc.logistic <- glm(mi ~ oral + age + smoke, family=binomial(link="logit"),data=thedat) 

## > summary(infarc.logistic)
## Call:
## glm(formula = mi ~ oral + age + smoke, family = binomial(link = "logit"), 
##     data = thedat)

## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.7019  -0.6045  -0.3215  -0.1562   2.5689  

## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -9.11405    1.75703  -5.187 2.13e-07 ***
## oral         1.97990    0.46968   4.215 2.49e-05 ***
## age          0.16256    0.04454   3.650 0.000262 ***
## smoke        1.81224    0.42938   4.221 2.44e-05 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## (Dispersion parameter for binomial family taken to be 1)

##     Null deviance: 208.20  on 199  degrees of freedom
## Residual deviance: 150.37  on 196  degrees of freedom
## AIC: 158.37

## Number of Fisher Scoring iterations: 5

beta <- infarc.logistic$coef

#  Calculate the odds ratio for smoke and standard error using the delta method

odds <- exp(beta[4])
grad.odds <- exp(beta[4])
se.deltamethod <- sqrt(grad.odds %*% vcov(infarc.logistic)[4,4] %*% grad.odds)

## > c(odds, se.deltamethod)
##    smoke          
## 6.124132 2.629598 
