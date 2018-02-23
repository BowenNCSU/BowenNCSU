#########################################################################
#
#   CHAPTER 5, EXAMPLE 1, Dental Study 
# 
#   Fitting linear mixed effects models uing the lmer() function in
#   the lme4 package; you can install.packages("lme4")
#
#   lme4 is apparently more computationally efficient than nlme().
#   It also implements generalized linear mixed effects models in the
#   glmer() function and artibrary nonlinear mixed effects models in
#   the nlmer() function (nlme has the nlme() function for the latter;
#   we discuss these models in Chapter 9).  However, it does not
#   provide a way to model within-individual variances or
#   within-individual correlation structures, so it is less flexible
#   for modeling.
# 
#   Here, we demonstrate its basic use with the dental data.  A basic
#   call looks like
#  
#   fit.object <- lmer(model formula,REML,...)
#
#   The key difference is the way the model is specified; instead of
#   there being a separate "random" statement, random effects
#   structure is specified as part of the formula, in parentheses, as
#   demonstrated below.  
#
#########################################################################

library(lme4)

#  Read in the data -- they are in the "long" form required by lme()

thedat <- read.table("dental.dat")
thedat <- thedat[,2:5]      #  remove the first column
colnames(thedat) <- c("id","age","distance","sex")

#  Total number of individuals

m <- max(thedat$id)

#  Create factor variables for use in lmer()

child <- factor(thedat$id)
gender <- factor(thedat$sex)

#  We fit only model (a), Common D matrix for both genders, default
#  diagonal wintin-child covariance matrix  R_i witih same variance
#  sigma^2 for each gender.  We use ML as in the lme() and SAS programs.

#  The random effects structure is specified in parentheses -- here,
#  we allow for random intercept and slope that are correlated

dental.lmer.a <- lmer(distance ~ -1 + gender + age:gender + (1 + age | child),
                      REML=FALSE,data=thedat)

#  Note that the model-based standard errors are correct

## > summary(dental.lmer.a)
## Linear mixed model fit by maximum likelihood  ['lmerMod']
## Formula: distance ~ -1 + gender + age:gender + (1 + age | child)
##    Data: thedat

##      AIC      BIC   logLik deviance df.resid 
##    443.8    465.3   -213.9    427.8      100 

## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.3360 -0.4154  0.0104  0.4917  3.8582 

## Random effects:
##  Groups   Name        Variance Std.Dev. Corr 
##  child    (Intercept) 4.55691  2.1347        
##           age         0.02376  0.1541   -0.60
##  Residual             1.71620  1.3100        
## Number of obs: 108, groups:  child, 27

## Fixed effects:
##             Estimate Std. Error t value
## gender0     17.37273    1.18202  14.697
## gender1     16.34063    0.98008  16.673
## gender0:age  0.47955    0.09980   4.805
## gender1:age  0.78437    0.08275   9.479

## Correlation of Fixed Effects:
##             gendr0 gendr1 gndr0:
## gender1      0.000              
## gender0:age -0.880  0.000       
## gender1:age  0.000 -0.880  0.000
## > sebeta.model.a
##     gender0     gender1 gender0:age gender1:age 
##  1.18202362  0.98008221  0.09980390  0.08275303 

beta.lmer.a <- fixef(dental.lmer.a)
b.lmer.a <- ranef(dental.lmer.a)
sigma2.lmer.a <- sigma(dental.lmer.a)^2

#  It is pretty unwieldy to extract the covariance matrix D of the
#  random effects.  We can look at the variances and correlation with

vc.a <- VarCorr(dental.lmer.a)

## > print(vc.a,comp="Variance")
##  Groups   Name        Variance Corr  
##  child    (Intercept) 4.556907       
##           age         0.023759 -0.603
##  Residual             1.716204   

#  All the covariance matrix stuff can be put in a data frame, from
#  which it can be extracted to form the matrices D, R_i, and V_i

vc.da <- as.dataframe(vc.a,order="lower.tri")

## > vc.da
##        grp        var1 var2        vcov      sdcor
## 1    child (Intercept) <NA>  4.55690750  2.1346914
## 2    child (Intercept)  age -0.19825328 -0.6025211
## 3    child         age <NA>  0.02375889  0.1541392
## 4 Residual        <NA> <NA>  1.71620388  1.3100396

#  One could presumably write a function to do this for arbitrary
#  numbers of random effects knowing how they are stored in such a
#  data frame; I am not doing that

D.lmer.a <- matrix(c(vc.da[1,4],vc.da[2,4],vc.da[2,4],vc[3,4]),2,2,byrow=TRUE)

## > D.lmer.a
##            [,1]        [,2]
## [1,]  4.5569075 -0.19825328
## [2,] -0.1982533  0.02375889

