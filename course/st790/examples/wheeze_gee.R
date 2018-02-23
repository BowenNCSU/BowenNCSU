#########################################################################
#
#   CHAPTER 8, EXAMPLE 6, Six Cities Study
# 
#   Fitting population-averaged models using the the geeglm() 
#   and geese() funtions in geepack.  You will need to install.packages("geepack").
#   The gee() function is not recommended for data with missing
#   values, as it apparently does not have a way to align the data to make
#   sure the correlation structure is estimated correctly. 
#
#   These functions are very similar; geese() seems to be more well-behaved.
#
#   The syntax for each is similar.  A general call looks like
#
#   fit.object <- gee(model formula, id, corstr, family, data, waves, scale.fix=FALSE)
#
#   model formula is a usual R-type model specification for a generalized
#   linear model -- see the documentation for glm()
#
#   corstr is a specification for the working correlation matrix
#
#   family is a specification for the scaled exponential family that
#   would be relevant under independence; the canonical link is the
#   default
#
#   waves is a variable that tells the function to which time point each
#   observation corresponds so that it can align the data properly for correlation
#   matrix estimation when the data are unbalanced.
#
#   scale.fix allows one to fix the scale parameter, which we want to
#   do here for binary data.  For some reason, scale.value=1 causes an error
#   with geeglm() but not witih geese().
#
#   Both of these implementations use moment methods to estimate the  
#   parameters in the working correlation structure.  It is not straightforward
#   to trick another R function into fitting the model with the covariance
#   parameters estimated using the quadratic estimating equation with
#   the Gaussian working structure.  This is possible using the gnls()
#   function in the nlme package, but requires some work on the part
#   of the user, so we don't show this.  
#
#   Robust sandwich standard errors are the default (std.error="san.se")
#
#   We fit a marginal model as discussed in Section 8.8 of the notes.  Given
#   concerns over the observational nature of this study and the
#   likelihood that the covariate smoking status is endogenous, it
#   might be saftest to report the fit with the independence working assumption.
#
###########################################################################

library(geepack)

#  Read in the data - they are in a complicated "wide" format

thedat.wide <- read.table("wheeze.dat",row.names=NULL,na.strings=c("."))
colnames(thedat.wide) <- c("id","city","a9","smk9","w9","a10","smk10","w10",
    "a11","smk11","w11","a12","smk12","w12")

thedat.long <- reshape(thedat.wide,varying=list(c("a9","a10","a11","a12"),
       c("smk9","smk10","smk11","smk12"),c("w9","w10","w11","w12")),                     
       v.names=c("age","smoke","wheeze"),times=9:12,timevar="age",direction="long")
thedat <- thedat.long[order(thedat.long$id),]

#  Fit the model using geese() -- the syntax is similar to gee() and geeglm().
#  Because there are missing values, we must be careful; also we must 
#  constrain the scale parameter to be fixed and = 1

#  Create a data set with the missing observations deleted

thedat.omit <- na.omit(thedat)
cty <- as.factor(thedat.omit$city)
smk <- as.factor(thedat.omit$smoke)

#  specify the reference level of each factor to match SAS

cty <- relevel(cty,ref="portage")
smk <- relevel(smk,ref="2")

#  for family=binomial, the "logit" link is the default

ind.geeglm <- geeglm(wheeze ~ cty+smk,id=id,family=binomial(link="logit"),corstr="independence",
              data=thedat.omit,scale.fix=TRUE,waves=age)
## > summary(ind.geeglm)

## Call:
## geeglm(formula = wheeze ~ cty + smk, family = binomial(link = "logit"), 
##     data = thedat.omit, id = id, waves = age, corstr = "independence", 
##     scale.fix = TRUE)

##  Coefficients:
##             Estimate Std.err Wald Pr(>|W|)
## (Intercept)   -0.456   0.407 1.26     0.26
## ctykingston    0.238   0.509 0.22     0.64
## smk0          -0.449   0.571 0.62     0.43
## smk1          -0.875   0.633 1.91     0.17

## Scale is fixed.

## Correlation: Structure = independenceNumber of clusters:   32   Maximum cluster size: 4 

#  link = identity in output refers to the "correlation link;" see documentation

cs.geeglm <- geeglm(wheeze ~ cty+smk,id=id,family=binomial,corstr="exchangeable",
              data=thedat.omit,scale.fix=TRUE,waves=age)
## > summary(cs.geeglm)
## Call:
## geeglm(formula = wheeze ~ cty + smk, family = binomial, data = thedat.omit, 
##     id = id, waves = age, corstr = "exchangeable", scale.fix = TRUE)

##  Coefficients:
##             Estimate Std.err  Wald Pr(>|W|)
## (Intercept)  -0.4772  0.4476 1.137    0.286
## ctykingston   0.2456  0.4978 0.243    0.622
## smk0         -0.4004  0.5783 0.479    0.489
## smk1         -0.8491  0.6757 1.579    0.209

## Scale is fixed.

## Correlation: Structure = exchangeable  Link = identity 

## Estimated Correlation Parameters:
##       Estimate Std.err
## alpha   0.1255 0.08128
## Number of clusters:   32   Maximum cluster size: 4 

ind.geese <- geese(wheeze ~ cty+smk,id=id,family=binomial,corstr="independence",
                  data=thedat.omit,scale.fix=TRUE,scale.value=1,waves=age)
## > summary(ind.geese)
## Call:
## geese(formula = wheeze ~ cty + smk, id = id, waves = age, data = thedat.omit, 
##     family = binomial, scale.fix = TRUE, scale.value = 1, corstr = "independence")

## Mean Model:
##  Mean Link:                 logit 
##  Variance to Mean Relation: binomial 

##  Coefficients:
##             estimate san.se  wald     p
## (Intercept)   -0.456  0.407 1.256 0.262
## ctykingston    0.238  0.509 0.219 0.640
## smk0          -0.449  0.571 0.620 0.431
## smk1          -0.875  0.633 1.910 0.167

## Scale is fixed.

## Correlation Model:
##  Correlation Structure:     independence 

## Returned Error Value:    0 
## Number of clusters:   32   Maximum cluster size: 4 

cs.geese <- geese(wheeze ~ cty+smk,id=id,family=binomial,corstr="exchangeable",
                  data=thedat.omit,scale.fix=TRUE,scale.value=1,waves=age)
## > summary(cs.geese)
## Call:
## geese(formula = wheeze ~ cty + smk, id = id, waves = age, data = thedat.omit, 
##     family = binomial, scale.fix = TRUE, scale.value = 1, corstr = "exchangeable")

## Mean Model:
##  Mean Link:                 logit 
##  Variance to Mean Relation: binomial 

##  Coefficients:
##               estimate    san.se      wald         p
## (Intercept) -0.4771897 0.4475741 1.1367165 0.2863472
## ctykingston  0.2456336 0.4977888 0.2434923 0.6216951
## smk0        -0.4003875 0.5782689 0.4794037 0.4886925
## smk1        -0.8490938 0.6756644 1.5792444 0.2088696

## Scale is fixed.

## Correlation Model:
##  Correlation Structure:     exchangeable 
##  Correlation Link:          identity 

##  Estimated Correlation Parameters:
##        estimate     san.se    wald         p
## alpha 0.1255118 0.08128246 2.38438 0.1225535

## Returned Error Value:    0 
## Number of clusters:   32   Maximum cluster size: 4

ar1.geeglm <- geeglm(wheeze ~ cty+smk,id=id,family=binomial(link="logit"),corstr="ar1",
              data=thedat.omit,scale.fix=TRUE,waves=age)
## > summary(ar1.geeglm)
## Call:
## geeglm(formula = wheeze ~ cty + smk, family = binomial, data = thedat.omit, 
##     id = id, waves = age, corstr = "ar1", scale.fix = TRUE)

##  Coefficients:
##             Estimate Std.err Wald Pr(>|W|)
## (Intercept)   -0.533   0.465 1.31     0.25
## ctykingston    0.270   0.486 0.31     0.58
## smk0          -0.389   0.585 0.44     0.51
## smk1          -0.709   0.671 1.12     0.29

## Scale is fixed.

## Correlation: Structure = ar1  Link = identity 

## Estimated Correlation Parameters:
##       Estimate Std.err
## alpha    0.242  0.0966
## Number of clusters:   32   Maximum cluster size: 4 

ar1.geese <- geese(wheeze ~ cty+smk,id=id,family=binomial(link="logit"),corstr="ar1",
              data=thedat.omit,scale.fix=TRUE,scale.value=1,waves=age)
## > summary(ar1.geese)
## Call:
## geese(formula = wheeze ~ cty + smk, id = id, waves = age, data = thedat.omit, 
##     family = binomial, scale.fix = TRUE, scale.value = 1, corstr = "ar1")

## Mean Model:
##  Mean Link:                 logit 
##  Variance to Mean Relation: binomial 

##  Coefficients:
##               estimate    san.se      wald         p
## (Intercept) -0.5326111 0.4647326 1.3134519 0.2517706
## ctykingston  0.2703892 0.4856557 0.3099713 0.5776978
## smk0        -0.3889886 0.5848397 0.4423846 0.5059737
## smk1        -0.7094551 0.6705642 1.1193583 0.2900567

## Scale is fixed.

## Correlation Model:
##  Correlation Structure:     ar1 
##  Correlation Link:          identity 

##  Estimated Correlation Parameters:
##        estimate     san.se     wald         p
## alpha 0.2422648 0.09656666 6.293991 0.0121148

## Returned Error Value:    0 
## Number of clusters:   32   Maximum cluster size: 4 

