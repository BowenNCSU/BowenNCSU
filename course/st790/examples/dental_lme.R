#########################################################################
#
#   CHAPTER 5, EXAMPLE 1, Dental Study 
# 
#   Fitting linear mixed effects models uing the lme() function in nlme
#
#   The nlme package is built in to some installations of R; if not,
#   you can install.packages("nlme")
# 
#   A general call to lme() looks like
#
#   fit.object <- lme(model formula, random, correlation, weights, data)
#
#   model formula is a usual R-type model formula specification 
#
#   correlation is a specification for the within-individual correlation
#   structure; the variable on the right specifies the factor
#   determining sets of observations that are assumed
#   independent/uncorrelated (observations are independent by 
#   "child" here).  In general, the specification is 
#   correlation = corTYPE(form = ~ 1 | groupvar).  Unfortunately, as
#   for gls(), it does not seem possible to have groupvar create separate
#   correlation structures for different groups (like "gender") also.
#
#   weights is a specification for the nature of the within-individual
#   variances; the variable on the right specifies a feature by which 
#   they can be different.  Here, using "age" specifies that the variances
#   on the diagonal of the overall covariance matrix can be different
#   across age (time); using "gender" specifies a different variance
#   for each gender (common for all times); and using "gender*age"
#   gives variances that change over age and are different between genders.
#   In gnenerl, the weights option allows one to make the variance on the diagonal
#   of the overall covariance model be different depending on a group
#   factor by weights =  varIdent(form = ~ 1 | groupvar ) or to change
#   over time by weights = varIdent(form = ~1 | timevar )
#
#   Visit https://stat.ethz.ch/R-manual/R-devel/library/nlme/html/lme.html
#   for more info and visit the corClasses and corStruct link for lists
#   of built-in correlation structures.  
#
#   Here, we fit the same models we fit in SAS except for the ones
#   that allow D and the within-child correlation structure to differ
#   by gender (because we can't).
#
#########################################################################

library(nlme)
library(magic)
library(Matrix)

#  Read in the data -- they are in the "long" form required by lme()

thedat <- read.table("dental.dat")
thedat <- thedat[,2:5]      #  remove the first column
colnames(thedat) <- c("id","age","distance","sex")

#  Total number of individuals

m <- max(thedat$id)

#  Create factor variables for use in lme()

child <- factor(thedat$id)
gender <- factor(thedat$sex)

#  Assume separate intercept and slope by gender and fit different
#  covariance structures using ML; REML is the default, so we have to 
#  add option method="ML"

#  First do OLS fit ignoring correlation

dental.ols <- lm(distance ~ -1 + gender + age:gender,data=thedat)

## > coef(dental.ols)
##     gender0     gender1 gender0:age gender1:age 
##  17.3727273  16.3406250   0.4795455   0.7843750 

#  (a) Common D matrix for both genders, default diagonal within-child                               
#      covariance matrix R_i with same variance sigma^2 for each
#      gender.  Thus, sigma^2 is the sum of realization process and
#      measurement error variance

dental.lme.a <- lme(distance ~ -1 + gender + age:gender,data=thedat,
              random = ~ age | child,method="ML")
beta.a <- fixed.effects(dental.lme.a)  #  beta, also fixef(dental.lme.a)
b.a <- random.effects(dental.lme.a)    #  posterior modes bi, also ranef(dental.lme.a)
sebeta.a <- summary(dental.lme.a)$tTable[,"Std.Error"]   #  these SEs
                                        #  are "off" by a factor very
                                        #  close to 1 
D.a <- getVarCov(dental.lme.a, type="random.effects")   #  D
sigma2.a <- dental.lme.a$sigma^2  # sigma^2
V.a <- getVarCov(dental.lme.a,type="marginal",individual=1)   # V_i
R.a <- getVarCov(dental.lme.a,type="conditional",individual=1)   # R_i

#  The attribute varFix returns the CORRECT model-based covariance matrix!
#  Thus, the square roots of its diagonal elements are model-based
#  standard errors (compare to SAS)

sebeta.model.a <-  sqrt(diag(dental.lme.a$varFix))

## > beta.a
##     gender0     gender1 gender0:age gender1:age 
##  17.3727273  16.3406250   0.4795455   0.7843750 
## > sebeta.a
##     gender0     gender1 gender0:age gender1:age 
##  1.20454038  0.99875212  0.10170510  0.08432942
## > sebeta.model.a
##     gender0     gender1 gender0:age gender1:age 
##  1.18202362  0.98008221  0.09980390  0.08275303 
## > sigma2.a
## [1] 1.716204
## > R.a
## child 1 
## Conditional variance covariance matrix
##        1      2      3      4
## 1 1.7162 0.0000 0.0000 0.0000
## 2 0.0000 1.7162 0.0000 0.0000
## 3 0.0000 0.0000 1.7162 0.0000
## 4 0.0000 0.0000 0.0000 1.7162
##   Standard Deviations: 1.31 1.31 1.31 1.31 
## > D.a
## Random effects variance covariance matrix
##             (Intercept)       age
## (Intercept)     4.55690 -0.198250
## age            -0.19825  0.023759
##   Standard Deviations: 2.1347 0.15414
## > V.a
## child 1 
## Marginal variance covariance matrix
##        1      2      3      4
## 1 4.6216 2.8891 2.8727 2.8563
## 2 2.8891 4.6839 3.0464 3.1251
## 3 2.8727 3.0464 4.9363 3.3938
## 4 2.8563 3.1251 3.3938 5.3788
##   Standard Deviations: 2.1498 2.1642 2.2218 2.3192 
## > random.effects(dental.lme.a)
##    (Intercept)          age
## 1  -0.68280789 -0.039971593
## 2  -0.45919650  0.071880564
## 3  -0.03101916  0.093013749
## 4   1.61181879  0.030833688
## 5   0.43846143 -0.042996959 ...

#  (b) Common D matrix for both genders, diagonal within-child                               
#      covariance matrix R_i with different variance for each gender

dental.lme.b <- lme(distance ~ -1 + gender + age:gender,data=thedat,
                    random = ~ age | child,weights = varIdent(form = ~ 1 | gender),
                    method="ML")
beta.b <- fixed.effects(dental.lme.b)  #  beta
sebeta.model.b <-  sqrt(diag(dental.lme.b$varFix))
b.b <- random.effects(dental.lme.b)    #  posterior modes bi
D.b <- getVarCov(dental.lme.b, type="random.effects")   #  D
R.b.1 <- getVarCov(dental.lme.b,type="conditional",individual=1)   # R_1
R.b.12 <- getVarCov(dental.lme.b,type="conditional",individual=12)   # R_12

#  As with gls(), when a weight statement is used, the standard
#  deviation for the first group is sigma, and those for other groups
#  are parameterized as this standard deviation x a factor. To get 
#  sigma^2_B, we must extract that factor from the varStruct object:
#
## > dental.lme.b$modelStruct$varStruct
## Variance function structure of class varIdent representing
##        0        1 
## 1.000000 2.431012  

sigma2vec.b <-
    matrix((1/unique(attributes(dental.lme.b$modelStruct$varStruct)$weights)*dental.lme.b$sigma)^2,nrow=1,byrow=TRUE)
colnames(sigma2vec.b) <- c("sigma2.b.G","sigma2.b.B")

V.b.1 <- getVarCov(dental.lme.b,type="marginal",individual=1)   # V_1 (girl)
V.b.12 <- getVarCov(dental.lme.b,type="marginal",individual=12)   # V_12 (boy)

## > beta.b
##     gender0     gender1 gender0:age gender1:age 
##  17.3727273  16.3406250   0.4795455   0.7843750 
## > sebeta.model.b
##     gender0     gender1 gender0:age gender1:age 
##  0.73864469  1.11140546  0.06179965  0.09722184 
## > D.b
## Random effects variance covariance matrix
##             (Intercept)       age
## (Intercept)     3.19860 -0.110360
## age            -0.11036  0.019766
##   Standard Deviations: 1.7885 0.14059 
## > sigma2vec.b
##      sigma2.b.G sigma2.b.B
## [1,]  0.4449132   2.629357
## > R.b.1
## child 1 
## Conditional variance covariance matrix
##         1       2       3       4
## 1 0.44491 0.00000 0.00000 0.00000
## 2 0.00000 0.44491 0.00000 0.00000
## 3 0.00000 0.00000 0.44491 0.00000
## 4 0.00000 0.00000 0.00000 0.44491
##   Standard Deviations: 0.66702 0.66702 0.66702 0.66702 
## > R.b.12
## child 12 
## Conditional variance covariance matrix
##        1      2      3      4
## 1 2.6294 0.0000 0.0000 0.0000
## 2 0.0000 2.6294 0.0000 0.0000
## 3 0.0000 0.0000 2.6294 0.0000
## 4 0.0000 0.0000 0.0000 2.6294
##   Standard Deviations: 1.6215 1.6215 1.6215 1.6215 
## > V.b.1
## child 1 
## Marginal variance covariance matrix
##        1      2      3      4
## 1 3.1428 2.7934 2.8889 2.9845
## 2 2.7934 3.4129 3.1426 3.3172
## 3 2.8889 3.1426 3.8412 3.6499
## 4 2.9845 3.3172 3.6499 4.4275
##   Standard Deviations: 1.7728 1.8474 1.9599 2.1042 
## > V.b.12
## child 12 
## Marginal variance covariance matrix
##        1      2      3      4
## 1 5.3272 2.7934 2.8889 2.9845
## 2 2.7934 5.5973 3.1426 3.3172
## 3 2.8889 3.1426 6.0256 3.6499
## 4 2.9845 3.3172 3.6499 6.6120
##   Standard Deviations: 2.3081 2.3659 2.4547 2.5714 

#  (c) Cannot fit a model with separate D matrices for each gender

#  (d)' We can't fit model (d), which specifies common within-child 
#      AR(1) and common-within child measurement error, because lme()
#      does not support this; this can only be done with spatial   
#      correlation structures like the exponential or Gaussian (see (e) below).
#      So instead, we illustrate specifying something other than a diagonal
#      R_i matrix by fitting a model with Common D matrix for both genders, 
#      common within-child AR(1), and no within-child measurement error

dental.lme.d <- lme(distance ~ -1 + gender + age:gender,data=thedat,
                    random = ~ age | child,
                    correlation=corAR1(form = ~ age | child),
                    method="ML")
beta.d <- fixed.effects(dental.lme.d)
b.d <- random.effects(dental.lme.d)
sebeta.d <- sqrt(diag(dental.lme.d$varFix))
D.d <- getVarCov(dental.lme.d, type="random.effects")   #  D
sigma2.d <- dental.lme.d$sigma^2  # sigma^2
V.d <- getVarCov(dental.lme.d,type="marginal",individual=1)   # V_i
R.d <- getVarCov(dental.lme.d,type="conditional",individual=1)   # R_i

#  The correlation parameter alpha, which is constrained to be |alpha|<=1,
#  is estimated to be 0 (= Phi1).  So this fit is the same as (a)

## > summary(dental.lme.d) ...
## Correlation Structure: ARMA(1,0)
##  Formula: ~age | child 
##  Parameter estimate(s):
## Phi1 
##    0 

#  Thus, we get a diagonal R_i

## > R.d
## child 1 
## Conditional variance covariance matrix
##        1      2      3      4
## 1 1.7162 0.0000 0.0000 0.0000
## 2 0.0000 1.7162 0.0000 0.0000
## 3 0.0000 0.0000 1.7162 0.0000
## 4 0.0000 0.0000 0.0000 1.7162
##   Standard Deviations: 1.31 1.31 1.31 1.31 

#  (e) Common D matrix for both genders, common within-child AR(1) 
#     (exponential correlation) and common within-child measurement error

dental.lme.e <- lme(distance ~ -1 + gender + age:gender,data=thedat,
                    random = ~ age | child,
                    correlation=corExp(form = ~ age | child , nugget=TRUE),
                    method="ML")
beta.e <- fixed.effects(dental.lme.e)
b.e <- random.effects(dental.lme.e)
sebeta.model.e <-  sqrt(diag(dental.lme.e$varFix))
D.e <- getVarCov(dental.lme.e, type="random.effects")   #  D
sigma2.e <- dental.lme.e$sigma^2  # sigma^2
R.e <- getVarCov(dental.lme.e,type="conditional",individual=1)  #  R_i
V.e <- getVarCov(dental.lme.e,type="marginal",individual=1)   # V_i

## > R.e
## child 1 
## Conditional variance covariance matrix
##            1          2          3          4
## 1 1.7162e+00 6.9306e-07 3.1672e-13 1.4473e-19
## 2 6.9306e-07 1.7162e+00 6.9306e-07 3.1672e-13
## 3 3.1672e-13 6.9306e-07 1.7162e+00 6.9306e-07
## 4 1.4473e-19 3.1672e-13 6.9306e-07 1.7162e+00
##   Standard Deviations: 1.31 1.31 1.31 1.31 

#  This is the associated correlation matrix

Rcorr.e <- cov2cor(simplify2array(getVarCov(dental.lme.e,type="conditional",individual=1)[[1]]))

#  The "range" and "nugget" parameters can be extracted from the fit.
#  The nugget parameter is equal to, in our notation, sigma^2_M/(sigma^2_P+sigma^2_M)

range.e <- exp(coef(dental.lme.e$modelStruct$corStruct))[1]
nugget.e <- exp(coef(dental.lme.e$modelStruct$corStruct))[2]

#  The implied correlations are then as follows:

## > c((1-nugget.e)*exp(-2/range.e),(1-nugget.e)*exp(-4/range.e),(1-nugget.e)*exp(-6/range.e))
##       nugget       nugget       nugget 
## 3.968420e-07 1.813500e-13 8.287382e-20 

#  The estimated correlation matrix is below, and has slightly
#  different correlations.  The ratio of those above to these is 
#  alomost exactly the same constant = the true model-based
#  standard errors/standard errors lme() reports in the summary()!

## > Rcorr.e
##              1            2            3            4
## 1 1.000000e+00 4.038363e-07 1.845462e-13 8.433446e-20
## 2 4.038363e-07 1.000000e+00 4.038363e-07 1.845462e-13
## 3 1.845462e-13 4.038363e-07 1.000000e+00 4.038363e-07
## 4 8.433446e-20 1.845462e-13 4.038363e-07 1.000000e+00

## > c((1-nugget.e)*exp(-2/range.e),(1-nugget.e)*exp(-4/range.e),(1-nugget.e)*exp(-6/range.e))/c(Rcorr.e[1,2],Rcorr.e[1,3],Rcorr.e[1,4])
##    nugget    nugget    nugget 
## 0.9826804 0.9826804 0.9826804 

#  Compare the fitted models via AIC and BIC; model (b), with
#  diagonal R_i with gender-specific within-child variances and common
#  random effects covariance matrix D is preferred

## > anova(dental.lme.a,dental.lme.b,dental.lme.d,dental.lme.e)
##              Model df      AIC      BIC    logLik   Test   L.Ratio p-value
## dental.lme.a     1  8 443.8060 465.2630 -213.9030                         
## dental.lme.b     2  9 424.0424 448.1816 -203.0212 1 vs 2 21.763562  <.0001
## dental.lme.d     3  9 445.8060 469.9451 -213.9030                         
## dental.lme.e     4 10 447.8060 474.6273 -213.9030 3 vs 4  0.000005  0.9982

#  Refit model (b) using REML and get 

dental.lme.b.reml <- lme(distance ~ -1 + gender + age:gender,data=thedat,
                    random = ~ age | child,weights = varIdent(form = ~ 1 | gender))
beta.b.reml <- fixed.effects(dental.lme.b.reml)  #  beta
sebeta.model.b.reml <-  sqrt(diag(dental.lme.b$varFix))
b.b.reml <- random.effects(dental.lme.b.reml)    #  posterior modes bi
D.b.reml <- getVarCov(dental.lme.b.reml, type="random.effects")   #  D

#  Random effects empirical Bayes estimates bhat_i for first 5 girls

## > b.b.reml[1:20,]
##    (Intercept)         age
## 1  -0.44921832 -0.07172941
## 2  -1.40366050  0.16096191
## 3  -1.07810345  0.19761263
## 4   1.76891772  0.03476773
## 5   1.05432275 -0.09938676
## 6  -0.64509606 -0.07587888
## 7  -0.09327961  0.03995026
## 8   2.16610833 -0.13534303
## 9  -0.12094371 -0.12428354
## 10 -3.09492662 -0.08314475
## 11  1.89587947  0.15647384
## 12  1.35617258  0.08962144
## 13 -0.89703274 -0.03952357
## 14 -0.36228439 -0.02199483
## 15  1.80008233 -0.04457382
## 16 -1.21632557 -0.03814274
## 17  0.95641574  0.01859871
## 18 -0.71791874 -0.02707062
## 19 -0.05077251 -0.08286598
## 20 -0.17798309  0.03011834

#  PA predicted values X_i betahat are produced by level=0; SS
#  predicted values X_i betahat + Z_i bhat_i are produced by level =
#  1; both are gotten by level = 0:1

## > fitted(dental.lme.b.reml,level=0:1)[1:20,] 
##       fixed    child
## 1  21.20909 20.18604
## 2  22.16818 21.00167
## 3  23.12727 21.81730
## 4  24.08636 22.63293
## 5  21.20909 21.09313
## 6  22.16818 22.37414
## 7  23.12727 23.65516
## 8  24.08636 24.93617
## 9  21.20909 21.71189
## 10 22.16818 23.06620
## 11 23.12727 24.42052
## 12 24.08636 25.77484
## 13 21.20909 23.25615
## 14 22.16818 24.28478
## 15 23.12727 25.31340
## 16 24.08636 26.34203
## 17 21.20909 21.46832
## 18 22.16818 22.22864
## 19 23.12727 22.98895
## 20 24.08636 23.74927

#  Can also extract both marginal (PA) residuals 
#  Y_i-X_i betahat and conditional (SS) residuals Y_i-X_i betahat -
#  Z_i bhat_i using the fitted() function -- level = 0 gives PA, level
#  = 1 gives SS (default), and level = 0:1 gives both; type = default
#  is "raw" residuals, type = 

#  PA and SS raw residuals; can get standardized "pearson"
#  residuals as in SAS PROC MIXED with type = option

## > residuals(dental.lme.b.reml,level=0:1)[1:20,] 
##         fixed       child
## 1  -0.2090909  0.81396273
## 2  -2.1681818 -1.00166935
## 3  -1.6272727 -0.31730143
## 4  -1.0863636  0.36706649
## 5  -0.2090909 -0.09312569
## 6  -0.6681818 -0.87414042
## 7   0.8727273  0.34484485
## 8   1.4136364  0.56383011
## 9  -0.7090909 -1.21188851
## 10  1.8318182  0.93379531
## 11  1.3727273  0.07947914
## 12  1.9136364  0.22516297
## 13  2.2909091  0.24384949
## 14  2.3318182  0.21522311
## 15  1.8727273 -0.31340327
## 16  2.4136364  0.15797036
## 17  0.2909091  0.03168042
## 18  0.8318182  0.77136303
## 19 -0.6272727 -0.48895436
## 20 -0.5863636 -0.24927175

#   PA and SS pearson residuals -- the PA residuals differ from those
#   obtained from PROC MIXED apparently due to different standardization

## > residuals(dental.lme.b.reml,level=0:1,type="pearson")[1:20,] 
##         fixed       child
## 1  -0.3138325  1.22170756
## 2  -3.2543064 -1.50344357
## 3  -2.4424354 -0.47624977
## 4  -1.6305644  0.55094403
## 5  -0.3138325 -0.13977589
## 6  -1.0028995 -1.31203055
## 7   1.3099095  0.51759072
## 8   2.1217805  0.84627403
## 9  -1.0643015 -1.81896949
## 10  2.7494455  1.40156884
## 11  2.0603785  0.11929326
## 12  2.8722495  0.33795565
## 13  3.4385124  0.36600296
## 14  3.4999144  0.32303654
## 15  2.8108475 -0.47039886
## 16  3.6227185  0.23710371
## 17  0.4366365  0.04755035
## 18  1.2485075  1.15776806
## 19 -0.9414975 -0.73389016
## 20 -0.8800954 -0.37414143

#  Autocorrelation function estimated from pearson SS residuals

## >  ACF(dental.lme.b.reml)
##   lag         ACF
## 1   0  1.00000000
## 2   1 -0.43568449
## 3   2 -0.07302713
## 4   3 -0.27813586

#  Estimated semivariogram based on pearson SS residuals.  An
#  approximately flat pattern reflects lack of serial correlation
#  (one can also plot this; we don't given it has only 3 lags)

## >  Variogram(dental.lme.b.reml, form = ~ age)
##      variog dist n.pairs
## 1 1.0219518    2      81
## 2 0.7441608    4      54
## 3 0.8075591    6      27

#  Plot SS residuals vs. predicted values

pdf("dental.residplot.pdf",width=10)
plot(dental.lme.b.reml, resid(. , type="p",level=1) ~ fitted(.,level=1) )
graphics.off()

#  QQ plot of SS residuals 

pdf("dental.qqplot.pdf",width=10)
qqnorm(dental.lme.b.reml, ~ resid(. , type="p",level=1),abline=c(0,1))
graphics.off()

#  One can also make QQ plots and histograms of the bhat_i themselves
#  to assess the normality of the random effects, but remember that 
#  these are "shrunken" so could be misleading.  

pdf("dental.qqplot.raneff.pdf",width=10)
qqnorm(dental.lme.b.reml, ~ ranef(.))
graphics.off()

#  histograms

pdf("dental.histo.raneff.pdf",width=10)
par(mfrow=c(1,2))
hist(b.b.reml[,1],xlab="Intercept Random Effect",main="Empirical Bayes Intercepts",freq=FALSE)
hist(b.b.reml[,2],xlab="Slope Random Effect",main="Empirical Bayes Slopes",freq=FALSE)
graphics.off()

