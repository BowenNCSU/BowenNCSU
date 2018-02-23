#########################################################################
#
#   CHAPTER 5, EXAMPLE 1, Dental Study 
# 
#   Fitting population-averaged models uing the gls() function in nlme
#
#   The nlme package is built in to some installations of R; if not,
#   you can install.packages(nlme)
# 
#   A general call to gls() looks like
#
#   fit.object <- gls(model formula, correlation, weights, data)
#
#   model formula is a usual R-type model specification 
#
#   correlation is a specification for the assumed overall correlation
#   structure; the variable on the right specifies the factor
#   determining sets of observations that are assumed
#   independent/uncorrelated (observations are independent by 
#   "child" here).  In general, the specification is 
#   correlation = corTYPE(form = ~ 1 | groupvar).  Unfortunately
#   it does not seem possible to have groupvar create separate
#   correlation structures for different groups (like "gender") also.
#
#   weights is a specification for the nature of the variances; the
#   variable on the right specifies a feature by which the variance
#   can be different.  Here, using "age" specifies that the variances
#   on the diagonal of the overall covariance matrix can be different
#   across age (time); using "gender" specifies a different variance
#   for each gender (common for all times); and using "gender*age"
#   gives variances that change over age and are different between genders.
#   In gnenerl, the weights option allows one to make the variance on the diagonal
#   of the overall covariance model be different depending on a group
#   factor by weights =  varIdent(form = ~ 1 | groupvar ) or to change
#   over time by weights = varIdent(form = ~1 | timevar )
#
#   Visit https://stat.ethz.ch/R-manual/R-devel/library/nlme/html/gls.html
#   for more info and visit the corClasses and corStruct link for lists
#   of built-in correlation structures.  As above, gls() DOES NOT allow
#   correlation structure to differ by levels of a group variable like
#   gender - one would have to write his/her own corStruct class to do 
#   this.  We don't do this here and restrict attention to specifications
#   that are the same for each group.
#
#   WARNING:  There is a MISTAKE in gls(), and it DOES NOT calculate
#   the model-based sampling covariance matrix of betahat correctly!
#   Thus, the model-based standard errors that gls() reports are ever
#   so slightly "off" from the correct values.  MOREOVER, gls() does
#   not calculate the robust sandwich covariance matrix, so does not
#   offer the option of getting robust (or empirical) standard
#   errors.  The function robust.cov below will compute correct
#   model-based standard errors along with robust sandwich standard
#   errors.  
#   
#########################################################################

#  Function to compute correct and robust standard errors from a gls()
#  model object

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

library(nlme)

#  Read in the data -- they are in the "long" form required by gls()

thedat <- read.table("dental.dat")
thedat <- thedat[,2:5]      #  remove the first column
colnames(thedat) <- c("id","age","distance","sex")

#  Total number of individuals

m <- max(thedat$id)

#  Create factor variables for use in gls()

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

#  (a) Common unstructured correlation with variances changing over time
#  for both genders -- we extract components of the fit.  Note that
#  gls() defines BIC dfferently from SAS (it uses the total number of
#  observations N while SAS MIXED uses the total number of individuals m                                      
#  The weights statement makes the variances on the diagonal differ over
#  time - the default with no weight statement is that they are the same
#  for all times

dental.un <- gls(distance ~ -1 + gender + age:gender,data=thedat,
                 correlation=corSymm(form =  ~ 1 | child),
                 weights = varIdent(form = ~ 1 | age),method="ML")
beta.un <- coef(dental.un)
sebeta.un <- summary(dental.un)$tTable[,"Std.Error"]
V.un <- getVarCov(dental.un)  
Gamma.un <- cov2cor(V.un)
vars.un <- diag(V.un)

#  Alternatively, can get the correlation and covariance matrix for
#  individual 1 as
#
#  Gamma.un <- corMatrix(dental.un$modelStruct$corStruct)[[1]]
#  parms.un <- coef(dental.un$model,unconstrained=FALSE)   
#  > parms.un
#  corStruct1   corStruct2   corStruct3   corStruct4   corStruct5   corStruct6 
#   0.5443349    0.6525649    0.5187568    0.5607205    0.7190314    0.7275871 
#  varStruct.10 varStruct.12 varStruct.14 
#   0.8759573    1.0807973    0.9497922 
#  > dental.un$sigma
#  [1] 2.262555
#  - the last 3 entries are multipliers of dental.un$sigma that can be
#  used to get the variances at each age as follows
#  vars.un <- ( c(1,parms.un[7:9])*dental.un$sigma)^2
#  Then get the covariance matrix as
#  V.un - diag(sqrt(vars.un))%*%Gamma.un%*%diag(sqrt(vars.un))

#  Call robust.cov to get corrected model-based SEs and compare to
#  those from gls().  Compare the results to those from SAS proc
#  mixed, which DOES compute the correct standard errors

u <- dental.un
robust.un <- robust.cov(u,m)
sebeta.un.corrected <- robust.un$se.model

## > rbind(beta.un,sebeta.un)
##             gender0    gender1 gender0:age gender1:age
## beta.un   17.425369 15.8422918  0.47636467  0.82680277
## sebeta.un  1.149863  0.9534164  0.09723211  0.08062061
## ## > sebeta.un.corrected
##     gender0     gender1 gender0:age gender1:age 
##  1.12836878  0.93559397  0.09541453  0.07911355 
## > V.un
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]
## [1,] 5.1192 2.4409 3.6105 2.5223
## [2,] 2.4409 3.9279 2.7175 3.0624
## [3,] 3.6105 2.7175 5.9798 3.8235
## [4,] 2.5223 3.0624 3.8235 4.6180
##   Standard Deviations: 2.2626 1.9819 2.4454 2.149 
## > Gamma.un
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]
## [1,] 1.00000 0.54433 0.65256 0.51876
## [2,] 0.54433 1.00000 0.56072 0.71903
## [3,] 0.65256 0.56072 1.00000 0.72759
## [4,] 0.51876 0.71903 0.72759 1.00000
##   Standard Deviations: 1 1 1 1 

##  In the rest of these fits, we do not calculate the corrected
##  standard errors.  We DO calculate them when we settle on a final
##  model below.

#  (c) Common compound symmetry 

dental.cs <- gls(distance ~ -1 + gender + age:gender,data=thedat,
                 correlation=corCompSymm(form = ~ 1 | child),method="ML")
beta.cs <- coef(dental.cs)
sebeta.cs <- summary(dental.cs)$tTable[,"Std.Error"]
V.cs <- getVarCov(dental.cs)  
Gamma.cs <- cov2cor(V.cs)     
vars.cs <- diag(V.cs)

#  Alternatively, 
#  Gamma.cs <- corMatrix(dental.cs$modelStruct$corStruct)[[1]]
#  V.cs <- dental.cs$sigma^2*Gamma.cs

#  Note that the OLS estimator and the estimator for beta obtained
#  from assuming any kind of compound symmetric structure are
#  NUMERICALLY IDENTICAL. This is no accident, as you'll learn in a
#  future homework.  Of course, the STANDARD ERRORS are different,
#  reflecting the assumed covariance structure.  

## > rbind(beta.cs,sebeta.cs)
##            gender0   gender1 gender0:age gender1:age
## beta.cs   17.37273 16.340625  0.47954545  0.78437500
## sebeta.cs  1.18365  0.981431  0.09406714  0.07799635
## > V.cs
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]
## [1,] 4.9052 3.0306 3.0306 3.0306
## [2,] 3.0306 4.9052 3.0306 3.0306
## [3,] 3.0306 3.0306 4.9052 3.0306
## [4,] 3.0306 3.0306 3.0306 4.9052
##   Standard Deviations: 2.2148 2.2148 2.2148 2.2148 
## > Gamma.cs
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]
## [1,] 1.00000 0.61783 0.61783 0.61783
## [2,] 0.61783 1.00000 0.61783 0.61783
## [3,] 0.61783 0.61783 1.00000 0.61783
## [4,] 0.61783 0.61783 0.61783 1.00000
##   Standard Deviations: 1 1 1 1 

#  (c)'Common compound symmetry with variance different by gender

dental.cs2 <- gls(distance ~ -1 + gender + gender:age,data=thedat,correlation=corCompSymm(form = ~ 1 | child),
                 weights = varIdent(form = ~ 1 | gender),method="ML")
beta.cs2 <- coef(dental.cs2)
sebeta.cs2 <- summary(dental.cs2)$tTable[,"Std.Error"]
V.cs2.girl <-  getVarCov(dental.cs2, individual=1)
V.cs2.boy <-  getVarCov(dental.cs2, individual=12)
Gamma.cs2 <- cov2cor(V.cs2.girl)

## > rbind(beta.cs2,sebeta.cs2)
##               gender0   gender1 gender0:age gender1:age
## beta.cs2   17.3727273 16.340625  0.47954545  0.78437500
## sebeta.cs2  0.8175862  1.163776  0.06125437  0.08719126
## >  V.cs2.girl
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]
## [1,] 2.8677 2.0728 2.0728 2.0728
## [2,] 2.0728 2.8677 2.0728 2.0728
## [3,] 2.0728 2.0728 2.8677 2.0728
## [4,] 2.0728 2.0728 2.0728 2.8677
##   Standard Deviations: 1.6934 1.6934 1.6934 1.6934 
## > V.cs2.boy
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]
## [1,] 8.4514 6.1088 6.1088 6.1088
## [2,] 6.1088 8.4514 6.1088 6.1088
## [3,] 6.1088 6.1088 8.4514 6.1088
## [4,] 6.1088 6.1088 6.1088 8.4514
##   Standard Deviations: 2.9071 2.9071 2.9071 2.9071 
## > Gamma.cs2
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]
## [1,] 1.00000 0.72281 0.72281 0.72281
## [2,] 0.72281 1.00000 0.72281 0.72281
## [3,] 0.72281 0.72281 1.00000 0.72281
## [4,] 0.72281 0.72281 0.72281 1.00000
##   Standard Deviations: 1 1 1 1 

#  (d) Common AR(1)

dental.ar1 <- gls(distance ~ -1 + gender + age:gender,data=thedat,
                  correlation=corAR1(form = ~ 1 | child),method="ML")
beta.ar1 <- coef(dental.ar1)
sebeta.ar1 <- summary(dental.ar1)$tTable[,"Std.Error"]
V.ar1 <- getVarCov(dental.ar1)  #  or corMatrix(dental.un$modelStruct$corStruct)[[1]]
Gamma.ar1 <- cov2cor(V.ar1)
vars.ar1 <- diag(V.ar1)

## > rbind(beta.ar1,sebeta.ar1)
##              gender0   gender1 gender0:age gender1:age
## beta.ar1   17.321720 16.591996   0.4837322   0.7695718
## sebeta.ar1  1.634509  1.355263   0.1409898   0.1169026
## > V.ar1
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]
## [1,] 4.8908 2.9693 1.8027 1.0944
## [2,] 2.9693 4.8908 2.9693 1.8027
## [3,] 1.8027 2.9693 4.8908 2.9693
## [4,] 1.0944 1.8027 2.9693 4.8908
##   Standard Deviations: 2.2115 2.2115 2.2115 2.2115 
## > Gamma.ar1
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]
## [1,] 1.00000 0.60712 0.36859 0.22378
## [2,] 0.60712 1.00000 0.60712 0.36859
## [3,] 0.36859 0.60712 1.00000 0.60712
## [4,] 0.22378 0.36859 0.60712 1.00000
##   Standard Deviations: 1 1 1 1 

#  (d)'  Common AR(1) with variance different by gender

dental.ar12 <- gls(distance ~ -1 + gender + age:gender,data=thedat,
             correlation=corAR1(form = ~ 1 | child),,method="ML",
             weights = varIdent(form = ~ 1 | gender))
beta.ar12<- coef(dental.ar12)
sebeta.ar12 <- summary(dental.ar12)$tTable[,"Std.Error"]
V.ar12.girl <-  getVarCov(dental.ar12, individual=1)
V.ar12.boy <-  getVarCov(dental.ar12, individual=12)
Gamma.ar12 <- cov2cor(V.ar12.girl)

## > rbind(beta.ar12,sebeta.ar12)
##               gender0   gender1 gender0:age gender1:age
## beta.ar12   17.313902 16.642771  0.48429625   0.7675776
## sebeta.ar12  1.130648  1.644926  0.09512683   0.1383955
## > V.ar12.girl
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]
## [1,] 2.8431 2.0824 1.5252 1.1172
## [2,] 2.0824 2.8431 2.0824 1.5252
## [3,] 1.5252 2.0824 2.8431 2.0824
## [4,] 1.1172 1.5252 2.0824 2.8431
##   Standard Deviations: 1.6861 1.6861 1.6861 1.6861 
## > V.ar12.boy
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]
## [1,] 8.7530 6.4111 4.6958 3.4394
## [2,] 6.4111 8.7530 6.4111 4.6958
## [3,] 4.6958 6.4111 8.7530 6.4111
## [4,] 3.4394 4.6958 6.4111 8.7530
##   Standard Deviations: 2.9585 2.9585 2.9585 2.9585 
## > Gamma.ar12
## Marginal variance covariance matrix
##         [,1]    [,2]    [,3]    [,4]
## [1,] 1.00000 0.73244 0.53648 0.39294
## [2,] 0.73244 1.00000 0.73244 0.53648
## [3,] 0.53648 0.73244 1.00000 0.73244
## [4,] 0.39294 0.53648 0.73244 1.00000
##   Standard Deviations: 1 1 1 1 


# (e) Common one-dependent

dental.1d <- gls(distance ~ -1 + gender +  age:gender,data=thedat,
                 correlation=corARMA(form = ~ 1 | child, q=1),method="ML")
beta.1d <- coef(dental.1d)
sebeta.1d <- summary(dental.1d)$tTable[,"Std.Error"]
V.1d <- getVarCov(dental.1d)  #  or corMatrix(dental.un$modelStruct$corStruct)[[1]]
Gamma.1d <- cov2cor(V.1d)
vars.1d <- diag(V.1d)

## > rbind(beta.1d,sebeta.1d)
##             gender0   gender1 gender0:age gender1:age
## beta.1d   17.303464 16.620793   0.4856185   0.7629024
## sebeta.1d  1.741097  1.443641   0.1540519   0.1277331
## > V.1d
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]
## [1,] 4.5294 1.6120 0.0000 0.0000
## [2,] 1.6120 4.5294 1.6120 0.0000
## [3,] 0.0000 1.6120 4.5294 1.6120
## [4,] 0.0000 0.0000 1.6120 4.5294
##   Standard Deviations: 2.1282 2.1282 2.1282 2.1282 
## > Gamma.1d
## Marginal variance covariance matrix
##        [,1]   [,2]   [,3]   [,4]
## [1,] 1.0000 0.3559 0.0000 0.0000
## [2,] 0.3559 1.0000 0.3559 0.0000
## [3,] 0.0000 0.3559 1.0000 0.3559
## [4,] 0.0000 0.0000 0.3559 1.0000
##   Standard Deviations: 1 1 1 1 

#########################################################################
#
#   SAS and R use different conventions to calculate AIC and BIC.  We
#   have used ML here, in which case AIC is the same but BIC differs as
#   noted above.  If we'd used REML, both AIC and BIC values are calculated
#   differently by SAS and R using different conventions regarding the
#   number of observations and number of parameters, so are not comparable, 
#   but can be compared within a single implementation (but not across
#   SAS and R).
#
#   We use the anova() function to print out these values for the above
#   models.  Both AIC and BIC support the common compound symmetric correlation
#   model with different variance for each gender (constant across time),
#   model (c)'.  This is as close as the models we can fit here with 
#   the generic gls() function can get to the model with different
#   compound correlation and variance for each gender we fit in SAS,
#   which was preferred among those fit in that program.
#
#   It seems that the data are suggesting that the correlation pattern
#   is close to compound symmetric, with heterogeneity of some sort 
#   between genders.
#
#########################################################################

## > anova(dental.un,dental.cs,dental.cs2,dental.ar1,dental.ar12,dental.1d)
##             Model df      AIC      BIC    logLik   Test  L.Ratio p-value
## dental.un       1 14 447.4770 485.0269 -209.7385                        
## dental.cs       2  6 440.6391 456.7318 -214.3195 1 vs 2  9.16201  0.3288
## dental.cs2      3  7 430.6521 449.4270 -208.3261 2 vs 3 11.98695  0.0005
## dental.ar1      4  6 452.6810 468.7738 -220.3405 3 vs 4 24.02890  <.0001
## dental.ar12     5  7 443.0267 461.8016 -214.5134 4 vs 5 11.65430  0.0006
## dental.1d       6  6 469.4166 485.5094 -228.7083 5 vs 6 28.38990  <.0001

#  We now fit full and reduced models with this covariance structure
#  and compare via likelihood ratio test

full <- gls(distance ~ -1 + gender + gender:age,data=thedat,correlation=corCompSymm(form = ~ 1 | child),
                 weights = varIdent(form = ~ 1 | gender),method="ML")
reduced <- gls(distance ~ -1 + gender + age,data=thedat,correlation=corCompSymm(form = ~ 1 | child),
               weights = varIdent(form = ~ 1 | gender),method="ML")

## > anova(full,reduced)
##         Model df      AIC      BIC    logLik   Test  L.Ratio p-value
## full        1  7 430.6521 449.4270 -208.3261                        
## reduced     2  6 436.7324 452.8252 -212.3662 1 vs 2 8.080332  0.0045

#  And of course we can test the difference in slopes directly using a
#  Wald-type test

full.alt <- gls(distance ~ gender + age + gender:age,data=thedat,correlation=corCompSymm(form = ~ 1 | child), weights = varIdent(form = ~ 1 | gender),method="ML")

#  Get the corrected standard errors

u <- full.alt
robust.full.alt <- robust.cov(u,m)
se.full.alt.corrected <- robust.full.alt$se.model

#  This is the table of coefficients and standard errors from summary(full.alt)

## Coefficients:
##                 Value Std.Error   t-value p-value
## (Intercept) 17.372727 0.8175862 21.248802  0.0000
## gender1     -1.032102 1.4222595 -0.725678  0.4697
## age          0.479545 0.0612544  7.828755  0.0000
## gender1:age  0.304830 0.1065571  2.860716  0.0051

#  Here are the corrected SEs - the conclusion about the difference in
#  slopes would not change

## > se.full.alt.corrected
## (Intercept)     gender1         age gender1:age 
##  0.80230287  1.39567285  0.06010933  0.10456520 



