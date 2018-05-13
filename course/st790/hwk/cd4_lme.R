#########################################################################
#
#   ACTG 193 CD4 data
#
#   1 = zidovudine alternating monthly with 400mg didanosine
#   2 = zidovudine plus 2.25mg of zalcitabine
#   3 = zidovudine plus 400mg of didanosine
#   4 = zidovudine plus 400mg of didanosine plus 400mg of nevirapine
#
#   We collapse 1,2,3 into a single "dual therapy" group and compare to
#   the "triple therapy" group 4
#
#   Homework 3 linear mixed model with splines
#
#########################################################################

library(nlme)
library(Matrix)
library(magic)

##  Read in the data -- they are in the "long" form of one record per
##  observation

thedat <- read.table("cd4.dat")
colnames(thedat) <- c("id","tx","age","gender","week","logcd4")

m <- length(unique(thedat$id))

weekplus <- (thedat$week-16)*(thedat$week-16 > 0)
trt <- as.factor(thedat$tx)
dual <- as.factor(as.numeric(thedat$tx==4))

cntl <- lmeControl(maxIter=200,msMaxIter=200,tolerance=1e-6,msTol=1e-7)

## We cannot posit models where the D matrix is diffrent for different
## treatment groups, as lme() does not support this. But we can have
## different variances for each group.  

##  Common D, common diagonal within-individual covariance

fit.a <- lme(logcd4 ~ week:dual + weekplus:dual, data=thedat,
              random = ~ week + weekplus | as.factor(id), method="ML", control=cntl)
beta.a <- fixed.effects(fit.a)
b.a <- random.effects(fit.a)
sebeta.a.incorrect <- summary(fit.a)$tTable[,"Std.Error"]  #  This gives incorrect values
sebeta.a <- sqrt(diag(fit.a$varFix))  #  This gives correct values
D.a <- getVarCov(fit.a,type="random.effects")
V.a <- getVarCov(fit.a,type="marginal",individual=1)
R.a <- getVarCov(fit.a,type="conditional",individual=1)
Dcorr.a <- cov2cor(D.a)
Vcorr.a <- cov2cor(simplify2array(getVarCov(fit.a,type="marginal",individual=1)[[1]]))

## The results agree with SAS proc mixed

## > cbind(beta.a,sebeta.a)
##                      beta.a    sebeta.a
## (Intercept)     2.941462232 0.025610961
## week:dual0     -0.007344209 0.001985066
## week:dual1      0.019507436 0.003335242
## dual0:weekplus -0.012036856 0.003171017
## dual1:weekplus -0.039774172 0.005337134
## > D.a
## Random effects variance covariance matrix
##             (Intercept)        week   weekplus
## (Intercept)   0.5850900  0.00726380 -0.0123560
## week          0.0072638  0.00091741 -0.0009105
## weekplus     -0.0123560 -0.00091050  0.0012253
##   Standard Deviations: 0.76491 0.030289 0.035005 
## > Dcorr.a
## Random effects variance covariance matrix
##             (Intercept)     week weekplus
## (Intercept)     1.00000  0.31352 -0.46145
## week            0.31352  1.00000 -0.85876
## weekplus       -0.46145 -0.85876  1.00000
##   Standard Deviations: 1 1 1 
## > fit.a$sigma^2
## [1] 0.3061655

## Common D, different diagonal within-individual covariance

fit.b <- lme(logcd4 ~ week:dual + weekplus:dual, data=thedat,
             random = ~ week + weekplus | as.factor(id),
             weights = varIdent(form = ~ 1 | dual),method="ML", control=cntl)
beta.b <- fixed.effects(fit.b)
b.b <- random.effects(fit.b)
sebeta.b.incorrect <- summary(fit.b)$tTable[,"Std.Error"]  #  This gives incorrect values
sebeta.b <- sqrt(diag(fit.b$varFix))  #  This gives correct values
D.b <- getVarCov(fit.b,type="random.effects")
V.b <- getVarCov(fit.b,type="marginal",individual=1)
R.b.0 <- getVarCov(fit.b,type="conditional",individual=1)
R.b.1 <- getVarCov(fit.b,type="conditional",individual=2)
Dcorr.b <- cov2cor(D.b)
Vcorr.b <- cov2cor(simplify2array(getVarCov(fit.b,type="marginal",individual=1)[[1]]))

## The results agree with SAS proc mixed

## > cbind(beta.b,sebeta.b)
##                      beta.b    sebeta.b
## (Intercept)     2.941282890 0.025620884
## week:dual0     -0.007326445 0.001995116
## week:dual1      0.019508349 0.003304294
## dual0:weekplus -0.012052277 0.003188510
## dual1:weekplus -0.039785382 0.005272845
## > D.b
## Random effects variance covariance matrix
##             (Intercept)        week    weekplus
## (Intercept)    0.586130  0.00715700 -0.01216500
## week           0.007157  0.00093006 -0.00092895
## weekplus      -0.012165 -0.00092895  0.00125200
##   Standard Deviations: 0.76559 0.030497 0.035384 
## > Dcorr.b
## Random effects variance covariance matrix
##             (Intercept)     week weekplus
## (Intercept)     1.00000  0.30653 -0.44905
## week            0.30653  1.00000 -0.86085
## weekplus       -0.44905 -0.86085  1.00000
##   Standard Deviations: 1 1 1 
## > R.b.0
## as.factor(id) 1 
## Conditional variance covariance matrix
##         1       2       3       4       5       6
## 1 0.30916 0.00000 0.00000 0.00000 0.00000 0.00000
## 2 0.00000 0.30916 0.00000 0.00000 0.00000 0.00000
## 3 0.00000 0.00000 0.30916 0.00000 0.00000 0.00000
## 4 0.00000 0.00000 0.00000 0.30916 0.00000 0.00000
## 5 0.00000 0.00000 0.00000 0.00000 0.30916 0.00000
## 6 0.00000 0.00000 0.00000 0.00000 0.00000 0.30916
##   Standard Deviations: 0.55602 0.55602 0.55602 0.55602 0.55602 0.55602 
## > R.b.1
## as.factor(id) 2 
## Conditional variance covariance matrix
##         1       2       3       4       5       6
## 1 0.29567 0.00000 0.00000 0.00000 0.00000 0.00000
## 2 0.00000 0.29567 0.00000 0.00000 0.00000 0.00000
## 3 0.00000 0.00000 0.29567 0.00000 0.00000 0.00000
## 4 0.00000 0.00000 0.00000 0.29567 0.00000 0.00000
## 5 0.00000 0.00000 0.00000 0.00000 0.29567 0.00000
## 6 0.00000 0.00000 0.00000 0.00000 0.00000 0.29567
##   Standard Deviations: 0.54375 0.54375 0.54375 0.54375 0.54375 0.54375 

## Common D, same within-individual variance and exponential
## correlation with no measurement error

fit.c <- lme(logcd4 ~ week:dual + weekplus:dual, data=thedat,
             random = ~ week + weekplus | as.factor(id),
             correlation=corExp(form = ~ week | as.factor(id), nugget=FALSE),
             method="ML", control=cntl)
beta.c <- fixed.effects(fit.c)
b.c <- random.effects(fit.c)
sebeta.c.incorrect <- summary(fit.c)$tTable[,"Std.Error"]  #  This gives incorrect values
sebeta.c <- sqrt(diag(fit.c$varFix))  #  This gives correct values
D.c <- getVarCov(fit.c,type="random.effects")
V.c <- getVarCov(fit.c,type="marginal",individual=1)
R.c <- getVarCov(fit.c,type="conditional",individual=1)
Dcorr.c <- cov2cor(D.c)
Vcorr.c <- cov2cor(simplify2array(getVarCov(fit.c,type="marginal",individual=1)[[1]]))
Rcorr.c <- cov2cor(simplify2array(getVarCov(fit.c,type="conditional",individual=1)[[1]]))

## Agrees for the most part with SAS -- the within correlation is
## driven to zero, suggesting that this model is overkill

## > cbind(beta.b,sebeta.b)
##                      beta.b    sebeta.b
## (Intercept)     2.941282890 0.025620884
## week:dual0     -0.007326445 0.001995116
## week:dual1      0.019508349 0.003304294
## dual0:weekplus -0.012052277 0.003188510
## dual1:weekplus -0.039785382 0.005272845
## > D.c
## Random effects variance covariance matrix
##             (Intercept)        week    weekplus
## (Intercept)   0.5850700  0.00726440 -0.01235600
## week          0.0072644  0.00091744 -0.00091061
## weekplus     -0.0123560 -0.00091061  0.00122560
##   Standard Deviations: 0.7649 0.030289 0.035009 
## > Dcorr.c
## Random effects variance covariance matrix
##             (Intercept)     week weekplus
## (Intercept)     1.00000  0.31355 -0.46144
## week            0.31355  1.00000 -0.85875
## weekplus       -0.46144 -0.85875  1.00000
##   Standard Deviations: 1 1 1 
## > R.c
## as.factor(id) 1 
## Conditional variance covariance matrix
##             1           2           3           4
## 1  3.0616e-01  3.6276e-73 3.6450e-149 3.6625e-225
## 2  3.6276e-73  3.0616e-01  3.0763e-77 3.0911e-153
## 3 3.6450e-149  3.0763e-77  3.0616e-01  3.0763e-77
## 4 3.6625e-225 3.0911e-153  3.0763e-77  3.0616e-01
## 5 1.1644e-310 9.8276e-239 9.7807e-163  9.7340e-87
## 6  0.0000e+00 2.6465e-309 2.6339e-233 2.6213e-157
##             5           6
## 1 1.1644e-310  0.0000e+00
## 2 9.8276e-239 2.6465e-309
## 3 9.7807e-163 2.6339e-233
## 4  9.7340e-87 2.6213e-157
## 5  3.0616e-01  8.2449e-72
## 6  8.2449e-72  3.0616e-01
##   Standard Deviations: 0.55332 0.55332 0.55332 0.55332 0.55332 0.55332 
## > Rcorr.c
##               1             2             3             4
## 1  1.000000e+00  1.184858e-72 1.190539e-148 1.196248e-224
## 2  1.184858e-72  1.000000e+00  1.004795e-76 1.009613e-152
## 3 1.190539e-148  1.004795e-76  1.000000e+00  1.004795e-76
## 4 1.196248e-224 1.009613e-152  1.004795e-76  1.000000e+00
## 5 3.803282e-310 3.209906e-238 3.194587e-162  3.179342e-86
## 6  0.000000e+00 8.644133e-309 8.602881e-233 8.561826e-157
##               5             6
## 1 3.803282e-310  0.000000e+00
## 2 3.209906e-238 8.644133e-309
## 3 3.194587e-162 8.602881e-233
## 4  3.179342e-86 8.561826e-157
## 5  1.000000e+00  2.692955e-71
## 6  2.692955e-71  1.000000e+00

## Common D, same within-individual variance and exponential
## correlation with measurement error

fit.d <- lme(logcd4 ~ week:dual + weekplus:dual, data=thedat,
             random = ~ week + weekplus | as.factor(id),
             correlation=corExp(form = ~ week | as.factor(id), nugget=TRUE),
             method="ML", control=cntl)
beta.d <- fixed.effects(fit.d)
b.d <- random.effects(fit.d)
sebeta.d.incorrect <- summary(fit.d)$tTable[,"Std.Error"]  #  This gives incorrect values
sebeta.d <- sqrt(diag(fit.d$varFix))  #  This gives correct values
D.d <- getVarCov(fit.d,type="random.effects")
V.d <- getVarCov(fit.d,type="marginal",individual=1)
R.d <- getVarCov(fit.d,type="conditional",individual=1)
Dcorr.d <- cov2cor(D.d)
Vcorr.d <- cov2cor(simplify2array(getVarCov(fit.d,type="marginal",individual=1)[[1]]))
Rcorr.d <- cov2cor(simplify2array(getVarCov(fit.d,type="conditional",individual=1)[[1]]))

## These results also agree with SAS

## > cbind(beta.d,sebeta.d)
##                      beta.d    sebeta.d
## (Intercept)     2.941462211 0.025611073
## week:dual0     -0.007344222 0.001985067
## week:dual1      0.019507424 0.003335246
## dual0:weekplus -0.012036862 0.003171022
## dual1:weekplus -0.039774166 0.005337144
## > D.d
## Random effects variance covariance matrix
##             (Intercept)        week    weekplus
## (Intercept)   0.5851000  0.00726380 -0.01235600
## week          0.0072638  0.00091744 -0.00091054
## weekplus     -0.0123560 -0.00091054  0.00122540
##   Standard Deviations: 0.76492 0.030289 0.035006 
## > Dcorr.d
## Random effects variance covariance matrix
##             (Intercept)     week weekplus
## (Intercept)     1.00000  0.31351 -0.46143
## week            0.31351  1.00000 -0.85876
## weekplus       -0.46143 -0.85876  1.00000
## > R.d
## as.factor(id) 1 
## Conditional variance covariance matrix
##             1           2           3           4
## 1  3.0616e-01  1.4038e-95 3.2610e-195 7.5751e-295
## 2  1.4038e-95  3.0616e-01 6.4478e-101 1.4978e-200
## 3 3.2610e-195 6.4478e-101  3.0616e-01 6.4478e-101
## 4 7.5751e-295 1.4978e-200 6.4478e-101  3.0616e-01
## 5  0.0000e+00 1.2225e-312 5.2627e-213 2.2655e-113
## 6  0.0000e+00  0.0000e+00 1.5981e-305 6.8796e-206
##             5           6
## 1  0.0000e+00  0.0000e+00
## 2 1.2225e-312  0.0000e+00
## 3 5.2627e-213 1.5981e-305
## 4 2.2655e-113 6.8796e-206
## 5  3.0616e-01  8.4288e-94
## 6  8.4288e-94  3.0616e-01
##   Standard Deviations: 0.55332 0.55332 0.55332 0.55332 0.55332 0.55332 
## > Rcorr.d
##               1             2             3             4
## 1  1.000000e+00  4.585075e-95 1.065105e-194 2.474220e-294
## 2  4.585075e-95  1.000000e+00 2.105995e-100 4.892190e-200
## 3 1.065105e-194 2.105995e-100  1.000000e+00 2.105995e-100
## 4 2.474220e-294 4.892190e-200 2.105995e-100  1.000000e+00
## 5  0.000000e+00 3.993052e-312 1.718933e-212 7.399683e-113
## 6  0.000000e+00  0.000000e+00 5.219876e-305 2.247058e-205
##               5             6
## 1  0.000000e+00  0.000000e+00
## 2 3.993052e-312  0.000000e+00
## 3 1.718933e-212 5.219876e-305
## 4 7.399683e-113 2.247058e-205
## 5  1.000000e+00  2.753040e-93
## 6  2.753040e-93  1.000000e+00

range.d <- exp(coef(fit.d$modelStruct$corStruct))[1]
nugget.d <- exp(coef(fit.d$modelStruct$corStruct))[2]

#  correlation parameter

## > (1-nugget.d)*exp(-1/range.d)
##       nugget 
## 3.151607e-13


## We could do still fancier models with different within individual
## variances and within-individual correlation, but it does not seem
## necessary - there does not seem to be any evidence of
## within-individual correlation.  So compare the models above:

## > anova(fit.a,fit.b,fit.c,fit.d)
##       Model df      AIC      BIC    logLik   Test
## fit.a     1 12 11912.12 11990.41 -5944.061       
## fit.b     2 13 11913.58 11998.40 -5943.790 1 vs 2
## fit.c     3 13 11914.12 11998.94 -5944.060       
## fit.d     4 14 11916.12 12007.46 -5944.061 3 vs 4
##         L.Ratio p-value
## fit.a                  
## fit.b 0.5418002  0.4617
## fit.c                  
## fit.d 0.0005038  0.9821

## AIC prefers larger models but prefers the common D/common
## within-indiv variance model a; BIC also prefers this model.  In
## SAS, model a is preferred by BIC, but AIC prefers the model with
## both different D and within-indiv variance, which we can't fit in
## R.  If using R, either model a would be a reasonable choice --
## we'll use model a for simplicity in further analyses.  Given the
## missing data, we need to recognize that standard errors and tests
## based on the usual asymptotic theory might be unreliable;
## likelihood ratio tests of nested models might be a better approach.
## Proceeding with Wald tests anyway just to get a sense...

#  Fit three models with successively more dependence on age and gender

#  age and gender in intercept

fit.a.1 <- lme(logcd4 ~ age + gender + week:dual + weekplus:dual, data=thedat,
              random = ~ week + weekplus | as.factor(id), method="ML", control=cntl)
beta.a.1 <- fixed.effects(fit.a.1)
b.a.1 <- random.effects(fit.a.1)
sebeta.a.1 <- sqrt(diag(fit.a.1$varFix))  #  This gives correct values
D.a.1 <- getVarCov(fit.a.1,type="random.effects")
Dcorr.a.1 <- cov2cor(D.a.1)

df <- nrow(thedat)-length(beta.a.1)
pvalue <- round(2*pt(-abs(beta.a.1/sebeta.a.1),df),4)

## > cbind(beta.a.1,sebeta.a.1,beta.a.1/sebeta.a.1,pvalue)
##                    beta.a.1  sebeta.a.1           pvalue
## (Intercept)     2.645538522 0.127850052 20.692510 0.0000
## age             0.009997378 0.003012036  3.319143 0.0009
## gender         -0.092720017 0.075281729 -1.231640 0.2181
## week:dual0     -0.007295563 0.001984526 -3.676224 0.0002
## week:dual1      0.019494089 0.003332681  5.849372 0.0000
## dual0:weekplus -0.012060066 0.003170775 -3.803507 0.0001
## dual1:weekplus -0.039783604 0.005335973 -7.455736 0.0000

#  age and gender in intercept and phase 1 slope

fit.a.2 <- lme(logcd4 ~ age + gender + week:dual + week:dual:age + week:dual:gender+
              weekplus:dual, data=thedat,
              random = ~ week + weekplus | as.factor(id), method="ML", control=cntl)
beta.a.2 <- fixed.effects(fit.a.2)
b.a.2 <- random.effects(fit.a.2)
sebeta.a.2 <- sqrt(diag(fit.a.2$varFix))  #  This gives correct values
D.a.2 <- getVarCov(fit.a.2,type="random.effects")
Dcorr.a.2 <- cov2cor(D.a.2)

df <- nrow(thedat)-length(beta.a.2)
pvalue <- round(2*pt(-abs(beta.a.2/sebeta.a.2),df),4)

## > cbind(beta.a.2,sebeta.a.2,beta.a.2/sebeta.a.2,pvalue)
##                        beta.a.2   sebeta.a.2           
## (Intercept)        2.638906e+00 0.1316912848 20.0385808
## age                9.371085e-03 0.0031058511  3.0172359
## gender            -5.826258e-02 0.0775726123 -0.7510715
## week:dual0        -1.028815e-02 0.0055744209 -1.8456004
## week:dual1         3.240412e-02 0.0090442849  3.5828284
## dual0:weekplus    -1.200098e-02 0.0031701823 -3.7855810
## dual1:weekplus    -3.936371e-02 0.0053371988 -7.3753506
## age:week:dual0     1.341437e-04 0.0001247604  1.0752102
## age:week:dual1    -3.988605e-05 0.0001942728 -0.2053095
## gender:week:dual0 -2.325016e-03 0.0031919292 -0.7284046
## gender:week:dual1 -1.302178e-02 0.0052761317 -2.4680546
##                   pvalue
## (Intercept)       0.0000
## age               0.0026
## gender            0.4526
## week:dual0        0.0650
## week:dual1        0.0003
## dual0:weekplus    0.0002
## dual1:weekplus    0.0000
## age:week:dual0    0.2823
## age:week:dual1    0.8373
## gender:week:dual0 0.4664
## gender:week:dual1 0.0136

fit.a.3 <- lme(logcd4 ~ age + gender + week:dual + week:dual:age + week:dual:gender+
              weekplus:dual + weekplus:dual:age + weekplus:dual:gender, data=thedat,
              random = ~ week + weekplus | as.factor(id), method="ML", control=cntl)
beta.a.3 <- fixed.effects(fit.a.3)
b.a.3 <- random.effects(fit.a.3)
sebeta.a.3 <- sqrt(diag(fit.a.3$varFix))  #  This gives correct values
D.a.3 <- getVarCov(fit.a.3,type="random.effects")
Dcorr.a.3 <- cov2cor(D.a.3)

df <- nrow(thedat)-length(beta.a.3)
pvalue <- round(2*pt(-abs(beta.a.3/sebeta.a.3),df),4)

## > cbind(beta.a.3,sebeta.a.3,beta.a.3/sebeta.a.3,pvalue)
##                            beta.a.3   sebeta.a.3
## (Intercept)            2.6619361097 0.1326210563
## age                    0.0089570280 0.0031284272
## gender                -0.0666413473 0.0781166283
## week:dual0            -0.0235797554 0.0104003149
## week:dual1             0.0270521839 0.0174973120
## dual0:weekplus         0.0131197862 0.0169731530
## dual1:weekplus        -0.0289450162 0.0277992162
## age:week:dual0         0.0004627801 0.0002462426
## age:week:dual1        -0.0001874799 0.0003982203
## gender:week:dual0     -0.0012813369 0.0062604343
## gender:week:dual1     -0.0006357206 0.0103553537
## age:dual0:weekplus    -0.0006206680 0.0003993344
## age:dual1:weekplus     0.0002826907 0.0006260142
## gender:dual0:weekplus -0.0019836895 0.0101913559
## gender:dual1:weekplus -0.0237917338 0.0171031467
##                                   pvalue
## (Intercept)           20.07174566 0.0000
## age                    2.86310896 0.0042
## gender                -0.85310066 0.3936
## week:dual0            -2.26721553 0.0234
## week:dual1             1.54607655 0.1221
## dual0:weekplus         0.77297283 0.4396
## dual1:weekplus        -1.04121699 0.2978
## age:week:dual0         1.87936615 0.0603
## age:week:dual1        -0.47079438 0.6378
## gender:week:dual0     -0.20467221 0.8378
## gender:week:dual1     -0.06139052 0.9511
## age:dual0:weekplus    -1.55425634 0.1202
## age:dual1:weekplus     0.45157240 0.6516
## gender:dual0:weekplus -0.19464432 0.8457
## gender:dual1:weekplus -1.39107348 0.1643

## There seems to be little evidence that gender has much to do with
## anything, and age seems to only have a strong effect in intercept.
## So fit a reduced model with no gender just to see if age emerges
## anywhere else

fit.a.4 <- lme(logcd4 ~ age  + week:dual + week:dual:age +
              weekplus:dual + weekplus:dual:age, data=thedat,
              random = ~ week + weekplus | as.factor(id), method="ML", control=cntl)
beta.a.4 <- fixed.effects(fit.a.4)
b.a.4 <- random.effects(fit.a.4)
sebeta.a.4 <- sqrt(diag(fit.a.4$varFix))  #  This gives correct values
D.a.4 <- getVarCov(fit.a.4,type="random.effects")
Dcorr.a.4 <- cov2cor(D.a.4)

df <- nrow(thedat)-length(beta.a.4)
pvalue <- round(2*pt(-abs(beta.a.4/sebeta.a.4),df),4)

## > cbind(beta.a.4,sebeta.a.4,beta.a.4/sebeta.a.4,pvalue)
##                        beta.a.4   sebeta.a.4           
## (Intercept)         2.6145472355 0.1201614105 21.7586264
## age                 0.0086635235 0.0031120548  2.7838596
## week:dual0         -0.0244311459 0.0094037536 -2.5980206
## week:dual1          0.0261761277 0.0154591074  1.6932496
## dual0:weekplus      0.0113870019 0.0152819310  0.7451285
## dual1:weekplus     -0.0470831263 0.0244141928 -1.9285146
## age:week:dual0      0.0004554178 0.0002444986  1.8626601
## age:week:dual1     -0.0001758824 0.0003974477 -0.4425296
## age:dual0:weekplus -0.0006237726 0.0003975095 -1.5692017
## age:dual1:weekplus  0.0001919947 0.0006249804  0.3072011
##                    pvalue
## (Intercept)        0.0000
## age                0.0054
## week:dual0         0.0094
## week:dual1         0.0905
## dual0:weekplus     0.4562
## dual1:weekplus     0.0538
## age:week:dual0     0.0626
## age:week:dual1     0.6581
## age:dual0:weekplus 0.1167
## age:dual1:weekplus 0.7587

## It seems that there's really not a lot going on here, except in the
## intercept.  Fit a final model.

fit.a.5 <- lme(logcd4 ~ age  + week:dual + weekplus:dual, data=thedat,
              random = ~ week + weekplus | as.factor(id), method="ML", control=cntl)
beta.a.5 <- fixed.effects(fit.a.5)
b.a.5 <- random.effects(fit.a.5)
sebeta.a.5 <- sqrt(diag(fit.a.5$varFix))  #  This gives correct values
D.a.5 <- getVarCov(fit.a.5,type="random.effects")
Dcorr.a.5 <- cov2cor(D.a.5)

df <- nrow(thedat)-length(beta.a.5)
pvalue <- round(2*pt(-abs(beta.a.5/sebeta.a.5),df),4)

## > cbind(beta.a.5,sebeta.a.5,beta.a.5/sebeta.a.5,pvalue)
##                    beta.a.5  sebeta.a.5           pvalue
## (Intercept)     2.578995830 0.115931769 22.245808 0.0000
## age             0.009606495 0.002997209  3.205147 0.0014
## week:dual0     -0.007305451 0.001984548 -3.681167 0.0002
## week:dual1      0.019499079 0.003332855  5.850564 0.0000
## dual0:weekplus -0.012071805 0.003170830 -3.807143 0.0001
## dual1:weekplus -0.039809533 0.005335812 -7.460819 0.0000

## Having established that there is evidence of an age effect in
## intercept, we can look at whether or not there is evidence that the
## phase 1 and phase 2 slopes are different between dual and triple therapy.
## We do this by fitting a reduced model with the same phase 1 and
## phase 2 slopes and by fitting the above model in an alternative parameterization

fit.a.6 <- lme(logcd4 ~ age  + week + weekplus, data=thedat,
              random = ~ week + weekplus | as.factor(id), method="ML", control=cntl)
beta.a.6 <- fixed.effects(fit.a.6)
b.a.6 <- random.effects(fit.a.6)
sebeta.a.6 <- sqrt(diag(fit.a.6$varFix))  #  This gives correct values
D.a.6 <- getVarCov(fit.a.6,type="random.effects")
Dcorr.a.6 <- cov2cor(D.a.6)

df <- nrow(thedat)-length(beta.a.6)
pvalue <- round(2*pt(-abs(beta.a.6/sebeta.a.6),df),4)

## > cbind(beta.a.6,sebeta.a.6,beta.a.6/sebeta.a.6,pvalue)
##                  beta.a.6  sebeta.a.6            pvalue
## (Intercept)  2.5757727839 0.115978258 22.2091005 0.0000
## age          0.0096891869 0.002998439  3.2314100 0.0012
## week        -0.0004294285 0.001757240 -0.2443767 0.8069
## weekplus    -0.0191660816 0.002760984 -6.9417570 0.0000

## > anova(fit.a.5,fit.a.6)
##         Model df      AIC      BIC    logLik   Test
## fit.a.5     1 13 11903.91 11988.72 -5938.952       
## fit.a.6     2 11 11957.73 12029.50 -5967.864 1 vs 2
##          L.Ratio p-value
## fit.a.5                 
## fit.a.6 57.82251  <.0001

## Evidence is strong that phase 1 and/or phase 2 slopes are different.

fit.a.5.alt <- lme(logcd4 ~ age  + week + week:dual + weekplus + weekplus:dual, data=thedat,
              random = ~ week + weekplus | as.factor(id), method="ML", control=cntl)
beta.a.5.alt <- fixed.effects(fit.a.5.alt)
b.a.5.alt <- random.effects(fit.a.5.alt)
sebeta.a.5.alt <- sqrt(diag(fit.a.5.alt$varFix))  #  This gives correct values
D.a.5.alt <- getVarCov(fit.a.5.alt,type="random.effects")
Dcorr.a.5.alt <- cov2cor(D.a.5.alt)

df <- nrow(thedat)-length(beta.a.5.alt)
pvalue <- round(2*pt(-abs(beta.a.5.alt/sebeta.a.5.alt),df),4)

## > cbind(beta.a.5.alt,sebeta.a.5.alt,beta.a.5.alt/sebeta.a.5.alt,pvalue)
##                beta.a.5.alt sebeta.a.5.alt          
## (Intercept)     2.578995833    0.115931769 22.245808
## age             0.009606495    0.002997209  3.205147
## week           -0.007305451    0.001984548 -3.681167
## weekplus       -0.012071805    0.003170830 -3.807143
## week:dual1      0.026804530    0.003839967  6.980406
## dual1:weekplus -0.027737728    0.006189390 -4.481496
##                pvalue
## (Intercept)    0.0000
## age            0.0014
## week           0.0002
## weekplus       0.0001
## week:dual1     0.0000
## dual1:weekplus 0.0000

## Going back to the final model in the original parameterization, we
## can get estimates of the mean rate of change during the second part
## of the study and for mean log CD4 at the end of the study

## Large sample covariance matrix based on expected information

Sig.a.5 <- fit.a.5$varFix
df <- nrow(thedat)-length(beta.a.5)

est.L <- function(L,beta,Sig,alpha,df){
    est <- L%*%beta
    se <- sqrt(L%*%Sig%*%t(L))
    crit.val <- abs(qt(alpha/2,df))
    CI <- c(est-crit.val*se,est+crit.val*se)
    result <- matrix(c(est,se,CI),1,4,byrow=TRUE)
    colnames(result) <- c("Estimate","Std Error","Lower CI", "Upper CI")
    return(result)
}

##  Estimate mean slope in the second part of the study for each group

## slope post-16 weeks, dual therapy

L <- matrix(c(0,0,1,0,1,0),nrow=1,byrow=TRUE)

## > est.L(L,beta.a.5,Sig.a.5,0.05,df)
##         Estimate   Std Error    Lower CI    Upper CI
## [1,] -0.01937726 0.001755837 -0.02281946 -0.01593505

## slope post-16 weeks, triple therapy

L <- matrix(c(0,0,0,1,0,1),nrow=1,byrow=TRUE)

## > est.L(L,beta.a.5,Sig.a.5,0.05,df)
##         Estimate   Std Error    Lower CI    Upper CI
## [1,] -0.02031045 0.002936487 -0.02606725 -0.01455366

##  Compare slopes in second part of study

L <- matrix(c(0,0,1,-1,1,-1),nrow=1,byrow=TRUE)

## > est.L(L,beta.a.5,Sig.a.5,0.05,df)
##          Estimate   Std Error     Lower CI    Upper CI
## [1,] 0.0009331975 0.003420377 -0.005772231 0.007638626

## Estimate mean log CD4 at 40 weeks; this depends on age because
## baseline depends on age.  We evaluate this at the average age for
## the study as an example

avgage <- round(mean(thedat$age[thedat$week<4]),2)
## > avgage
## [1] 37.77

## Dual therapy

L <- matrix(c(1,avgage,40,0,40-16,0),nrow=1,byrow=TRUE)

## >  est.L(L,beta.a.5,Sig.a.5,0.05,df)
##      Estimate  Std Error Lower CI Upper CI
## [1,] 2.359892 0.04338566 2.274837 2.444947

## Triple therapy

L <- matrix(c(1,avgage,0,40,0,40-16),nrow=1,byrow=TRUE)

## > est.L(L,beta.a.5,Sig.a.5,0.05,df)
##      Estimate  Std Error Lower CI Upper CI
## [1,] 2.766368 0.06867822 2.631728 2.901007

## Difference Triple - Dual

L <- matrix(c(0,0,-40,40,-(40-16),(40-16)),nrow=1,byrow=TRUE)

## > est.L(L,beta.a.5,Sig.a.5,0.05,df)
##       Estimate  Std Error  Lower CI  Upper CI
## [1,] 0.4064757 0.07718521 0.2551591 0.5577924
