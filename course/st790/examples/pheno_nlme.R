############################################################
#  
#  Chapter 9, Pharmacokinetic study of phenobarbital
#
#  Fit nonlinear mixed effects model with within-individual         
#  power-of-mean variance function using nlme
#
############################################################

library(nlme)

############################################################
#  
#  The data are part of the nlme distribution -- the data
#  set is Phenobarb.  Here is the first two infant
#
## > Phenobarb[1:30,]
## Grouped Data: conc ~ time | Subject
##    Subject  Wt Apgar ApgarInd  time dose conc
## 1        1 1.4     7     >= 5   0.0 25.0   NA
## 2        1 1.4     7     >= 5   2.0   NA 17.3
## 3        1 1.4     7     >= 5  12.5  3.5   NA
## 4        1 1.4     7     >= 5  24.5  3.5   NA
## 5        1 1.4     7     >= 5  37.0  3.5   NA
## 6        1 1.4     7     >= 5  48.0  3.5   NA
## 7        1 1.4     7     >= 5  60.5  3.5   NA
## 8        1 1.4     7     >= 5  72.5  3.5   NA
## 9        1 1.4     7     >= 5  85.3  3.5   NA
## 10       1 1.4     7     >= 5  96.5  3.5   NA
## 11       1 1.4     7     >= 5 108.5  3.5   NA
## 12       1 1.4     7     >= 5 112.5   NA 31.0
#
############################################################

############################################################
#
#  The one compartment model with intravenous administration
#  in (9.30) of the notes is available in nlme() as 
#  phenoModel - see the calls to nlme() below.  The model 
#  is parameterized in terms of log Cl_i and log V_i.
#  The diligent student might want to try his/her hand at 
#  programming this model.
#
#  The model is programmed to expect NA values in the dose and 
#  conc variables; the na.action=na.pass option tells nlme() to 
#  include the NAs so that the function can identify, dosing
#  and sampling times, and the naPattern option tells nlme
#  that values of conc that are not NAs are to be retained as the
#  response variable and NAs are to be ignored.
#  
#  We fit three models: (i) with no among-individual covariates
#  in either parameter, (ii) with the among-individual covariate
#  birthweight in both parameters, and (iii) with both among-individual
#  covariates birthweight and indicator of whether or not apgar
#  score >= 5 in both parameters.
#  
#  After each fit, we plot the empirical Bayes estimates of random
#  effects for log Cl_i and log V_i against covariates to look at 
#  possible systematic relationships.  We use the plotting capability
#  in nlme() to make these plots.
#
############################################################

# (i) No covariates

pheno.mlfit.nocov <- nlme(conc ~ phenoModel(Subject,time,dose,lCl,lV),
          data=Phenobarb,               
          fixed = list(lCl ~ 1, lV ~ 1),
          random = list(lCl ~ 1, lV ~ 1),
          groups = ~ Subject,
          weights = varPower(fixed=1), 
          start = list(fixed = c(-5,0.3)),
          na.action=na.pass, naPattern = ~ !is.na(conc),                
          method="ML",verbose=TRUE)
## > summary(pheno.mlfit.nocov)
## Nonlinear mixed-effects model fit by maximum likelihood
##   Model: conc ~ phenoModel(Subject, time, dose, lCl, lV) 
##  Data: Phenobarb 
##        AIC      BIC    logLik
##   986.6177 1004.878 -487.3089

## Random effects:
##  Formula: list(lCl ~ 1, lV ~ 1)
##  Level: Subject
##  Structure: General positive-definite, Log-Cholesky parametrization
##          StdDev    Corr 
## lCl      0.5088706 lCl  
## lV       0.4006982 0.951
## Residual 0.1129108      

## Variance function:
##  Structure: Power of variance covariate
##  Formula: ~fitted(.) 
##  Parameter estimates:
## power 
##     1 
## Fixed effects: list(lCl ~ 1, lV ~ 1) 
##         Value  Std.Error DF   t-value p-value
## lCl -4.988225 0.07610856 95 -65.54091       0
## lV   0.332197 0.05433713 95   6.11363       0
##  Correlation: 
##    lCl  
## lV 0.751

## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -2.88501446 -0.31261344  0.08948398  0.36911612  2.42826168 

## Number of Observations: 155
## Number of Groups: 59

#  Extract the empirical Bayes estimates of random effects
#  for each parameter and plot against covariates

pdf("pheno.nocov.pdf",width=8)
nocov.raneff <- ranef(pheno.mlfit.nocov, augFrame=TRUE)
plot(nocov.raneff,form = lCl~ Wt + ApgarInd )
plot(nocov.raneff,form = lV~ Wt + ApgarInd )
graphics.off()

#  (ii)  Add birthweight as a linear covariate in each parameter
#   and plot empirical Bayes estimates to see if the apparent
#  systematic relationship persists.

pheno.mlfit.wtonly <- nlme(conc ~ phenoModel(Subject,time,dose,lCl,lV),
          data=Phenobarb,               
          fixed = list(lCl ~ Wt, lV ~ Wt),
          random = list(lCl ~ 1, lV ~ 1),
          groups = ~ Subject,
          weights = varPower(fixed=1), 
          start = list(fixed = c(-5,0,0.3,0)),
          na.action=na.pass, naPattern = ~ !is.na(conc),                
          method="ML",verbose=TRUE)                        
## > summary(pheno.mlfit.wtonly)
## Nonlinear mixed-effects model fit by maximum likelihood
##   Model: conc ~ phenoModel(Subject, time, dose, lCl, lV) 
##  Data: Phenobarb 
##        AIC      BIC    logLik
##   883.2707 907.6181 -433.6353

## Random effects:
##  Formula: list(lCl ~ 1, lV ~ 1)
##  Level: Subject
##  Structure: General positive-definite, Log-Cholesky parametrization
##                 StdDev    Corr  
## lCl.(Intercept) 0.2292397 lC.(I)
## lV.(Intercept)  0.1409428 0.763 
## Residual        0.1140465       

## Variance function:
##  Structure: Power of variance covariate
##  Formula: ~fitted(.) 
##  Parameter estimates:
## power 
##     1 
## Fixed effects: list(lCl ~ Wt, lV ~ Wt) 
##                     Value  Std.Error DF   t-value p-value
## lCl.(Intercept) -5.959398 0.12403306 93 -48.04685       0
## lCl.Wt           0.639351 0.07691382 93   8.31256       0
## lV.(Intercept)  -0.491629 0.05710098 93  -8.60981       0
## lV.Wt            0.538492 0.03466386 93  15.53468       0
##  Correlation: 
##                lC.(I) lCl.Wt lV.(I)
## lCl.Wt         -0.926              
## lV.(Intercept)  0.090 -0.053       
## lV.Wt          -0.052  0.043 -0.913

## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -2.78605290 -0.37208927  0.04861949  0.42820401  2.19427910 

## Number of Observations: 155
## Number of Groups: 59

pdf("pheno.cov.pdf",width=8)
wtonly.raneff <- ranef(pheno.mlfit.wtonly, augFrame=TRUE)
plot(wtonly.raneff,form = lCl.(Intercept) ~ Wt + ApgarInd )
plot(wtonly.raneff,form = lV.(Intercept) ~ Wt + ApgarInd )
graphics.off()

#  (iii) Add Apgar score category to the previous model (ii) 
#  in each parameter.  It does not appear that this is necessary

pheno.mlfit.wtapgar <- nlme(conc ~ phenoModel(Subject,time,dose,lCl,lV),
          data=Phenobarb,               
          fixed = list(lCl ~ Wt + ApgarInd, lV ~ Wt + ApgarInd),
          random = list(lCl ~ 1, lV ~ 1),
          groups = ~ Subject,
          weights = varPower(fixed=1), 
          start = list(fixed = c(-5,0,0,0.3,0,0)),
          na.action=na.pass, naPattern = ~ !is.na(conc),                
          method="ML",verbose=TRUE)                        
## > summary(pheno.mlfit.wtapgar)
## Nonlinear mixed-effects model fit by maximum likelihood
##   Model: conc ~ phenoModel(Subject, time, dose, lCl, lV) 
##  Data: Phenobarb 
##        AIC      BIC    logLik
##   887.2203 917.6546 -433.6102

## Random effects:
##  Formula: list(lCl ~ 1, lV ~ 1)
##  Level: Subject
##  Structure: General positive-definite, Log-Cholesky parametrization
##                 StdDev    Corr  
## lCl.(Intercept) 0.2290444 lC.(I)
## lV.(Intercept)  0.1406943 0.763 
## Residual        0.1141014       

## Variance function:
##  Structure: Power of variance covariate
##  Formula: ~fitted(.) 
##  Parameter estimates:
## power 
##     1 
## Fixed effects: list(lCl ~ Wt + ApgarInd, lV ~ Wt + ApgarInd) 
##                      Value  Std.Error DF    t-value p-value
## lCl.(Intercept)  -5.932316 0.21800171 91 -27.212245  0.0000
## lCl.Wt            0.636990 0.08098199 91   7.865828  0.0000
## lCl.ApgarInd>= 5 -0.026936 0.16165674 91  -0.166622  0.8680
## lV.(Intercept)   -0.483991 0.07800085 91  -6.204950  0.0000
## lV.Wt             0.539216 0.03503128 91  15.392416  0.0000
## lV.ApgarInd>= 5  -0.010669 0.06462644 91  -0.165081  0.8692
##  Correlation: 
##                  lC.(I) lCl.Wt lC.AI5 lV.(I) lV.Wt 
## lCl.Wt           -0.740                            
## lCl.ApgarInd>= 5 -0.817  0.277                     
## lV.(Intercept)   -0.070  0.047  0.074              
## lV.Wt             0.043  0.018 -0.092 -0.651       
## lV.ApgarInd>= 5   0.064 -0.092  0.007 -0.672 -0.037

## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -2.78831468 -0.37414083  0.05443072  0.42344909  2.18863006 

## Number of Observations: 155
## Number of Groups: 59

#  Likelihood ratio tests of these nested models suggest that while including
#  birthweigh is important, Apgar score category does not add further value.
#  There is strong evidience to suggest that both log Cl_i and log V_i
#  are systematically associated with birthweight, but no evidence of
#  such an association with Apgar score.

## >  anova(pheno.mlfit.nocov,pheno.mlfit.wtonly, pheno.mlfit.wtapgar)
##                     Model df      AIC       BIC    logLik   Test   L.Ratio
## pheno.mlfit.nocov       1  6 986.6177 1004.8783 -487.3089                 
## pheno.mlfit.wtonly      2  8 883.2707  907.6181 -433.6353 1 vs 2 107.34701
## pheno.mlfit.wtapgar     3 10 887.2203  917.6546 -433.6102 2 vs 3   0.05036
##                     p-value
## pheno.mlfit.nocov          
## pheno.mlfit.wtonly   <.0001
## pheno.mlfit.wtapgar  0.9751

wtapgar.raneff <- ranef(pheno.mlfit.wtapgar, augFrame=TRUE)
plot(wtapgar.raneff,form = lCl.(Intercept) ~ Wt + ApgarInd )
plot(wtapgar.raneff,form = lV.(Intercept) ~ Wt + ApgarInd )

