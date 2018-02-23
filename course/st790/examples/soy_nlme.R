############################################################
#  
#  Chapter 9, Soybean study
#
#  Fit nonlinear mixed effects model with within-individual         
#  power-of-mean variance function using nlme
#
############################################################

library(nlme)

############################################################
#  
#  The data are part of the nlme distribution -- the data
#  set is Soybean.  Here are the first 2 plots.
#
#  NOTE:  The data that come with nlme are not exactly the 
#  same as the data used in the original analyses of this
#  study; the response has been rounded to 5 decimal places, 
#  whereas in the original data set provided by the investigators
#  the response was available to 6 decimal places.  Remarkably,
#  this is sufficient to yield substantial differences in the
#  estimates of some parameters, in particular covariance
#  parameters.
#
############################################################

## > Soybean[1:16,]
## Grouped Data: weight ~ Time | Plot
##      Plot Variety Year Time   weight
## 1  1988F1       F 1988   14  0.10600
## 2  1988F1       F 1988   21  0.26100
## 3  1988F1       F 1988   28  0.66600
## 4  1988F1       F 1988   35  2.11000
## 5  1988F1       F 1988   42  3.56000
## 6  1988F1       F 1988   49  6.23000
## 7  1988F1       F 1988   56  8.71000
## 8  1988F1       F 1988   63 13.35000
## 9  1988F1       F 1988   70 16.34170
## 10 1988F1       F 1988   77 17.75083
## 11 1988F2       F 1988   14  0.10400
## 12 1988F2       F 1988   21  0.26900
## 13 1988F2       F 1988   28  0.77800
## 14 1988F2       F 1988   35  2.12000
## 15 1988F2       F 1988   42  2.93000
## 16 1988F2       F 1988   49  5.29000
	
############################################################
#  
#  Define the mean response function.  Here, we program    
#  both the function and its derivatives with respect to  
#  each parameter.  If analytic derivatives are not       
#  provided by the user, nlme will use numeric derivatives
#
#  Use the reparameterization of the model in the course notes
#  so that the scale parameter is more likely to be normally
#  distributed.  We also scale by 100 to put all the parameters
#  on an equal magnitude basis.
#
############################################################	
	
logistf <- function(b1,b2,b3,day){
	denom <- 1+exp(-b3*(day-b2))
	denom2 <- denom*denom
	f1 <- b1/denom

	qfunc <- f1

#  analytical dervivatives

   qgrad <- array(0,c(length(day),3),list(NULL,c("b1","b2","b3")))
   qgrad[,"b1"] <- (1/denom)
   qgrad[,"b2"] <- -(b3*b1*(denom-1)/denom2)
   qgrad[,"b3"] <- ((day-b2)*(b1)*(denom-1)/denom2)
   attr(qfunc,"gradient") <- qgrad
	qfunc
}

#  Initial model with no covariates

soybean.mlfit.nocov <- nlme(weight ~ logistf(b1,b2,b3,Time),
    fixed = list(b1 ~ 1, b2 ~ 1, b3 ~ 1),
    random = list(b1 ~ 1,b2 ~ 1,b3 ~ 1),
    groups = ~ Plot,
    weights = varPower(0.5),                      
    data = Soybean,
    start = list(fixed = c(20,50,0.125)),
    method="ML",verbose=TRUE)

## > summary(soybean.mlfit0)
## Nonlinear mixed-effects model fit by maximum likelihood
##   Model: weight ~ logistf(b1, b2, b3, Time) 
##  Data: Soybean 
##        AIC    BIC    logLik
##   744.3487 788.58 -361.1744

## Random effects:
##  Formula: list(b1 ~ 1, b2 ~ 1, b3 ~ 1)
##  Level: Plot
##  Structure: General positive-definite, Log-Cholesky parametrization
##          StdDev       Corr         
## b1       3.517176e+00 b1     b2    
## b2       8.663002e-07  0.977       
## b3       5.117989e-03 -1.000 -0.977
## Residual 2.286027e-01              

## Variance function:
##  Structure: Power of variance covariate
##  Formula: ~fitted(.) 
##  Parameter estimates:
##     power 
## 0.9413516 
## Fixed effects: list(b1 ~ 1, b2 ~ 1, b3 ~ 1) 
##       Value Std.Error  DF   t-value p-value
## b1 17.17029 0.6281987 362  27.33258       0
## b2 51.99796 0.3998056 362 130.05812       0
## b3  0.13186 0.0016793 362  78.52096       0
##  Correlation: 
##    b1     b2    
## b2  0.458       
## b3 -0.600 -0.758

## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -2.4647498 -0.5928431 -0.1271307  0.5572828  4.0636373 

## Number of Observations: 412
## Number of Groups: 48 

#  Note that in the above fit the estimated random effects standard
#  deviations for b2 and b3 are very small, and the correlations suggest
#  that the edtimated D may be very close to being singular -- we return
#  to this below

#  Get random effects estimates and plot against year-genotype to look for
#  systematic associations

pdf("soybean.nocov.pdf",width=8)
plot(ranef(soybean.mlfit.nocov,augFrame=TRUE),form = ~ Year*Variety,layout=c(3,1))
dev.off()

#  Parameterize the model as specified with different means for each
#  combination of year and genotype for each of the 3
#  individual-specific parameters

soybean.mlfit.cov <- nlme(weight ~ logistf(b1,b2,b3,Time),
    fixed = list(b1 ~ -1 + Variety:Year,
    b2 ~ -1  + Variety:Year,
    b3 ~ -1 + Variety:Year),
    random = list(b1 ~ 1,b2 ~ 1,b3 ~ 1),
    groups = ~ Plot,
    weights = varPower(0.5),                      
    data = Soybean,
    start = list(fixed = c(20,20,20,20,20,20,
	50,50,50,50,50,50,
	0.125,0.125,0.125,0.125,0.125,0.125)),
    method="ML",verbose=TRUE)
## > summary(soybean.mlfit.cov)
## Nonlinear mixed-effects model fit by maximum likelihood
##   Model: weight ~ logistf(b1, b2, b3, Time) 
##  Data: Soybean 
##        AIC      BIC    logLik
##   635.5266 740.0732 -291.7633

## Random effects:
##  Formula: list(b1 ~ 1, b2 ~ 1, b3 ~ 1)
##  Level: Plot
##  Structure: General positive-definite, Log-Cholesky parametrization
##                StdDev       Corr         
## b1.(Intercept) 1.031965e+00 b1.(I) b2.(I)
## b2.(Intercept) 4.242654e-12  0.000       
## b3.(Intercept) 7.766293e-07  0.000 -0.001
## Residual       2.180105e-01              

## Variance function:
##  Structure: Power of variance covariate
##  Formula: ~fitted(.) 
##  Parameter estimates:
##     power 
## 0.9415936 
## Fixed effects: list(b1 ~ -1 + Variety:Year, b2 ~ -1 + Variety:Year, b3 ~ -1 +      Variety:Year) 
##                         Value Std.Error  DF  t-value p-value
## b1.VarietyF:Year1988 19.12377 1.1270778 347 16.96757       0
## b1.VarietyP:Year1988 21.35907 1.1673480 347 18.29709       0
## b1.VarietyF:Year1989 10.45746 0.6558918 347 15.94388       0
## b1.VarietyP:Year1989 18.02238 0.9744041 347 18.49580       0
## b1.VarietyF:Year1990 15.88834 0.9343355 347 17.00497       0
## b1.VarietyP:Year1990 17.35646 0.9925907 347 17.48602       0
## b2.VarietyF:Year1988 54.50107 1.0293330 347 52.94795       0
## b2.VarietyP:Year1988 53.81469 1.0142704 347 53.05753       0
## b2.VarietyF:Year1989 52.23501 0.9365423 347 55.77432       0
## b2.VarietyP:Year1989 51.64393 0.9280962 347 55.64502       0
## b2.VarietyF:Year1990 49.84900 0.9283547 347 53.69608       0
## b2.VarietyP:Year1990 48.58598 0.9455062 347 51.38621       0
## b3.VarietyF:Year1988  0.12467 0.0032134 347 38.79819       0
## b3.VarietyP:Year1988  0.12336 0.0032344 347 38.13911       0
## b3.VarietyF:Year1989  0.14076 0.0038362 347 36.69223       0
## b3.VarietyP:Year1989  0.13812 0.0037722 347 36.61472       0
## b3.VarietyF:Year1990  0.13637 0.0039169 347 34.81505       0
## b3.VarietyP:Year1990  0.13404 0.0040456 347 33.13162       0

## Number of Observations: 412
## Number of Groups: 48 

#  Within-individual scale parameter (approx CV)
## >  soybean.mlfit$sigma      
## 0.218008 

#  Plot random effects estimates from this model to see if systematic
#  associations are taken into account

pdf("soybean.cov.pdf",width=8)
plot(ranef(soybean.mlfit.cov,augFrame=TRUE),form = ~ Year*Variety,layout=c(3,1))
dev.off()

#  As in the no covariate model, the estimated covariance matrix D has 
#  estimated variances that are almost zero.  It seems prudent to
#  introduce some sort of approximation; e.g., taking either b2 or b3
#  or both to have no associated random effect.  Both have very small
#  estimated random effects and estimated standard deviations that are
#  very small relative to the size of the fixed effect. Pinheiro and Bates
#  (2000) end up treating both b2 ("soybean half-life") and b3
#  ("growth rate parameter") as having no associated random effect.
#  We do this here for illustration

soybean.reduced.cov <- nlme(weight ~ logistf(b1,b2,b3,Time),
    fixed = list(b1 ~ -1 + Variety:Year,
    b2 ~ -1  + Variety:Year,
    b3 ~ -1 + Variety:Year),
    random = list(b1 ~ 1),
    groups = ~ Plot,
    weights = varPower(0.5),                      
    data = Soybean,
    start = list(fixed = c(20,20,20,20,20,20,
	50,50,50,50,50,50,
	0.125,0.125,0.125,0.125,0.125,0.125)),
    method="ML",verbose=TRUE)
## > summary(soybean.reduced.cov)
## Nonlinear mixed-effects model fit by maximum likelihood
##   Model: weight ~ logistf(b1, b2, b3, Time) 
##  Data: Soybean 
##        AIC      BIC    logLik
##   625.5176 709.9591 -291.7588

## Random effects:
##  Formula: b1 ~ 1 | Plot
##         b1.(Intercept) Residual
## StdDev:       1.032192 0.218009

## Variance function:
##  Structure: Power of variance covariate
##  Formula: ~fitted(.) 
##  Parameter estimates:
##     power 
## 0.9416202 
## Fixed effects: list(b1 ~ -1 + Variety:Year, b2 ~ -1 + Variety:Year, b3 ~ -1 +      Variety:Year) 
##                         Value Std.Error  DF  t-value p-value
## b1.VarietyF:Year1988 19.12255 1.1270538 347 16.96684       0
## b1.VarietyP:Year1988 21.35791 1.1673238 347 18.29647       0
## b1.VarietyF:Year1989 10.45518 0.6558214 347 15.94212       0
## b1.VarietyP:Year1989 18.01947 0.9742605 347 18.49553       0
## b1.VarietyF:Year1990 15.88508 0.9342223 347 17.00353       0
## b1.VarietyP:Year1990 17.35800 0.9928298 347 17.48336       0
## b2.VarietyF:Year1988 54.50003 1.0293233 347 52.94743       0
## b2.VarietyP:Year1988 53.81353 1.0142493 347 53.05750       0
## b2.VarietyF:Year1989 52.23181 0.9364445 347 55.77673       0
## b2.VarietyP:Year1989 51.64005 0.9279427 347 55.65004       0
## b2.VarietyF:Year1990 49.84744 0.9283797 347 53.69294       0
## b2.VarietyP:Year1990 48.58899 0.9457007 347 51.37883       0
## b3.VarietyF:Year1988  0.12468 0.0032135 347 38.79816       0
## b3.VarietyP:Year1988  0.12336 0.0032345 347 38.13913       0
## b3.VarietyF:Year1989  0.14077 0.0038364 347 36.69248       0
## b3.VarietyP:Year1989  0.13814 0.0037726 347 36.61585       0
## b3.VarietyF:Year1990  0.13637 0.0039171 347 34.81349       0
## b3.VarietyP:Year1990  0.13402 0.0040453 347 33.13051       0

## Number of Observations: 412
## Number of Groups: 48 

#  Within-individual scale parameter (approx CV)
## > soybean.reduced.cov$sigma
## [1] 0.218009

#  Based on this fit, can test main effects and interactions of year
#  and genotype for asymptotic growth value (beta1) and growth rate
#  constant (beta3).  I do this by constructing the tests myself.

Sigbeta <- soybean.reduced.cov$varFix
beta <- soybean.reduced.cov$coef$fixed

#  Function to contstruct F test statistics

f.test <- function(L,beta,Sig,dendf){
    numdf <- nrow(L)
    F <- t(L%*%beta)%*%solve(L%*%Sig%*%t(L))%*%(L%*%beta)/numdf
    p.value <- pf(F,numdf,dendf,lower.tail=FALSE)
    return(list(F=F,p.value=p.value,F.pvalue=c(F,round(p.value,4))))
}

#  Test statistic for main effect of genotype in
#  asymptotic growth value

L <- matrix(c(c(-1/3,1/3,-1/3,1/3,-1/3,1/3),rep(0,12)),nrow=1)
## > f.test(L,beta,Sigbeta,347)$F.pvalue
## [1] 22.63391  0.00000

#  Test for main effect of year (weather) in asymptotic growth value

L <- matrix(c(1/2,1/2,-1/2,-1/2,0,0,rep(0,12),
            1/2,1/2,0,0,-1/2,-1/2,rep(0,12)),nrow=2,byrow=TRUE)
## > f.test(L,beta,Sigbeta,347)$F.pvalue
## [1] 18.86112  0.00000

#  Test for interaction of year and genotype in asymptotic growth value

L <- matrix(c(1/2,-1/2,-1/2,1/2,0,0,rep(0,12),
            1/2,-1/2,0,0,-1/2,1/2,rep(0,12)),nrow=2,byrow=TRUE)
## >  f.test(L,beta,Sigbeta,347)$F.pvalue
## [1] 7.132952 0.000900

#  Test statistic for main effect of genotype in
#  growth rate 

L <- matrix(c(rep(0,12),c(-1/3,1/3,-1/3,1/3,-1/3,1/3)),nrow=1)
## > f.test(L,beta,Sigbeta,347)$F.pvalue
## [1] 0.5084351 0.4763000

#  Test for main effect of year (weather) in growth rate

L <- matrix(c(rep(0,12),1/2,1/2,-1/2,-1/2,0,0,
            rep(0,12),1/2,1/2,0,0,-1/2,-1/2),nrow=2,byrow=TRUE)
## > f.test(L,beta,Sigbeta,347)$F.pvalue
## [1] 11.08105  0.00000

#  Test for interaction of year and genotype in growth rate

L <- matrix(c(rep(0,12),1/2,-1/2,-1/2,1/2,0,0,
            rep(0,12),1/2,-1/2,0,0,-1/2,1/2),nrow=2,byrow=TRUE)
## > f.test(L,beta,Sigbeta,347)$F.pvalue
## [1] 0.02087844 0.97930000



