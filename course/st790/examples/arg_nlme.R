############################################################
#
#  Chapter 9, Argartroban pharmacokinetic data             
#                                                          
#  Fit nonlinear mixed effects model with within-individual         
#  power-of-mean variance function using nlme
#
############################################################

library(nlme)

############################################################
#                                                          
#  Read in the data 
#                                                          
############################################################

thedat <- read.table("argconc.dat")
colnames(thedat) <- c("obsno","id","dose","time","conc")

############################################################
#                                                          
#  Define the within-individual conditional mean model.
#  Here, we specify the function and its matrix of partial
#  derivatives.  If the derivatives are not provided, 
#  nlme will use numeric derivatives.  Analytic derivatives
#  generally speed up the calculation.
#
############################################################	
	
pkmodel <- function(t,b1,b2,dose){
	tinf <- 240
	cl <- exp(b1)
	v <- exp(b2)
        t1 <- t<=tinf
        t2 <- tinf*(1-t1)+t1*t
        f1 <- (dose/cl)*(1-exp(-cl*t2/v))*exp(-cl*(1-t1)*(t-tinf)/v)

#  compute analytical dervivatives -- gradient matrix X

        t3 <- (1-t1)*(t-tinf)
	temp1 <- (dose/cl)*exp(-cl*t3/v)
        temp2 <- (dose/(v^2))*exp(-cl*t3/v)
        
   gradmat <- array(0,c(length(t),2),list(NULL,c("b1","b2")))
   gradmat[,"b1"] <- temp1*(-(1-exp(-cl*t2/v))*(1/cl + t3/v)+(t2/v)*exp(-cl*t2/v))*cl 
   gradmat[,"b2"] <- temp2*((1-exp(-cl*t2/v))*t3-exp(-cl*t2/v)*t2)*v
 
   attr(f1,"gradient") <- gradmat
   f1
}

############################################################
#     
#   nlme() uses the method in the course notes based on the
#   linearizing the model about current empirical Bayes
#   estimates of the random effects.  Estimation of the 
#   covariance parameters is by using the normal likelihood
#   (method="ML", the default) or using a REML version 
#   (method="REML")

#   verbose=TRUE prints out history of each iteration
#
#   start provides starting values for the fixed effects beta
#
#   The second stage model is specified by fixed = and 
#   random = ; here, there are no among-individual covariates
#
#   The weights statement specifies the within-individual 
#   variance model -- here it is taken to be the power-of-mean
#   model with starting value for the power delta = 0.5
#
############################################################

arg.mlfit1 <- nlme(conc ~ pkmodel(time,b1,b2,dose),
    fixed = list(b1 ~ 1,b2 ~1),
    random = list(b1 ~ 1,b2 ~ 1),
    groups = ~ id,
    weights = varPower(0.5),
    start = list(fixed = c(-6.0,-2.0)),
    data = thedat,                  
    method="ML",verbose=TRUE)

## > summary(arg.mlfit1)
## Nonlinear mixed-effects model fit by maximum likelihood
##   Model: conc ~ pkmodel(time, b1, b2, dose) 
##  Data: thedat 
##        AIC      BIC    logLik
##   5738.425 5767.568 -2862.212

## Random effects:
##  Formula: list(b1 ~ 1, b2 ~ 1)
##  Level: id
##  Structure: General positive-definite, Log-Cholesky parametrization
##          StdDev      Corr 
## b1        0.37168388 b1   
## b2        0.06753305 0.268
## Residual 20.41960002      

## Variance function:
##  Structure: Power of variance covariate
##  Formula: ~fitted(.) 
##  Parameter estimates:
##     power 
## 0.2432885 
## Fixed effects: list(b1 ~ 1, b2 ~ 1) 
##        Value  Std.Error  DF   t-value p-value
## b1 -5.432545 0.06230331 437 -87.19513       0
## b2 -1.917953 0.02513017 437 -76.32073       0
##  Correlation: 
##    b1   
## b2 0.156

## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -2.48366041 -0.42424579  0.03283586  0.57935560  9.38411228 

## Number of Observations: 475
## Number of Groups: 37

## scale parameter

## > arg.mlfit1$sigma        
## 20.4196
## > arg.mlfit1$sigma^2     
## 416.9601

## random effects estimates (first 8 individuals)
## >  arg.mlfit1$coef$random
## $id
##              b1           b2
## 1   0.264321672  0.018991541
## 2   0.105589914  0.003912753
## 3   0.356888627  0.012955221
## 4   0.251467974  0.020897031
## 5   0.283323599  0.021841625
## 6   0.598379352  0.031227550
## 7  -0.382613189 -0.033197546
## 8   0.553777073  0.033120548

#  It takes some contortions to extract the matrix D -- from the output
#  above, we can compute it ourselves as

D <- diag(c(0.37168388,0.06753305),2,2)%*%matrix(c(1,0.268,0.268,1),2,2)%*%diag(c(0.37168388,0.06753305),2,2)
## > D
##             [,1]        [,2]
## [1,] 0.138148907 0.006727054
## [2,] 0.006727054 0.004560713

#  Or we can use the VarCorr() function to extract this stuff -- still
#  a pain; why there is not a straightforward way to extract D from
#  an nlme object is a mystery -- what you get with open source
 #  software

## > VarCorr(arg.mlfit1)
## id = pdLogChol(list(b1 ~ 1,b2 ~ 1)) 
##          Variance     StdDev      Corr 
## b1       1.381489e-01  0.37168388 b1   
## b2       4.560714e-03  0.06753305 0.268
## Residual 4.169601e+02 20.41960002   

D.T <- diag(VarCorr(arg.mlfit1)[1:2,2])
## > D.T
##           [,1]       [,2]
## [1,] 0.3716839 0.00000000
## [2,] 0.0000000 0.06753305

D.corr <- as.numeric(VarCorr(arg.mlfit1)[2,3])
## > D.corr
## [1] 0.268
D.Gam <- diag(2)
D.Gam[1,2]<-D.corr; D.Gam[2,1]<-D.corr
## > D.Gam
##       [,1]  [,2]
## [1,] 1.000 0.268
## [2,] 0.268 1.000

D.arg.mlfit1 <- D <- D.T%*%D.Gam%*%D.T
## > D.arg.mlfit1
##             [,1]        [,2]
## [1,] 0.138148907 0.006727054
## [2,] 0.006727054 0.004560713


############################################################
#     
#   It is also possible to specify a within-indiviudal correlation 
#   model using the correlation= option.  This is not really 
#   that pleasing here, as the main source of within-individual
#   variation is thought to be measurement error, but in principle
#   this could be done.  The times here are not equally spaced,
#   so using something like the exponential correlation model would
#   be necessary.  This would entail adding
#
#   correlation=corExp(form ~ time | id, nugget=FALSE)
#
#   to the nlme() call.  This does not end well - it generates an
#   error, which probably reflects that it is not really needed.
#
############################################################


