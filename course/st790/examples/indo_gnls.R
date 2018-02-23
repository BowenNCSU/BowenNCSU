#########################################################################
#                                                                       #
#  CHAPTER 7 - Example of using the GLS iterative algorithm to fit      #
#  a nonlinear model with variance function                             #
#                                                                       #
#  Data on the pharmacokinetics of indomethacin on a single subject     #  
#                                                                       #
#  Using the gnls() function in the nlme package in R to fit a          #
#  nonlinear mean model with a variance function with a variance        #
#  function with unknown parameters to be estimated.                    #
#                                                                       #                                      
#  Here, we fit the biexponential model for the mean E(Y_j | x_j) =     #
#                                                                       #
#     f(x,b) = exp(b1) exp(-exp(b2) x) + exp(b3) exp(-exp(b4) x),       #
#                                                                       #
#  which is a common PK model for a two-compartment representation of   #
#  the body, parameterized to ensure positivity of all the estimates.   #
#  We fit the power of the mean variance function                       #
#                                                                       #
#      var(Y_j | x_j) = sigma^2 f(x,b)^{2 delta}                        #
#                                                                       #
#  where delta is estimated by solving the quadratic estimating         #
#  equation (7.21) in the notes.                                        #
#                                                                       #
#  The gnls() function carries out the iterative gls algorithm in       #
#  Section 7.3 in the notes.                                            #
#                                                                       #
#########################################################################

library(nlme)

#########################################################################
#                                                                       #
#  The data are blood concentrations from a single subject; time is     #
#  in hours                                                             #
#                                                                       #                                        
#########################################################################

indodat <- read.table("indo.dat")
colnames(indodat) <- c("time","conc")

#########################################################################
#                                                                       #
#  Define the bioexponential mean function f and the gradient matrix    #
#  of its partial derivatives with respect to beta (n x p) as an        #
#  attribute.  The gnls() function will know to use analytic            #
#  derivatives when it spots the presence of the attribute "gradient"   #
#  along with the function.  Alternatively, one could leave this off    #
#  and let the function calculate numerical derivatives.  One can also  #
#  use the deriv() function to have R do this numerically               #
#                                                                       #
#########################################################################

biexp <- function(time,b1,b2,b3,b4){
   eb1 <- exp(b1) 
   eb2 <- exp(b2)
   eb3 <- exp(b3)
   eb4 <- exp(b4)
   indof <- eb1*exp(-eb2*time)+eb3*exp(-eb4*time)

#  compute analytical dervivatives -- create the gradient matrix X(beta)
	
   indograd <- array(0,c(length(time),4),list(NULL,c("b1","b2","b3","b4")))
   indograd[,"b1"] <- eb1*exp(-eb2*time)
   indograd[,"b2"] <- -eb1*eb2*time*exp(-eb2*time)
   indograd[,"b3"] <- eb3*exp(-eb4*time)
   indograd[,"b4"] <- -eb3*eb4*time*exp(-eb4*time)
   attr(indof,"gradient") <- indograd
   indof
}

#########################################################################
#                                                                       #
#  Call the gnls() function to implement the fit.  A starting value     #
#  for beta must be provided for a nonlinear model such as this one     #
#                                                                       #                                        
#########################################################################

indo.gls <- gnls(conc ~ biexp(time,b1,b2,b3,b4),start = c(b1=0.69,b2=0.69,b3=-1.6,b4=-1.6),
                 weights = varPower(1,form = ~ fitted(.)),data=indodat)

#  power = delta 

## >   summary(indo.gls)
## Generalized nonlinear least squares fit
##   Model: conc ~ biexp(time, b1, b2, b3, b4) 
##   Data: indodat 
##         AIC       BIC   logLik
##   -30.74966 -28.36229 21.37483

## Variance function:
##  Structure: Power of variance covariate
##  Formula: ~fitted(.) 
##  Parameter estimates:
##     power 
## 0.8189284 

## Coefficients:
##         Value Std.Error   t-value p-value
## b1  1.2365449 0.1768386  6.992507   2e-04
## b2  0.9700257 0.1379033  7.034101   2e-04
## b3 -1.4290571 0.2321700 -6.155219   5e-04
## b4 -1.7426471 0.2632090 -6.620772   3e-04

##  Correlation: 
##    b1    b2    b3   
## b2 0.825            
## b3 0.272 0.621      
## b4 0.228 0.545 0.932

## Standardized residuals:
##        Min         Q1        Med         Q3        Max 
## -0.9964536 -0.5850906 -0.3832389  0.5557971  1.3143110 

## Residual standard error: 0.1313473 
## Degrees of freedom: 11 total; 7 residual
