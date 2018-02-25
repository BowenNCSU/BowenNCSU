##############################################################################
#
#   Using gls() to fit the dental data with artificial missingness
#   and unstructured covariance matrix 
#
##############################################################################

library(nlme)

#  read in the data

dent1 <- read.table("dental_dropout_R.dat")
colnames(dent1) <- c("num","child","age","distance","sex")

#  create time variable for each individual as required for the CorSymm
#  correlation structure

time <- rep(seq(1,4,1),max(dent1$child))

#  create class variable for gender

gender <- factor(dent1$sex)

dental.un <- gls(distance ~ -1 + gender + age:gender,correlation=corSymm(form = ~ time | factor(child)),
                weights = varIdent(form = ~ 1|age),method="ML",na.action=na.omit,data=dent1)

## > summary(dental.un)
## Generalized least squares fit by maximum likelihood
##   Model: distance ~ -1 + gender + age:gender 
##   Data: dent1 
##       AIC      BIC   logLik
##   375.028 409.7107 -173.514

## Correlation Structure: General
##  Formula: ~time | factor(child) 
##  Parameter estimate(s):
##  Correlation: 
##   1     2     3    
## 2 0.538            
## 3 0.615 0.487      
## 4 0.729 0.771 0.673
## Variance function:
##  Structure: Different standard deviations per stratum
##  Formula: ~1 | age 
##  Parameter estimates:
##         8        10        12        14 
## 1.0000000 0.9205712 1.0853978 0.9548465 

## Coefficients:
##                 Value Std.Error   t-value p-value
## gender0     16.896144 0.9872384 17.114553       0
## gender1     16.174496 0.8221987 19.672249       0
## gender0:age  0.536666 0.0845572  6.346778       0
## gender1:age  0.792489 0.0704744 11.245069       0

##  Correlation: 
##             gendr0 gendr1 gndr0:
## gender1      0.000              
## gender0:age -0.807  0.000       
## gender1:age  0.000 -0.807  0.000

## Standardized residuals:
##         Min          Q1         Med          Q3         Max 
## -2.44389127 -0.55645921  0.05390317  0.67555783  2.20952657 

## Residual standard error: 2.256406 
## Degrees of freedom: 88 total; 84 residual

#  This is the SAS code using proc mixed to produce the same analysis

#  data dent2; infile 'dental_dropout_sas.dat';
#    input obsno child age distance gender;
#  run;

#  proc mixed method=ml data=dent2;
#    class  child;
#    model distance = gender gender*age / noint solution chisq ;
#    repeated / type=un subject=child r=2 rcorr=2;
#  run;






