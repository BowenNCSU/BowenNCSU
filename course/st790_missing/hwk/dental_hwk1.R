## ST 790, Homework 1, Problems 4 and 5

## Dental data with induced missingness and LOCF

library(gee)

#  Actual, full data

dent1 <- read.table("dental.dat")
colnames(dent1) <- c("num","child","age","distance","gender")

gee.full <- gee(distance ~ age + gender+ age*gender, id=child, family=gaussian,
corstr="exchangeable",data=dent1)

## > summary(gee.full)

##  GEE:  GENERALIZED LINEAR MODELS FOR DEPENDENT DATA
##  gee S-function, version 4.13 modified 98/01/27 (1998) 

## Model:
##  Link:                      Identity 
##  Variance to Mean Relation: Gaussian 
##  Correlation Structure:     Exchangeable 

## Call:
## gee(formula = distance ~ age + gender + age * gender, id = child, 
##     data = dent1, family = gaussian, corstr = "exchangeable")

## Summary of Residuals:
##        Min         1Q     Median         3Q        Max 
## -5.6156250 -1.3218750 -0.1681818  1.3299006  5.2468750 


## Coefficients:
##               Estimate Naive S.E.    Naive z Robust S.E.
## (Intercept) 17.3727273 1.19173082 14.5777276   0.7252063
## age          0.4795455 0.09502472  5.0465336   0.0631326
## gender      -1.0321023 1.54810375 -0.6666881   1.3777851
## age:gender   0.3048295 0.12344073  2.4694405   0.1168673
##               Robust z
## (Intercept) 23.9555672
## age          7.5958453
## gender      -0.7491025
## age:gender   2.6083390

## Estimated Scale Parameter:  5.093818
## Number of Iterations:  1

## Working Correlation
##           [,1]      [,2]      [,3]      [,4]
## [1,] 1.0000000 0.6100109 0.6100109 0.6100109
## [2,] 0.6100109 1.0000000 0.6100109 0.6100109
## [3,] 0.6100109 0.6100109 1.0000000 0.6100109
## [4,] 0.6100109 0.6100109 0.6100109 1.0000000

#  Fake data with induced missingness

dent2 <- read.table("dental_dropout_R.dat")
colnames(dent2) <- c("num","child","age","distance","gender")

gee.avail <- gee(distance ~ age + gender+ age*gender, id=child, family=gaussian,
corstr="exchangeable",na.action=na.omit,data=dent2)

## > summary(gee.avail)

##  GEE:  GENERALIZED LINEAR MODELS FOR DEPENDENT DATA
##  gee S-function, version 4.13 modified 98/01/27 (1998) 

## Model:
##  Link:                      Identity 
##  Variance to Mean Relation: Gaussian 
##  Correlation Structure:     Exchangeable 

## Call:
## gee(formula = distance ~ age + gender + age * gender, id = child, 
##     data = dent2, na.action = na.omit, family = gaussian, corstr = "exchangeable")

## Summary of Residuals:
##        Min         1Q     Median         3Q        Max 
## -5.6572339 -1.2993507  0.0217077  1.3651028  5.2428124 


## Coefficients:
##               Estimate Naive S.E.    Naive z Robust S.E.
## (Intercept) 16.6280091  1.3188736 12.6077349  0.95569671
## age          0.5671342  0.1136168  4.9916411  0.08451198
## gender      -0.1706828  1.7140404 -0.0995792  1.36942624
## age:gender   0.2078543  0.1476876  1.4073914  0.11846312
##               Robust z
## (Intercept) 17.3988347
## age          6.7106955
## gender      -0.1246382
## age:gender   1.7545906

## Estimated Scale Parameter:  5.191848
## Number of Iterations:  3

## Working Correlation
##           [,1]      [,2]      [,3]      [,4]
## [1,] 1.0000000 0.6080346 0.6080346 0.6080346
## [2,] 0.6080346 1.0000000 0.6080346 0.6080346
## [3,] 0.6080346 0.6080346 1.0000000 0.6080346
## [4,] 0.6080346 0.6080346 0.6080346 1.0000000

#  Reconfigure full data and get means and test  

full <- dent1[,2:5]
full.wide <- reshape(full,v.names="distance",idvar="child",
                       timevar="age",direction="wide")

## > full.wide[1:10,]
##    child gender distance.8 distance.10 distance.12 distance.14
## 1      1      0       21.0        20.0        21.5        23.0
## 5      2      0       21.0        21.5        24.0        25.5
## 9      3      0       20.5        24.0        24.5        26.0
## 13     4      0       23.5        24.5        25.0        26.5
## 17     5      0       21.5        23.0        22.5        23.5
## 21     6      0       20.0        21.0        21.0        22.5
## 25     7      0       21.5        22.5        23.0        25.0
## 29     8      0       23.0        23.0        23.5        24.0
## 33     9      0       20.0        21.0        22.0        21.5
## 37    10      0       16.5        19.0        19.0        19.5

#  t.test at age 14

## > t.test(full.wide$distance.14[full.wide$gender==0],full.wide$distance.14[full.wide$gender==1],var.equal=TRUE,mu=0,conf.level=0.95)

## 	Two Sample t-test

## data:  full.wide$distance.14[full.wide$gender == 0] and full.wide$distance.14[full.wide$gender == 1]
## t = -3.8623, df = 25, p-value = 0.000705
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -5.179034 -1.576648
## sample estimates:
## mean of x mean of y 
##  24.09091  27.46875 

mar <- dent2[,2:5]
mar.wide <-  reshape(mar,v.names="distance",idvar="child",
                       timevar="age",direction="wide")

## > mar.wide[1:10,]
##    child gender distance.8 distance.10 distance.12 distance.14
## 1      1      0       21.0        20.0        21.5          NA
## 5      2      0       21.0        21.5        24.0        25.5
## 9      3      0       20.5        24.0        24.5        26.0
## 13     4      0       23.5        24.5        25.0        26.5
## 17     5      0       21.5        23.0        22.5        23.5
## 21     6      0       20.0        21.0          NA          NA
## 25     7      0       21.5        22.5        23.0        25.0
## 29     8      0       23.0          NA          NA          NA
## 33     9      0       20.0        21.0        22.0        21.5
## 37    10      0       16.5        19.0          NA          NA

#  create LOCF data by brute force

locf <- mar.wide
locf[is.na(mar.wide[,4]),4]<-locf[is.na(mar.wide[,4]),3]
locf[is.na(mar.wide[,5]),5]<-locf[is.na(mar.wide[,5]),4]
locf[is.na(mar.wide[,6]),6]<-locf[is.na(mar.wide[,6]),5]

## >  locf[1:10,]
##      child gender distance.8 distance.10 distance.12 distance.14
## 1      1      0       21.0        20.0        21.5        21.5
## 5      2      0       21.0        21.5        24.0        25.5
## 9      3      0       20.5        24.0        24.5        26.0
## 13     4      0       23.5        24.5        25.0        26.5
## 17     5      0       21.5        23.0        22.5        23.5
## 21     6      0       20.0        21.0        21.0        21.0
## 25     7      0       21.5        22.5        23.0        25.0
## 29     8      0       23.0        23.0        23.0        23.0
## 33     9      0       20.0        21.0        22.0        21.5
## 37    10      0       16.5        19.0        19.0        19.0

#  t-test with LOCF data

## > t.test(locf$distance.14[locf$gender==0],locf$distance.14[locf$gender==1],var.equal=TRUE,mu=0,conf.level=0.95)

## 	Two Sample t-test

## data:  locf$distance.14[locf$gender == 0] and locf$distance.14[locf$gender == 1]
## t = -2.0152, df = 25, p-value = 0.05475
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -4.62416162  0.05029798
## sample estimates:
## mean of x mean of y 
##  23.68182  25.96875 



