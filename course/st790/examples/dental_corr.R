#########################################################################
#
#   CHAPTER 2, EXAMPLE 1, Dental Study 
# 
#   Explore overall and within-child patterns of correlation in the
#   dental study data using R
#
#########################################################################

#  You may need to install the reshape package to run this code 
#
#  > install.packages("reshape")
#
#  We use the reshape function below to reconfigure the dental data to 
#  facilitate computation of summary statistics and plotting

#  Read in the data -- they are in the form of one record per observation
#  The gender indicator = 0 for girls and = 1 for boys

thedat <- read.table("dental.dat")
thedat <- thedat[,2:5]      #  remove the first column
colnames(thedat) <- c("id","age","distance","gender")

#  Total number of individuals

m <- max(thedat$id)

#########################################################################
# 
#  POPULATION-AVERAGED PERSPECTIVE
#
#########################################################################

#  Calculate girl- and boy-specific overall sample covariance and correlation
#  matrices and inspect to get a sense of overall pattern of variation and
#  correlation

#  Reconfigure the data set from one record per observation to one record
#  per child using the reshape function (so from "long" to "wide" format 
#  in R-speak) 

thedat.wide <- reshape(thedat,v.names="distance",idvar="id",
                       timevar="age",direction="wide")

#  Extract boys and girls separately

girlsdat <- thedat.wide[thedat.wide$gender==0,]
boysdat <- thedat.wide[thedat.wide$gender==1,]

#  Numbers of girls and boys

m.girls <- nrow(girlsdat)
m.boys <- nrow(boysdat)

#  Sample mean vectors -- we could instead fit a straight line to
#  the data on all m = 27 children and calculate the means as the predicted
#  values at 8, 10, 12, and 14 years

girls.mean <- apply(girlsdat[,3:6],2,mean)
boys.mean <- apply(boysdat[,3:6],2,mean)

#> round(girls.mean,3)
# distance.8 distance.10 distance.12 distance.14 
#     21.182      22.227      23.091      24.091 
#> round(boys.mean,3)
# distance.8 distance.10 distance.12 distance.14 
#     22.875      23.812      25.719      27.469 

#  Calculate the sample covariance matrix and correlation matrix for
#  each gender

girls.cov <- cov(girlsdat[,3:6])
boys.cov <- cov(boysdat[,3:6])

girls.corr <- cor(girlsdat[,3:6])  #  or use cov2cor(girls.cov)
boys.corr <- cor(boysdat[,3:6])

#> round(girls.cov,3)
#            distance.8 distance.10 distance.12 distance.14
#distance.8       4.514       3.355       4.332       4.357
#distance.10      3.355       3.618       4.027       4.077
#distance.12      4.332       4.027       5.591       5.466
#distance.14      4.357       4.077       5.466       5.941

#> round(boys.cov,3)
#            distance.8 distance.10 distance.12 distance.14
#distance.8       6.017       2.292       3.629       1.612
#distance.10      2.292       4.562       2.194       2.810
#distance.12      3.629       2.194       7.032       3.241
#distance.14      1.612       2.810       3.241       4.349

#> round(girls.corr,3)
#            distance.8 distance.10 distance.12 distance.14
#distance.8       1.000       0.830       0.862       0.841
#distance.10      0.830       1.000       0.895       0.879
#distance.12      0.862       0.895       1.000       0.948
#distance.14      0.841       0.879       0.948       1.000

#> round(boys.corr,3)
#            distance.8 distance.10 distance.12 distance.14
#distance.8       1.000       0.437       0.558       0.315
#distance.10      0.437       1.000       0.387       0.631
#distance.12      0.558       0.387       1.000       0.586
#distance.14      0.315       0.631       0.586       1.000

#  Center and scale to create scatterplot matrices using the scale and
#  pairs functions

girls.sd <- sqrt(diag(girls.cov))
boys.sd <- sqrt(diag(boys.cov))

girls.centscale <- scale(girlsdat[,3:6],center=girls.mean,scale=girls.sd)
boys.centscale <- scale(boysdat[,3:6],center=boys.mean,scale=boys.sd)

#  Now create the scatterplots

agelabs <- c("Age 8", "Age 10", "Age 12", "Age 14")

pdf("girlsscatter.pdf",width=8)
pairs(girls.centscale,label=agelabs,oma=c(5,5,5,5),main="Girls")
graphics.off()

pdf("boysscatter.pdf",width=8)
pairs(boys.centscale,label=agelabs,oma=c(5,5,5,5),main="Boys")
graphics.off()

#  Assuming a common pattern of overall covariance/correlation for 
#  boys and girls (but different overall means), calculate the pooled
#  covariance matrix and associated correlation matrix

cov.pooled <- ( (m.girls-1)*girls.cov + (m.boys-1)*boys.cov )/(m-2)
corr.pooled <- cov2cor(cov.pooled)

#> round(cov.pooled,3)
#            distance.8 distance.10 distance.12 distance.14
#distance.8       5.415       2.717       3.910       2.710
#distance.10      2.717       4.185       2.927       3.317
#distance.12      3.910       2.927       6.456       4.131
#distance.14      2.710       3.317       4.131       4.986

#> round(corr.pooled,3)
#            distance.8 distance.10 distance.12 distance.14
#distance.8       1.000       0.571       0.661       0.522
#distance.10      0.571       1.000       0.563       0.726
#distance.12      0.661       0.563       1.000       0.728
#distance.14      0.522       0.726       0.728       1.000

#  Autocorrelations at lags 1, 2, and 3 under stationarity for girls

ac.lag1 <- cor(c(girls.centscale[,1],girls.centscale[,2],girls.centscale[,3]),
           c(girls.centscale[,2],girls.centscale[,3],girls.centscale[,4]))
ac.lag2 <- cor(c(girls.centscale[,1],girls.centscale[,2]),
           c(girls.centscale[,3],girls.centscale[,4]))
ac.lag3 <- cor(girls.centscale[,1],girls.centscale[,4])

#> round(c(ac.lag1,ac.lag2,ac.lag3),3)
#[1] 0.891 0.871 0.841

#  Create lag plots for girls

pdf("girlslag.pdf",width=8)
par(mfrow=c(2,2))

plot(c(girls.centscale[,1],girls.centscale[,2],girls.centscale[,3]),
     c(girls.centscale[,2],girls.centscale[,3],girls.centscale[,4]),
     xlab="at age",ylab="at age + 2 years",xlim=c(-2,2),ylim=c(-2,2))
title("Lag 1 (2 years)")

plot(c(girls.centscale[,1],girls.centscale[,2]),
     c(girls.centscale[,3],girls.centscale[,4]),xlab="at age",
     ylab="at age + 4 years",xlim=c(-2,2),ylim=c(-2,2))
title("Lag 2 (4 years)")

plot(girls.centscale[,1],girls.centscale[,4],xlab="at age",
     ylab="at age + 6 years",xlim=c(-2,2),ylim=c(-2,2))
title("Lag 3 (6 years)")
graphics.off()

#  Autocorrelations at lags 1, 2, and 3 under stationarity for boys
#  (could be plotted as above)

ac.lag1.boys <- cor(c(boys.centscale[,1],boys.centscale[,2],boys.centscale[,3]),
           c(boys.centscale[,2],boys.centscale[,3],boys.centscale[,4]))
ac.lag2.boys <- cor(c(boys.centscale[,1],boys.centscale[,2]),
           c(boys.centscale[,3],boys.centscale[,4]))
ac.lag3.boys <- cor(boys.centscale[,1],boys.centscale[,4])

#> round(c(ac.lag1.boys,ac.lag2.boys,ac.lag3.boys),3)
#[1] 0.470 0.594 0.315

#########################################################################
# 
#  SUBJECT-SPECIFIC PERSPECTIVE
#
#########################################################################

#  Fit child-specific straight lines, and calculate child-specific
#  residuals.  Then calculate autocorrelations on the lagged residuals
#  and plot

#  Illustrate for boys -- same for girls (not shown)

#  Put data for each boy in a list and fit simple linear regression to each,
#  and extract the residuals and predicted values for each boy

boys.regdat <- thedat[thedat$gender==1,1:3]

boys.list <- as.list(1:m.boys)
for (i in 1:m.boys){
    boys.list[[i]] <- boys.regdat[boys.regdat$id==m.girls+i,]
}

boys.fits <- lapply(boys.list,FUN=function(u){lm(distance ~ age,data=u)})
boys.residuals <- sapply(boys.fits,residuals)
boys.pred <- sapply(boys.fits,fitted.values)

#  Assume that the within-child response variance (aggregate variance from 
#  realization and measurement error) is constant and the same for
#  all boys.  Estimate it by pooling the residuals across all boys
#  and then scale the residuals by its square root (the estimated within-child
#  response standard deviation).  The degrees of freedom in the denominator
#  are computed as the total number of observations minus the total
#  number of parameters estimated.  

boys.var <- sum(boys.residuals^2)/(m.boys*4-m.boys*2)
boys.stdresid <- boys.residuals/sqrt(boys.var)

#  Make a standard residual plot of the pooled residuals (standardized
#  residuals vs predicted values from all boys on the same plot)

pdf("boysresidplot.pdf",width=8)
matplot(boys.pred,boys.stdresid,xlab="predicted values",
        ylab="standardized residual",ylim=c(-3,3),pch=rep(1,m.boys))
abline(h=0)
title("Pooled, Standardized Least Squares Residuals, Boys")
graphics.off()

#  Calculate the estimated autocorrelation function

boys.lagresid <- t(boys.stdresid)

ac.resid1 <- cor(c(boys.lagresid[,1],boys.lagresid[,2],boys.lagresid[,3]),
     c(boys.lagresid[,2],boys.lagresid[,3],boys.lagresid[,4]))
ac.resid2 <- cor(c(boys.lagresid[,1],boys.lagresid[,2]),
                 c(boys.lagresid[,3],boys.lagresid[,4]))
ac.resid3 <- cor(boys.lagresid[,1],boys.lagresid[,4])

#> round(c(ac.resid1,ac.resid2,ac.resid3),3)
#[1] -0.685  0.144  0.290

#  Plot the pooled lagged residuals

pdf("boyslag.pdf",width=8)
par(mfrow=c(2,2))

plot(c(boys.lagresid[,1],boys.lagresid[,2],boys.lagresid[,3]),
     c(boys.lagresid[,2],boys.lagresid[,3],boys.lagresid[,4]),xlab="at age",
     ylab="at age + 2 years",xlim=c(-3,3),ylim=c(-3,3))
title("Lag 1 (2 years)")

plot(c(boys.lagresid[,1],boys.lagresid[,2]),c(boys.lagresid[,3],boys.lagresid[,4]),
     xlab="at age",ylab="at age + 4 years",xlim=c(-3,3),ylim=c(-3,3))
title("Lag 2 (4 years)")

plot(boys.lagresid[,1],boys.lagresid[,4],xlab="at age",ylab="at age + 6 years",
     xlim=c(-3,3),ylim=c(-3,3))
title("Lag 3 (6 years)")
graphics.off()


