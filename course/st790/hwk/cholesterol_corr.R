#########################################################################
#
#   Cholesterol data - high-dose chenodiol vs. placebo
#
#   1 = high-dose chenodiol, 2 = placebo   
#
#########################################################################

#  Read in the data -- they are already in the long form of one
#  observation per record


thedat.wide <- read.table("cholesterol.dat",row.names=NULL)
colnames(thedat.wide) <- c("trt","id","month0","month6","month12","month20","month24")

#  Total number of individuals

m <- nrow(thedat.wide)

#########################################################################
#
#   Make separate plots by treatment group
#
#########################################################################

#  Extract the three groups and get summary statistics

onedat <- thedat.wide[thedat.wide$trt==1,]
twodat <- thedat.wide[thedat.wide$trt==2,]

m.one <- nrow(onedat)
m.two <- nrow(twodat)

one.mean <- apply(onedat[,3:7],2,mean)
two.mean <- apply(twodat[,3:7],2,mean)

## > one.mean
##   month0   month6  month12  month20  month24 
## 226.7778 249.6111 252.6111 253.1389 256.7222 
## > two.mean
##   month0   month6  month12  month20  month24 
## 236.6452 243.3226 244.5484 261.9032 257.4839 

one.cov <- cov(onedat[,3:7])
two.cov <- cov(twodat[,3:7])

## > one.cov
##           month0   month6   month12   month20   month24
## month0  1962.463 1302.197 1150.8540  952.3460 1009.2794
## month6  1302.197 1715.216 1109.2159 1023.4270 1199.3746
## month12 1150.854 1109.216 1553.9016  696.8556 1265.5746
## month20  952.346 1023.427  696.8556 1147.6087  866.6111
## month24 1009.279 1199.375 1265.5746  866.6111 2545.6921
## > two.cov
##           month0   month6  month12  month20  month24
## month0  3080.437 2342.718 2158.734 2404.831 2086.777
## month6  2342.718 2755.492 2261.284 2392.099 2123.539
## month12 2158.734 2261.284 2267.723 2184.922 1828.959
## month20 2404.831 2392.099 2184.922 2666.957 2012.982
## month24 2086.777 2123.539 1828.959 2012.982 2439.191

sqrt(diag(one.cov))
sqrt(diag(two.cov))

## > sqrt(diag(one.cov))
##   month0   month6  month12  month20  month24 
## 44.29970 41.41516 39.41956 33.87637 50.45485 
## > sqrt(diag(two.cov))
##   month0   month6  month12  month20  month24 
## 55.50168 52.49278 47.62061 51.64259 49.38817 

one.cor <- cov2cor(one.cov)
two.cor <- cov2cor(two.cov)

## > one.cor
##            month0    month6   month12   month20   month24
## month0  1.0000000 0.7097680 0.6590338 0.6345956 0.4515519
## month6  0.7097680 1.0000000 0.6794303 0.7294584 0.5739744
## month12 0.6590338 0.6794303 1.0000000 0.5218361 0.6363163
## month20 0.6345956 0.7294584 0.5218361 1.0000000 0.5070192
## month24 0.4515519 0.5739744 0.6363163 0.5070192 1.0000000
## > two.cor
##            month0    month6   month12   month20   month24
## month0  1.0000000 0.8041079 0.8167669 0.8390164 0.7612846
## month6  0.8041079 1.0000000 0.9046082 0.8824122 0.8191013
## month12 0.8167669 0.9046082 1.0000000 0.8884498 0.7776534
## month20 0.8390164 0.8824122 0.8884498 1.0000000 0.7892396
## month24 0.7612846 0.8191013 0.7776534 0.7892396 1.0000000

#  Scatterplot matrices

one.centscale <- scale(onedat[,3:7],center=one.mean,scale=sqrt(diag(one.cov)))
two.centscale <- scale(twodat[,3:7],center=two.mean,scale=sqrt(diag(two.cov)))

monthlab <- c("Month 0","Month 6","Month 12","Month 20","Month 24")

pdf("trt1.scatter.pdf",width=8)
pairs(one.centscale,label=monthlab,oma=c(5,5,5,5),main="Trt 1")
graphics.off()

pdf("trt2.scatter.pdf",width=8)
pairs(two.centscale,label=monthlab,oma=c(5,5,5,5),main="Trt 2")
graphics.off()

#  Pooled covariance/correlation

cov.pooled <- ( (m.one-1)*one.cov + (m.two-1)*two.cov )/(m-2)
corr.pooled <- cov2cor(cov.pooled)

## > cov.pooled
##           month0   month6  month12  month20  month24
## month0  2478.451 1782.437 1616.030 1622.724 1506.586
## month6  1782.437 2195.344 1640.940 1655.122 1625.912
## month12 1616.030 1640.940 1883.357 1383.655 1525.598
## month20 1622.724 1655.122 1383.655 1848.846 1395.705
## month24 1506.586 1625.912 1525.598 1395.705 2496.538
## > corr.pooled
##            month0    month6   month12   month20   month24
## month0  1.0000000 0.7641399 0.7479848 0.7580608 0.6056681
## month6  0.7641399 1.0000000 0.8070027 0.8215392 0.6945071
## month12 0.7479848 0.8070027 1.0000000 0.7414999 0.7035658
## month20 0.7580608 0.8215392 0.7414999 1.0000000 0.6496422
## month24 0.6056681 0.6945071 0.7035658 0.6496422 1.0000000

#  Autocorrelation function and lag plots -- could write a generic function
#  for length n but too lazy

ac.func <- function(centscale){
   ac.lag1 <- cor(c(centscale[,1],centscale[,2],centscale[,3],centscale[,4]),
        c(centscale[,2],centscale[,3],centscale[,4],centscale[,5]))
   ac.lag2 <- cor(c(centscale[,1],centscale[,2],centscale[,3]),
        c(centscale[,3],centscale[,4],centscale[,5]))
   ac.lag3 <- cor(c(centscale[,1],centscale[,2]),c(centscale[,4],centscale[,5]))
   ac.lag4 <- cor(centscale[,1],centscale[,5])
   ac <- c(ac.lag1,ac.lag2,ac.lag3,ac.lag4)
   ac
}

##  The time points are NOT EQUALLY SPACED; thus, calculation of the
##  autocorrelation function is not really relevant, either at the
##  population or individual levels. Accordingly, the following
##  calculations are not really appropriate here.

ac.one <- ac.func(one.centscale)
ac.two <- ac.func(two.centscale)

## > ac.one
## [1] 0.6045134 0.6749362 0.6042850 0.4515519
## > ac.two
## [1] 0.8466014 0.8256108 0.8290589 0.7612846

lag.plot <- function(plotname,plottitle,centscale){
    pdf(plotname,width=8)
    par(mfrow=c(2,2), oma=c(0,0,2,0))

    plot(c(centscale[,1],centscale[,2],centscale[,3],centscale[,4]),
        c(centscale[,2],centscale[,3],centscale[,4],centscale[,5]), 
        xlab="at week",ylab="at week + 1",xlim=c(-2.5,2.5),ylim=c(-2.5,2.5))
    title("Lag 1")

   plot(c(centscale[,1],centscale[,2],centscale[,3]),
        c(centscale[,3],centscale[,4],centscale[,5]),xlab="at month",
        ylab="at month + 2",xlim=c(-2,2),ylim=c(-2,2))
    title("Lag 2")

    plot(c(centscale[,1],centscale[,2]),c(centscale[,4],centscale[,5]),
        xlab="at week",ylab="at month + 3",xlim=c(-2,2),ylim=c(-2,2))
    title("Lag 3")
    
    plot(centscale[,1],centscale[,5],xlab="at month",
        ylab="at month + 4",xlim=c(-2,2),ylim=c(-2,2))
    title("Lag 4")

    title(plottitle,outer=TRUE)
    graphics.off()
}
    
lag.plot("trt1.lag.pdf","Trt 1",one.centscale)
lag.plot("trt2.lag.pdf","Trt 2",two.centscale)

#  Within-individual correlation

#  Reconfigure as one observation per record ("long" format) and sort
#  in id order

thedat.long <- reshape(thedat.wide,varying=c("month0","month6","month12","month20","month24"),v.names="chol",idvar="id",times=c(0,6,12,20,24),timevar="month",direction="long")

thedat.long <- thedat.long[order(thedat.long$id),]

within.reg <- function(regdat,residplotname,lagplotname,plottitle){

    ids <- unique(regdat[,2])
    m <- length(ids)
    reg.list <- as.list(1:m)
    for (i in 1:m){
          reg.list[[i]] <- regdat[regdat$id==ids[i],]
    }  

    reg.fits <- lapply(reg.list,FUN=function(u){lm(chol ~ month,data=u)})
    reg.residuals <- sapply(reg.fits,residuals)
    reg.pred <- sapply(reg.fits,fitted.values)

    reg.var <- sum(reg.residuals^2)/(m*5-m*2)
    reg.stdresid <- reg.residuals/sqrt(reg.var)

    pdf(residplotname,width=8)
       matplot(reg.pred,reg.stdresid,xlab="predicted values",
        ylab="standardized residual",ylim=c(-3,3),pch=rep(1,m))
      abline(h=0)
      title(paste("Pooled, Standardized Least Squares Residuals,",plottitle))
    graphics.off()

#  Calculate the estimated autocorrelation function

    reg.lagresid <- t(reg.stdresid)
    acf <- ac.func(reg.lagresid)

#  Plot the pooled lagged residuals

    lag.plot(lagplotname,plottitle,reg.lagresid)

    acf
}       

one.regdat <- thedat.long[thedat.long$trt==1,]
two.regdat <- thedat.long[thedat.long$trt==2,]

ac.within.one <-  within.reg(one.regdat,"trt1.resid.pdf","trt1.withinlag.pdf","Trt 1")
ac.within.two <- within.reg(two.regdat,"trt2.resid.pdf","trt2.withinlag.pdf","Trt 2")

## > ac.within.one
## [1] -0.52531241 -0.16856328 -0.04672705  0.39680658
## > ac.within.two
## [1] -0.42004984 -0.37501296 -0.04677818  0.49305693


