#########################################################################
#
#   Weightloss data - three weightloss programs
#
#   1 = Control
#   2 = Diet + Exercise
#   3 = Diet Alone
#
#########################################################################

#  Read in the data -- they are already in the long form of one
#  observation per record

thedat.wide <- read.table("weightloss.dat",row.names=NULL)
colnames(thedat.wide) <-  c("id","month0","month1","month2","month3","month4","program")

#  Total number of individuals

m <- nrow(thedat.wide)

#########################################################################
#
#   Make separate plots by treatment group
#
#########################################################################

#  Extract the three groups and get summary statistics

onedat <- thedat.wide[thedat.wide$program==1,]
twodat <- thedat.wide[thedat.wide$program==2,]
threedat <- thedat.wide[thedat.wide$program==3,]

m.one <- nrow(onedat)
m.two <- nrow(twodat)
m.three <- nrow(threedat)

one.mean <- apply(onedat[,2:6],2,mean)
two.mean <- apply(twodat[,2:6],2,mean)
three.mean <- apply(threedat[,2:6],2,mean)

## > one.mean
##   month0   month1   month2   month3   month4 
## 248.4235 233.4824 240.5324 249.7176 248.5471 
## > two.mean
##   month0   month1   month2   month3   month4 
## 250.4714 216.3107 200.4821 167.0607 163.5286 
## > three.mean
##   month0   month1   month2   month3   month4 
## 246.8974 217.1000 218.0763 218.5684 209.0553

one.cov <- cov(onedat[,2:6])
two.cov <- cov(twodat[,2:6])
three.cov <- cov(threedat[,2:6])

## > one.cov
##          month0   month1   month2   month3   month4
## month0 288.2443 269.2783 219.6750 269.4853 238.7822
## month1 269.2783 383.8609 277.1973 337.8621 231.8051
## month2 219.6750 277.1973 324.0344 280.1430 247.9312
## month3 269.4853 337.8621 280.1430 435.2573 258.7870
## month4 238.7822 231.8051 247.9312 258.7870 334.9729
## > two.cov
##          month0   month1   month2   month3   month4
## month0 502.4806 372.5192 473.4102 460.9755 529.1542
## month1 372.5192 457.3358 405.2617 356.4049 438.4404
## month2 473.4102 405.2617 633.3052 478.0889 529.9320
## month3 460.9755 356.4049 478.0889 536.9684 493.2786
## month4 529.1542 438.4404 529.9320 493.2786 701.6266
## > three.cov
##          month0   month1   month2   month3   month4
## month0 437.9992 347.6230 422.3480 388.0267 384.8047
## month1 347.6230 413.9070 402.5338 354.3354 342.7157
## month2 422.3480 402.5338 540.4089 401.9533 426.5246
## month3 388.0267 354.3354 401.9533 497.1390 402.5894
## month4 384.8047 342.7157 426.5246 402.5894 519.2906

sqrt(diag(one.cov))
sqrt(diag(two.cov))
sqrt(diag(three.cov))

## > sqrt(diag(one.cov))
##   month0   month1   month2   month3   month4 
## 16.97776 19.59237 18.00095 20.86282 18.30226 
## > sqrt(diag(two.cov))
##   month0   month1   month2   month3   month4 
## 22.41608 21.38541 25.16556 23.17258 26.48823 
## > sqrt(diag(three.cov))
##   month0   month1   month2   month3   month4 
## 20.92843 20.34471 23.24670 22.29661 22.78795 

one.cor <- cov2cor(one.cov)
two.cor <- cov2cor(two.cov)
three.cor <- cov2cor(three.cov)

## > one.cor
##           month0    month1    month2    month3    month4
## month0 1.0000000 0.8095321 0.7187944 0.7608197 0.7684520
## month1 0.8095321 1.0000000 0.7859709 0.8265699 0.6464445
## month2 0.7187944 0.7859709 1.0000000 0.7459527 0.7525420
## month3 0.7608197 0.8265699 0.7459527 1.0000000 0.6777424
## month4 0.7684520 0.6464445 0.7525420 0.6777424 1.0000000
## > two.cor
##           month0    month1    month2    month3    month4
## month0 1.0000000 0.7770901 0.8392115 0.8874498 0.8911885
## month1 0.7770901 1.0000000 0.7530285 0.7192033 0.7739982
## month2 0.8392115 0.7530285 1.0000000 0.8198375 0.7949881
## month3 0.8874498 0.7192033 0.8198375 1.0000000 0.8036462
## month4 0.8911885 0.7739982 0.7949881 0.8036462 1.0000000
## > three.cor
##           month0    month1    month2    month3    month4
## month0 1.0000000 0.8164328 0.8681057 0.8315455 0.8068606
## month1 0.8164328 1.0000000 0.8511178 0.7811317 0.7392261
## month2 0.8681057 0.8511178 1.0000000 0.7754886 0.8051515
## month3 0.8315455 0.7811317 0.7754886 1.0000000 0.7923520
## month4 0.8068606 0.7392261 0.8051515 0.7923520 1.0000000

#  Scatterplot matrices

one.centscale <- scale(onedat[,2:6],center=one.mean,scale=sqrt(diag(one.cov)))
two.centscale <- scale(twodat[,2:6],center=two.mean,scale=sqrt(diag(two.cov)))
three.centscale <- scale(threedat[,2:6],center=three.mean,scale=sqrt(diag(three.cov)))

monthlab <- c("Month 0","Month 3","Month 6","Month 9","Month 12")

pdf("program1.scatter.pdf",width=8)
pairs(one.centscale,label=monthlab,oma=c(5,5,5,5),main="Program 1")
graphics.off()

pdf("program2.scatter.pdf",width=8)
pairs(two.centscale,label=monthlab,oma=c(5,5,5,5),main="Program 2")
graphics.off()

pdf("program3.scatter.pdf",width=8)
pairs(three.centscale,label=monthlab,oma=c(5,5,5,5),main="Program 3")
graphics.off()

#  Pooled covariance/correlation

cov.pooled <- ( (m.one-1)*one.cov + (m.two-1)*two.cov + (m.three-1)*three.cov )/(m-3)
corr.pooled <- cov2cor(cov.pooled)

## > cov.pooled
##          month0   month1   month2   month3   month4
## month0 405.0001 327.8995 367.6106 368.0035 375.3067
## month1 327.8995 415.7736 360.6528 349.3071 331.6282
## month2 367.6106 360.6528 492.6547 381.7051 394.5495
## month3 368.0035 349.3071 381.7051 487.1730 378.9103
## month4 375.3067 331.6282 394.5495 378.9103 507.3379
## > corr.pooled
##           month0    month1    month2    month3    month4
## month0 1.0000000 0.7990699 0.8229798 0.8284815 0.8279615
## month1 0.7990699 1.0000000 0.7968739 0.7761353 0.7220619
## month2 0.8229798 0.7968739 1.0000000 0.7791392 0.7891900
## month3 0.8284815 0.7761353 0.7791392 1.0000000 0.7621600
## month4 0.8279615 0.7220619 0.7891900 0.7621600 1.0000000

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

ac.one <- ac.func(one.centscale)
ac.two <- ac.func(two.centscale)
ac.three <- ac.func(three.centscale)

## > ac.one
## [1] 0.7547995 0.7659688 0.7036321 0.7684520
## > ac.two
## [1] 0.7884006 0.7844676 0.8307240 0.8911885
## > ac.three
## [1] 0.8088478 0.8181296 0.7853858 0.8068606


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
    
lag.plot("program1.lag.pdf","Program 1",one.centscale)
lag.plot("program2.lag.pdf","Program 2",two.centscale)
lag.plot("program3.lag.pdf","Program 3",three.centscale)

#  Within-individual correlation

#  Reconfigure as one observation per record ("long" format) and sort
#  in id order

thedat.long <- reshape(thedat.wide,varying=c("month0","month1","month2","month3","month4"),
            v.names="weight",idvar="id",times=c(0,3,6,9,12),timevar="month",direction="long")

thedat.long <- thedat.long[order(thedat.long$id),]

within.reg <- function(regdat,residplotname,lagplotname,plottitle){

    ids <- unique(regdat[,1])
    m <- length(ids)
    reg.list <- as.list(1:m)
    for (i in 1:m){
          reg.list[[i]] <- regdat[regdat$id==ids[i],]
    }  

    reg.fits <- lapply(reg.list,FUN=function(u){lm(weight ~ month,data=u)})
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

#  Return autocorrelation function and pooled within-subject variance estimate
    
    return(list(acf=acf,reg.var=reg.var))
}       

one.regdat <- thedat.long[thedat.long$program==1,]
two.regdat <- thedat.long[thedat.long$program==2,]
three.regdat <- thedat.long[thedat.long$program==3,]

ac.within.one <-  within.reg(one.regdat,"program1.resid.pdf","program1.withinlag.pdf","Program 1")
ac.within.two <- within.reg(two.regdat,"program2.resid.pdf","program2.withinlag.pdf","Program 2")
ac.within.three <- within.reg(three.regdat,"program3.resid.pdf","program3.withinlag.pdf","Program 3")

## > ac.within.one$acf
## [1] -0.45783236 -0.31617228 -0.05973742  0.71891553
## > ac.within.one$reg.var
## [1] 145.8541
## > ac.within.two$acf
## [1] -0.58128590 -0.02545299 -0.38791865  0.43458802
## > ac.within.two$reg.var
## [1] 204.2702
## > ac.within.three$acf
## [1] -0.4538525 -0.4253651  0.1069522  0.4866333
## > ac.within.three$reg.var
## [1] 185.0339

